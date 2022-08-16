import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs,
)
import pandas as pd
import numpy as np
import nelpy as nel
import os
from itertools import combinations
from ripple_heterogeneity.place_cells import maps


def get_epochs(basepath):
    # load session epoch data
    epoch_df = loading.load_epoch(basepath)

    # compress repeated sleep sessions
    epoch_df = compress_repeated_epochs.main(epoch_df)
    # search for pre task post epochs
    idx, _ = functions.find_pre_task_post(epoch_df.environment)
    if idx is None:
        return None, None

    epoch_df = epoch_df[idx]
    epochs = nel.EpochArray(
        [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
        label=epoch_df.environment.values,
    )
    return epochs, epoch_df


def find_valid_pairs(st, cm, allowed_pairs):
    # locate pairs of brain regions to include in the analysis
    x = np.arange(0, st.data.shape[0])
    pairs = np.array(list(combinations(x, 2)))

    # remove pairs that compare the same brain region
    pairs = pairs[
        cm.brainRegion.iloc[pairs[:, 0]].values
        != cm.brainRegion.iloc[pairs[:, 1]].values
    ]

    # refine to only include allowed pairs
    keep_idx = []
    for pair in allowed_pairs:
        keep_idx.append(
            (cm.brainRegion.iloc[pairs[:, 0]].values == pair[0])
            & (cm.brainRegion.iloc[pairs[:, 1]].values == pair[1])
        )
        keep_idx.append(
            (cm.brainRegion.iloc[pairs[:, 0]].values == pair[1])
            & (cm.brainRegion.iloc[pairs[:, 1]].values == pair[0])
        )
    pairs = pairs[np.vstack(keep_idx).any(axis=0)]
    return pairs


def remove_middle(cm, pairs):
    pairs = pairs[cm.deepSuperficial.iloc[pairs[:, 0]] != "middle"]
    pairs = pairs[cm.deepSuperficial.iloc[pairs[:, 1]] != "middle"]
    return pairs


def run(
    basepath,
    putativeCellType="Pyr",  # type of cell to use for the analysis
    brainRegions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to include
    allowed_pairs=[("CA1", "MEC"), ("CA1", "PFC")],  # allowed pairs of brain regions
    rip_exp_start=0.1,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.1,  # ripple expansion stop, in seconds, how much to expand ripples
    binsize=0.005,
    nbins=200,
    bst_ds_sleep=0.05,
    bst_ds_wake=0.125,
):

    epochs, epoch_df = get_epochs(basepath)
    if epochs is None:
        return None

    # load in spike data
    st, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegions
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)
    # simplify brain region names
    cm.loc[
        cm.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"), "brainRegion"
    ] = "MEC"
    cm.loc[cm.brainRegion.str.contains("CA1"), "brainRegion"] = "CA1"

    if st.isempty:
        return None
    if st.n_active < 2:
        return None

    # load ripples and apply ripple expansion
    ripples_df = loading.load_ripples_events(basepath)
    ripples = (
        nel.EpochArray(np.array([ripples_df.peaks, ripples_df.peaks]).T)
        .expand(rip_exp_start, direction="start")
        .expand(rip_exp_stop, direction="stop")
    )

    # load states to restrict analysis
    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict["NREMstate"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    position_df = loading.load_animal_behavior(basepath)
    if position_df.x.isna().all():
        return None
    bad_idx = np.isnan(position_df.x)
    position_df = position_df[~bad_idx]
    
    pos = nel.AnalogSignalArray(
        data=np.array(position_df[["x", "y"]].values.T),
        timestamps=position_df.time.values,
        step=position_df.sr.iloc[0],
        fs=position_df.sr.iloc[0],
    )

    # check if there is enough data to do the analysis
    try:
        for ep, env_label, ep_label in zip(
            epochs, epoch_df.environment.values, ["pre", "task", "post"]
        ):
        # we need spiking data for each epoch within nrem and wake
            if (ep_label == "pre") | (ep_label == "post"):
                if st[ripples][nrem_epochs][ep].isempty:
                    return
            else:
                # for wake, we also need position data
                if st[ripples][wake_epochs][ep].isempty | pos[epochs[1]].isempty:
                    return
    except:
        return

    # locate pairs of brain regions to include in the analysis
    pairs = find_valid_pairs(st, cm, allowed_pairs)
    pairs = remove_middle(cm, pairs)

    # set up vars for results
    ccgs = pd.DataFrame()
    label_df = pd.DataFrame()

    for ep, env_label, ep_label in zip(
        epochs, epoch_df.environment.values, ["pre", "task", "post"]
    ):
        # if epoch is pre or post sleep, actually restrict to NREM
        if (ep_label == "pre") | (ep_label == "post"):
            restriction = ripples[nrem_epochs][ep]
            current_st = st[restriction]
            ccgs_temp = functions.pairwise_cross_corr(
                current_st.data,
                return_index=False,
                binsize=binsize,
                nbins=nbins,
                pairs=pairs,
            )
            bst = current_st.bin(ds=bst_ds_sleep).data
            corr_pearson, pval_pearson, _ = functions.pairwise_corr(
                bst, method="pearson", pairs=pairs
            )
            corr_spearman, pval_spearman, _ = functions.pairwise_corr(
                bst, method="spearman", pairs=pairs
            )
            spatial_corr = np.zeros_like(corr_pearson) * np.nan
            spatial_info_ref = np.zeros_like(corr_pearson) * np.nan
            spatial_info_tar = np.zeros_like(corr_pearson) * np.nan
            peak_rate_ref = np.zeros_like(corr_pearson) * np.nan
            peak_rate_tar = np.zeros_like(corr_pearson) * np.nan
            n_spikes_ref = current_st.n_spikes[pairs[:, 0]]
            n_spikes_tar = current_st.n_spikes[pairs[:, 1]]
        else:
            restriction = wake_epochs[ep]
            current_st = st[restriction]

            ccgs_temp = functions.pairwise_cross_corr(
                current_st.data,
                return_index=False,
                binsize=binsize,
                nbins=nbins,
                pairs=pairs,
            )
            bst = current_st.bin(ds=bst_ds_wake).data
            corr_pearson, pval_pearson, _ = functions.pairwise_corr(
                bst, method="pearson", pairs=pairs
            )
            corr_spearman, pval_spearman, _ = functions.pairwise_corr(
                bst, method="spearman", pairs=pairs
            )
            # do pairwise spatial correlation of ratemaps
            spatial_maps = maps.SpatialMap(pos[epochs[1]], current_st, dim=2)
            spatial_corr = functions.pairwise_spatial_corr(
                spatial_maps.tc.ratemap, pairs=pairs
            )

            spatial_info = spatial_maps.tc.spatial_information()
            spatial_info_ref = spatial_info[pairs[:, 0]]
            spatial_info_tar = spatial_info[pairs[:, 1]]

            peak_rate = [np.nanmax(x) for x in spatial_maps.tc.ratemap]
            peak_rate = np.array(peak_rate)
            peak_rate_ref = peak_rate[pairs[:, 0]]
            peak_rate_tar = peak_rate[pairs[:, 1]]

            n_spikes_ref = current_st.n_spikes[pairs[:, 0]]
            n_spikes_tar = current_st.n_spikes[pairs[:, 1]]

        ccgs = pd.concat([ccgs, ccgs_temp], axis=1, ignore_index=True)

        label_df_temp = pd.DataFrame()
        label_df_temp["corr_pearson"] = corr_pearson
        label_df_temp["pval_pearson"] = pval_pearson
        label_df_temp["corr_spearman"] = corr_spearman
        label_df_temp["pval_spearman"] = pval_spearman
        label_df_temp["spatial_corr"] = spatial_corr
        label_df_temp["spatial_info_ref"] = spatial_info_ref
        label_df_temp["spatial_info_tar"] = spatial_info_tar
        label_df_temp["peak_rate_ref"] = peak_rate_ref
        label_df_temp["peak_rate_tar"] = peak_rate_tar
        label_df_temp["n_spikes_ref"] = n_spikes_ref
        label_df_temp["n_spikes_tar"] = n_spikes_tar
        label_df_temp["deepSuperficial_ref"] = cm.iloc[
            pairs[:, 0]
        ].deepSuperficial.values
        label_df_temp["deepSuperficial_target"] = cm.iloc[
            pairs[:, 1]
        ].deepSuperficial.values
        label_df_temp["ref_UID"] = cm.iloc[pairs[:, 0]].UID.values
        label_df_temp["target_UID"] = cm.iloc[pairs[:, 1]].UID.values
        label_df_temp["ref_region"] = cm.iloc[pairs[:, 0]].brainRegion.values
        label_df_temp["target_region"] = cm.iloc[pairs[:, 1]].brainRegion.values
        label_df_temp["ref_id"] = pairs[:, 0]
        label_df_temp["target_id"] = pairs[:, 1]
        label_df_temp["epoch"] = ep_label
        label_df_temp["environment"] = env_label
        label_df_temp["basepath"] = basepath

        label_df = pd.concat([label_df, label_df_temp], ignore_index=True)

    results = {"label_df": label_df, "ccgs": ccgs}
    return results


def load_results(save_path, verbose=False):
    """
    load_results: load results from a pickle file
    """

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    label_df = pd.DataFrame()
    ccgs = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results_ = pickle.load(f)
        if results_ is None:
            continue
        label_df = pd.concat([label_df, results_["label_df"]], ignore_index=True)
        ccgs = pd.concat([ccgs, results_["ccgs"]], axis=1, ignore_index=True)

    return label_df, ccgs
