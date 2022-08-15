import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs
)
import pandas as pd
import numpy as np
import nelpy as nel
import os

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

def run(
    basepath,
    putativeCellType="Pyr",  # type of cell to use for the analysis
    brainRegions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to include
    rip_exp_start=0.5,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.5,  # ripple expansion stop, in seconds, how much to expand ripples
    binsize=0.005,
    nbins=200,
):

    epochs, epoch_df = get_epochs(basepath)
    if epochs is None:
        return None

    # load in spike data
    st, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegions
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

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

    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict["NREMstate"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    ccgs = pd.DataFrame()
    label_df = pd.DataFrame()

    # check if there is enough data to do the analysis
    try:
        for ep, env_label, ep_label in zip(
            epochs, epoch_df.environment.values, ["pre", "task", "post"]
        ):
            if (ep_label == "pre") | (ep_label == "post"):
                if st[ripples][nrem_epochs][ep].isempty:
                    return
            else:
                if st[ripples][wake_epochs][ep].isempty:
                    return
    except:
        return
        
    for ep, env_label, ep_label in zip(
        epochs, epoch_df.environment.values, ["pre", "task", "post"]
    ):
        if (ep_label == "pre") | (ep_label == "post"):
            ccgs_temp, pairs = functions.pairwise_cross_corr(
                st[ripples][nrem_epochs][ep].data, return_index=True, binsize=binsize, nbins=nbins
            )
        else:
            ccgs_temp, pairs = functions.pairwise_cross_corr(
                st[ripples][wake_epochs][ep].data, return_index=True, binsize=binsize, nbins=nbins
            )

        ccgs_temp = pd.concat([ccgs,ccgs_temp],axis=1, ignore_index=True)

        label_df_temp = pd.DataFrame()

        label_df_temp["deepSuperficial_ref"] = cm.iloc[pairs[:, 0]].deepSuperficial.values
        label_df_temp["deepSuperficial_target"] = cm.iloc[pairs[:, 1]].deepSuperficial.values
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
        ccgs = pd.concat([ccgs, results_["ccgs"]],axis=1, ignore_index=True)

    return label_df, ccgs
