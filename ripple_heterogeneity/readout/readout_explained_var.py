import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel
from quantities import s
import quantities as pq
from neo.core import SpikeTrain
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import correlation_coefficient
from ripple_heterogeneity.utils import compress_repeated_epochs
import itertools


def get_corrcoef(st, epoch, bin_size=0.05):
    """
    Calculate correlation coefficient for epoch
    input:
        st: nel.SpikeTrain object
        epoch: nel.Epoch object
        bin_size: bin size in seconds
    output:
        corrcoef_r: correlation matrix
    """
    spk_train_list = st[epoch]
    if spk_train_list.isempty:
        return None
    bst = spk_train_list.bin(ds=bin_size)

    return np.corrcoef(bst.data)


def get_cells(
    basepath, ref="CA1", target="PFC", ref_sublayer="Deep", putativeCellType="Pyr"
):
    """
    Get cells from a specific region and sublayer (Deep or Superficial)
    input:
        basepath: path to ripple_heterogeneity folder
        ref: reference region (e.g. "CA1")
        target: target region (e.g. "PFC")
        ref_sublayer: reference sublayer (e.g. "Deep")
    output:
        st: neo.SpikeTrain object
        cell_metrics: pandas dataframe with cell metrics
    """
    # get cells from these regions
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion=[ref, target], putativeCellType=putativeCellType
    )
    # classify deep and superficial cells by distance to pyr layer
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)
    # re-label any ca1 to ca1
    cell_metrics.loc[cell_metrics.brainRegion.str.contains(ref), "brainRegion"] = ref
    # restrict to deep superficial
    idx = (
        cell_metrics.deepSuperficial.str.contains(ref_sublayer)
        | cell_metrics.brainRegion.str.contains(target)
    ).values
    cell_metrics = cell_metrics[idx]
    st = st[:, np.where(idx)[0].astype(int) + 1]
    return st, cell_metrics


def remov_within_reg_corr(U, brainRegion):
    """
    Remove correlations within a brain region
    input:
        U: correlation matrix (n x n)
        brainRegion: pandas column, brain region to remove correlations from (e.g. "cell_metrics.brainRegion")
    output:
        U: correlation matrix with correlations within brain region removed
    """
    for i in range(len(brainRegion)):
        for j in range(len(brainRegion)):
            if brainRegion.iloc[i] == brainRegion.iloc[j]:
                U[i, j] = np.nan
    return U


def get_explained_var(
    st,
    beh_epochs,
    cell_metrics,
    state_epoch,
    task_binsize=0.125,
    restrict_task=False,
    theta_epochs=None,
):
    """
    Calculate explained variance
    input:
        st: neo.SpikeTrain object
        beh_epochs: nel.EpochArray object with 3 epochs: sleep, task, sleep
        cell_metrics: pandas dataframe with cell metrics
    output:
        EV: explained variance
        rEV: reverse explained variance
    """

    # get correlation matrix per epoch

    # restrict spike times to epoch (could be many things)
    st_restrict = st[state_epoch]

    # restrict to spike times to theta epochs
    if theta_epochs is not None:
        st = st[theta_epochs]

    # pre task
    corrcoef_r_pre = get_corrcoef(st_restrict, beh_epochs[0], bin_size=0.05)

    # task
    if restrict_task:
        corrcoef_r_beh = get_corrcoef(st_restrict, beh_epochs[1], bin_size=0.05)
    else:
        corrcoef_r_beh = get_corrcoef(st, beh_epochs[1], bin_size=task_binsize)

    # post task
    corrcoef_r_post = get_corrcoef(st_restrict, beh_epochs[2], bin_size=0.05)

    # get uids for ref and target cells
    c = np.array(list(itertools.product(cell_metrics.UID.values, repeat=2)))
    ref_uid = c[:, 0]
    target_uid = c[:, 1]

    if corrcoef_r_pre is None or corrcoef_r_beh is None or corrcoef_r_post is None:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    # remove correlations within region
    corrcoef_r_pre = remov_within_reg_corr(corrcoef_r_pre, cell_metrics.brainRegion)
    corrcoef_r_beh = remov_within_reg_corr(corrcoef_r_beh, cell_metrics.brainRegion)
    corrcoef_r_post = remov_within_reg_corr(corrcoef_r_post, cell_metrics.brainRegion)

    # remove upper triangle correlations
    corrcoef_r_pre[np.tril_indices(corrcoef_r_pre.shape[0], 1)] = np.nan
    corrcoef_r_beh[np.tril_indices(corrcoef_r_beh.shape[0], 1)] = np.nan
    corrcoef_r_post[np.tril_indices(corrcoef_r_post.shape[0], 1)] = np.nan

    # flatten and calculate cross-correlation between epochs
    corr_df = pd.DataFrame(
        {
            "r_pre": corrcoef_r_pre.flatten(),
            "r_beh": corrcoef_r_beh.flatten(),
            "r_post": corrcoef_r_post.flatten(),
        }
    ).corr()
    # pull out specific between epoch correlations
    beh_pos = corr_df.loc["r_beh", "r_post"]
    beh_pre = corr_df.loc["r_beh", "r_pre"]
    pre_pos = corr_df.loc["r_pre", "r_post"]

    # calculate explained variance
    EV = (
        (beh_pos - beh_pre * pre_pos) / np.sqrt((1 - beh_pre**2) * (1 - pre_pos**2))
    ) ** 2
    # calculate reverse explained variance
    rEV = (
        (beh_pre - beh_pos * pre_pos) / np.sqrt((1 - beh_pos**2) * (1 - pre_pos**2))
    ) ** 2

    return (
        EV,
        rEV,
        corrcoef_r_pre.flatten(),
        corrcoef_r_beh.flatten(),
        corrcoef_r_post.flatten(),
        ref_uid,
        target_uid,
    )


def run(
    basepath,  # path to data folder
    reference_region="CA1",  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    putativeCellType="Pyr",  # cell type
    min_cells=5,  # minimum number of cells per region
    restrict_task=False,  # restrict restriction_type to task epochs (ex. ripples in task (True) vs. all task (False))
    restriction_type="ripples",  # "ripples" or "NREMstate" or "barrage"
    ripple_expand=0.05,  # in seconds, how much to expand ripples
    task_binsize=0.125,  # in seconds, bin size for task epochs
    restrict_task_to_theta=True,  # restrict task to theta epochs
):
    # locate epochs
    ep_df = loading.load_epoch(basepath)
    ep_df = compress_repeated_epochs.main(ep_df, epoch_name="sleep")

    # locate pre task post structure
    idx, _ = functions.find_pre_task_post(ep_df.environment)
    if idx is None:
        return None

    ep_df = ep_df[idx]
    beh_epochs = nel.EpochArray(np.array([ep_df.startTime, ep_df.stopTime]).T)

    # choose which times to restrict to
    if restriction_type == "ripples":
        ripples = loading.load_ripples_events(basepath)
        ripple_epochs = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)
        restrict_epochs = ripple_epochs.expand(ripple_expand)
    elif restriction_type == "NREMstate":
        state_dict = loading.load_SleepState_states(basepath)
        restrict_epochs = nel.EpochArray(state_dict["NREMstate"])
    elif restriction_type == "barrage":
        barrage = loading.load_barrage_events(basepath)
        barrage_epochs = nel.EpochArray(np.array([barrage.start, barrage.stop]).T)
        restrict_epochs = barrage_epochs.expand(ripple_expand)
    else:
        raise ValueError("restriction_type must be 'ripples' or 'NREMstate'")

    if restrict_task_to_theta:
        state_dict = loading.load_SleepState_states(basepath)
        theta_epochs = nel.EpochArray(state_dict["THETA"])
    else:
        theta_epochs = None

    # needs exactly 3 epochs for analysis
    if ep_df.shape[0] != 3:
        return None

    evs = []
    revs = []
    sublayers = []
    regions = []
    n_ca1s = []
    n_targets = []
    pairwise_corr = []
    pairwise_corr_epoch = []
    pairwise_corr_region = []
    pairwise_corr_sublayer = []
    pairwise_corr_ref_uid = []
    pairwise_corr_target_uid = []
    for region in target_regions:
        for sublayer in ["Deep", "Superficial"]:
            st, cell_metrics = get_cells(
                basepath,
                ref=reference_region,
                target=region,
                ref_sublayer=sublayer,
                putativeCellType=putativeCellType,
            )
            n_ca1 = cell_metrics.brainRegion.str.contains("CA1").sum()
            n_target = cell_metrics.brainRegion.str.contains(region).sum()
            if st.isempty | (n_ca1 < min_cells) | (n_target < min_cells):
                continue
            (
                ev,
                rev,
                cor_pre,
                cor_beh,
                cor_post,
                ref_uid,
                target_uid,
            ) = get_explained_var(
                st,
                beh_epochs,
                cell_metrics,
                restrict_epochs,
                task_binsize,
                restrict_task,
                theta_epochs,
            )
            evs.append(ev)
            revs.append(rev)
            sublayers.append(sublayer)
            regions.append(region)
            n_ca1s.append(n_ca1)
            n_targets.append(n_target)

            # store pairwise correlations
            pairwise_corr.append(np.hstack([cor_pre, cor_beh, cor_post]))
            pairwise_corr_ref_uid.append(np.hstack([ref_uid, ref_uid, ref_uid]))
            pairwise_corr_target_uid.append(
                np.hstack([target_uid, target_uid, target_uid])
            )
            pairwise_corr_epoch.append(
                np.hstack(
                    [
                        ["pre"] * len(cor_pre),
                        ["task"] * len(cor_beh),
                        ["post"] * len(cor_post),
                    ]
                )
            )
            pairwise_corr_region.append(
                np.hstack(
                    [
                        [region] * len(cor_pre),
                        [region] * len(cor_beh),
                        [region] * len(cor_post),
                    ]
                )
            )
            pairwise_corr_sublayer.append(
                np.hstack(
                    [
                        [sublayer] * len(cor_pre),
                        [sublayer] * len(cor_beh),
                        [sublayer] * len(cor_post),
                    ]
                )
            )

    ev_df = pd.DataFrame()
    ev_df["region"] = regions
    ev_df["sublayer"] = sublayers
    ev_df["ev"] = evs
    ev_df["rev"] = revs
    ev_df["n_ca1"] = n_ca1s
    ev_df["n_target"] = n_targets
    ev_df["basepath"] = basepath

    pairwise_corr_df = pd.DataFrame()
    if len(pairwise_corr) == 0:
        pairwise_corr_df["corr"] = pairwise_corr
        pairwise_corr_df["epoch"] = pairwise_corr_epoch
        pairwise_corr_df["region"] = pairwise_corr_region
        pairwise_corr_df["sublayer"] = pairwise_corr_sublayer
        pairwise_corr_df["ref_uid"] = pairwise_corr_ref_uid
        pairwise_corr_df["target_uid"] = pairwise_corr_target_uid
        pairwise_corr_df["basepath"] = basepath
    else:
        pairwise_corr_df["corr"] = np.hstack(pairwise_corr)
        pairwise_corr_df["epoch"] = np.hstack(pairwise_corr_epoch)
        pairwise_corr_df["region"] = np.hstack(pairwise_corr_region)
        pairwise_corr_df["sublayer"] = np.hstack(pairwise_corr_sublayer)
        pairwise_corr_df["ref_uid"] = np.hstack(pairwise_corr_ref_uid)
        pairwise_corr_df["target_uid"] = np.hstack(pairwise_corr_target_uid)
        pairwise_corr_df["basepath"] = basepath

    results = {"ev_df": ev_df, "pairwise_corr_df": pairwise_corr_df}
    return results


def load_results(save_path, verbose=False):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    ev_df = pd.DataFrame()
    pairwise_corr_df = pd.DataFrame()
    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        ev_df = pd.concat([ev_df, results["ev_df"]], ignore_index=True)
        pairwise_corr_df = pd.concat(
            [pairwise_corr_df, results["pairwise_corr_df"]], ignore_index=True
        )
    return ev_df, pairwise_corr_df
