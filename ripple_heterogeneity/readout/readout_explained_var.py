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


def get_corrcoef(st, epoch, bin_size=0.50):
    """
    Calculate correlation coefficient for epoch
    input:
        st: nel.SpikeTrain object
        epoch: nel.Epoch object
        bin_size: bin size in seconds
    output:
        corrcoef_r: correlation matrix
    """
    spiketrain = []
    for spk in st[epoch].data:
        spiketrain.append(
            SpikeTrain(spk, t_start=epoch.start, t_stop=epoch.stop, units="s")
        )
    corrcoef = correlation_coefficient(
        BinnedSpikeTrain(spiketrain, bin_size=bin_size * pq.s)
    )
    return corrcoef


def get_cells(basepath, ref="CA1", target="PFC", ref_sublayer="Deep"):
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
    st, cell_metrics = loading.load_spikes(basepath, brainRegion=[ref, target])
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


def get_explained_var(st, beh_epochs, cell_metrics):
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
    corrcoef_r_pre = get_corrcoef(st, beh_epochs[0])
    corrcoef_r_beh = get_corrcoef(st, beh_epochs[1])
    corrcoef_r_post = get_corrcoef(st, beh_epochs[2])

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

    return EV, rEV


def run(
    basepath, reference_region="CA1", target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"]
):
    # locate epochs
    ep_df = loading.load_epoch(basepath)
    ep_df = compress_repeated_epochs.main(ep_df, epoch_name="sleep")
    # idx = functions.find_epoch_pattern(ep_df.environment, ["sleep", "linear", "sleep"])
    # ep_df = ep_df[idx[0]]
    beh_epochs = nel.EpochArray(np.array([ep_df.startTime, ep_df.stopTime]).T)

    if ep_df.shape[0] != 3:
        return None

    evs = []
    revs = []
    sublayers = []
    regions = []
    n_ca1s = []
    n_targets = []
    for region in target_regions:
        for sublayer in ["Deep", "Superficial"]:
            st, cell_metrics = get_cells(
                basepath, ref=reference_region, target=region, ref_sublayer=sublayer
            )
            n_ca1 = cell_metrics.brainRegion.str.contains("CA1").sum()
            n_target = cell_metrics.brainRegion.str.contains(region).sum()
            if st.isempty | (n_ca1 < 5) | (n_target < 5):
                continue
            ev, rev = get_explained_var(st, beh_epochs, cell_metrics)
            evs.append(ev)
            revs.append(rev)
            sublayers.append(sublayer)
            regions.append(region)
            n_ca1s.append(n_ca1)
            n_targets.append(n_target)

    results = pd.DataFrame()
    results["region"] = regions
    results["sublayer"] = sublayers
    results["ev"] = evs
    results["rev"] = revs
    results["n_ca1"] = n_ca1s
    results["n_target"] = n_targets
    results["basepath"] = basepath

    return results


def load_results(save_path):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    df = pd.DataFrame()
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        df = pd.concat([df, results], ignore_index=True)
    return df
