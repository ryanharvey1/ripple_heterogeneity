import glob
import itertools
import os
import pickle
from tqdm import tqdm
import numpy as np
import pandas as pd
from ripple_heterogeneity.readout import assembly_multi_region
from ripple_heterogeneity.utils import functions, loading
from ripple_heterogeneity.place_cells import maps
import nelpy as nel
import copy

def get_pairs(curr_assem):
    x = np.arange(0, curr_assem.shape[0])
    pairs = np.array(list(itertools.combinations(x, 2)))

    # add ref and tar metadata
    label_df = pd.DataFrame()
    label_df["idx_ref"] = pairs[:, 0]
    label_df["idx_tar"] = pairs[:, 1]

    label_df["UID_ref"] = curr_assem.UID.iloc[pairs[:, 0]].values
    label_df["UID_tar"] = curr_assem.UID.iloc[pairs[:, 1]].values

    label_df["brainRegion_ref"] = curr_assem.brainRegion.iloc[pairs[:, 0]].values
    label_df["brainRegion_tar"] = curr_assem.brainRegion.iloc[pairs[:, 1]].values

    label_df["deepSuperficial_ref"] = curr_assem.deepSuperficial.iloc[
        pairs[:, 0]
    ].values
    label_df["deepSuperficial_tar"] = curr_assem.deepSuperficial.iloc[
        pairs[:, 1]
    ].values

    label_df["is_member_sig_ref"] = curr_assem.is_member_sig.iloc[pairs[:, 0]].values
    label_df["is_member_sig_tar"] = curr_assem.is_member_sig.iloc[pairs[:, 1]].values

    # relabel as simple keys
    label_df.loc[
        label_df.brainRegion_ref.str.contains("CA1"), "brainRegion_ref"
    ] = "CA1"
    label_df.loc[
        label_df.brainRegion_tar.str.contains("CA1"), "brainRegion_tar"
    ] = "CA1"

    label_df.loc[
        label_df.brainRegion_ref.str.contains("EC5|EC4|EC2|EC3|EC1"), "brainRegion_ref"
    ] = "MEC"
    label_df.loc[
        label_df.brainRegion_tar.str.contains("EC5|EC4|EC2|EC3|EC1"), "brainRegion_tar"
    ] = "MEC"

    # remove within region comparisons
    label_df = label_df.query("brainRegion_ref != brainRegion_tar")

    # make sure comparisons are cross region
    idx = (
        ((label_df.brainRegion_ref == "CA1") & (label_df.brainRegion_tar == "MEC"))
        | ((label_df.brainRegion_ref == "MEC") & (label_df.brainRegion_tar == "CA1"))
    ) | (
        ((label_df.brainRegion_ref == "CA1") & (label_df.brainRegion_tar == "PFC"))
        | ((label_df.brainRegion_ref == "PFC") & (label_df.brainRegion_tar == "CA1"))
    )

    label_df = label_df[idx]

    # put cortex always as target
    idx = label_df.brainRegion_tar.str.contains("CA1")
    idx_tar = label_df.loc[idx, "idx_tar"].values
    label_df.loc[idx, "idx_tar"] = label_df.loc[idx, "idx_ref"]
    label_df.loc[idx, "idx_ref"] = idx_tar

    UID_tar = label_df.loc[idx, "UID_tar"].values
    label_df.loc[idx, "UID_tar"] = label_df.loc[idx, "UID_ref"]
    label_df.loc[idx, "UID_ref"] = UID_tar

    brainRegion_tar = label_df.loc[idx, "brainRegion_tar"].values
    label_df.loc[idx, "brainRegion_tar"] = label_df.loc[idx, "brainRegion_ref"]
    label_df.loc[idx, "brainRegion_ref"] = brainRegion_tar

    deepSuperficial_tar = label_df.loc[idx, "deepSuperficial_tar"].values
    label_df.loc[idx, "deepSuperficial_tar"] = label_df.loc[idx, "deepSuperficial_ref"]
    label_df.loc[idx, "deepSuperficial_ref"] = deepSuperficial_tar

    is_member_sig_tar = label_df.loc[idx, "is_member_sig_tar"]
    label_df.loc[idx, "is_member_sig_tar"] = label_df.loc[idx, "is_member_sig_ref"]
    label_df.loc[idx, "is_member_sig_ref"] = is_member_sig_tar

    # discard middle cells
    label_df = label_df[
        (label_df.brainRegion_ref == "CA1")
        & (label_df.deepSuperficial_ref.str.contains("Deep|Superficial"))
    ]
    return label_df


def run(basepath,binsize=0.005, nbins=200):

    with open(basepath, "rb") as f:
        results = pickle.load(f)
    if results is None:
        return None

    # pull in previous results from assembly multi region analysis
    prop_df, assembly_df, _ = assembly_multi_region.compile_results_df(results)
    m1 = results["react"]

    # if there is only single epoch, that must by task
    if m1.epoch_df.name.shape[0] == 1:
        task_idx = 0
    # if there are not exactly 3 epochs, find the longest task
    elif m1.epoch_df.name.shape[0] != 3:
        epoch_df = m1.epoch_df.reset_index()
        epoch_df = epoch_df.query("environment != 'sleep'")
        epoch_df["duration"] = epoch_df.stopTime.values - epoch_df.startTime.values
        task_idx = int(epoch_df.sort_values("duration", ascending=False).index[0])
    # if there is exactly 3 epochs, the center will be the task
    else:
        task_idx = 1

    # load position
    position_df = loading.load_animal_behavior(results['react'].basepath)
    position_df_no_nan = position_df.query("not x.isnull() & not y.isnull()")

    # if there is no position, skip session
    if position_df_no_nan.shape[0] == 0:
        return None

    # put position into position array
    pos = nel.PositionArray(
        data=position_df_no_nan[["x", "y"]].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )

    if pos[m1.epochs[task_idx]].isempty:
        return None

    # calculate tuning curves
    tc = maps.SpatialMap(pos[m1.epochs[task_idx]], m1.st[m1.epochs[task_idx]], dim=2)
    # tc.shuffle_spatial_information()
    
    label_df = pd.DataFrame()
    ccgs = pd.DataFrame()

    current_st = m1.st[m1.epochs[task_idx]]

    # deep copy ratemap and nan low occ
    ratemap = copy.deepcopy(tc.tc.ratemap)
    ratemap[:,tc.tc.occupancy < 0.1] = np.nan

    for assembly_n in assembly_df.assembly_n.unique():
        curr_assem = assembly_df.query("assembly_n == @assembly_n")

        label_df_ = get_pairs(curr_assem)

        label_df_["assembly_n"] = assembly_n

        spatial_corr = functions.pairwise_spatial_corr(
            ratemap, return_index=False, pairs=label_df_[["idx_ref", "idx_tar"]].values
        )
        label_df_["spatial_corr"] = spatial_corr

        label_df_["spatial_info_ref"] = tc.tc.spatial_information()[label_df_.idx_ref]
        label_df_["spatial_info_tar"] = tc.tc.spatial_information()[label_df_.idx_tar]

        # label_df_["spatial_pval_ref"] = tc.spatial_information_pvalues[label_df_.idx_ref]
        # label_df_["spatial_pval_tar"] = tc.spatial_information_pvalues[label_df_.idx_tar]

        label_df_["spatial_sparsity_ref"] = tc.tc.spatial_sparsity()[label_df_.idx_ref]
        label_df_["spatial_sparsity_tar"] = tc.tc.spatial_sparsity()[label_df_.idx_tar]

        label_df_["peak_rate_ref"] = tc.tc.ratemap.max(axis=1).max(axis=1)[label_df_.idx_ref]
        label_df_["peak_rate_tar"] = tc.tc.ratemap.max(axis=1).max(axis=1)[label_df_.idx_tar]

        label_df_["n_spikes_ref"] = current_st.n_events[label_df_.idx_ref]
        label_df_["n_spikes_tar"] = current_st.n_events[label_df_.idx_tar]
        
        ccgs_ = functions.pairwise_cross_corr(
            current_st.data,
            binsize=binsize,
            nbins=nbins,
            pairs=label_df_[["idx_ref", "idx_tar"]].values,
        )
        ccgs = pd.concat([ccgs, ccgs_], axis=1, ignore_index=True)

        label_df = pd.concat([label_df, label_df_], ignore_index=True)

    label_df["basepath"] = m1.basepath

    results = {"ccgs": ccgs, "label_df": label_df}

    return results


def load_results(save_path, verbose=False):

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    ccgs = pd.DataFrame()
    label_df = pd.DataFrame()

    for session in tqdm(sessions):
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        ccgs = pd.concat([ccgs, results["ccgs"]], axis=1, ignore_index=True)
        label_df = pd.concat([label_df, results["label_df"]], ignore_index=True)

    return ccgs, label_df
