import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
)
import pandas as pd
import numpy as np
import nelpy as nel
import os


def run(
    basepath,
    putativeCellType="Pyr",  # type of cell to use for the analysis
    brainRegions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to include
    rip_exp_start=0.2,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.2,  # ripple expansion stop, in seconds, how much to expand ripples
    binsize=0.005,
    nbins=200,
):

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

    ccgs, pairs = functions.pairwise_cross_corr(
        st[ripples].data, return_index=True, binsize=binsize, nbins=nbins
    )

    label_df = pd.DataFrame()
    label_df["deepSuperficial_ref"] = cm.iloc[pairs[:, 0]].deepSuperficial.values
    label_df["deepSuperficial_target"] = cm.iloc[pairs[:, 1]].deepSuperficial.values
    label_df["ref_UID"] = cm.iloc[pairs[:, 0]].UID.values
    label_df["target_UID"] = cm.iloc[pairs[:, 1]].UID.values
    label_df["ref_region"] = cm.iloc[pairs[:, 0]].brainRegion.values
    label_df["target_region"] = cm.iloc[pairs[:, 1]].brainRegion.values
    label_df["ref_id"] = pairs[:, 0]
    label_df["target_id"] = pairs[:, 1]
    label_df["basepath"] = basepath

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
