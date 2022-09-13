from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs
)
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
import pandas as pd
import numpy as np
import nelpy as nel
import os
from sklearn.linear_model import LinearRegression
import glob
import pickle
import logging
import random


def get_cross_region_assemblies(
    basepath,
    regions=None,  # brain regions to load
    deepSuperficial=None,
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,
    epoch=None,
):
    m1 = assembly_reactivation.AssemblyReact(
        basepath,
        brainRegion=regions[0] + "|" + regions[1],
        putativeCellType=putativeCellType,
        weight_dt=weight_dt,
        z_mat_dt=z_mat_dt,
    )
    st, cell_metrics = loading.load_spikes(
        basepath,
        brainRegion=regions[0] + "|" + regions[1],
        putativeCellType=putativeCellType,
    )
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)
    m1.st = st.iloc[
        :,
        cell_metrics.deepSuperficial.str.contains(deepSuperficial + "|unknown").values,
    ]
    m1.cell_metrics = cell_metrics[
        cell_metrics.deepSuperficial.str.contains(deepSuperficial + "|unknown").values
    ]
    # check for any cortex cells
    if (
        not m1.cell_metrics.brainRegion.str.contains(regions[0]).any()
        & m1.cell_metrics.brainRegion.str.contains(regions[1]).any()
    ):
        return None, None

    m1.get_weights(epoch=epoch)
    if len(m1.patterns) == 0:
        return None, None

    return m1

def locate_task_epoch(epoch_df):
    epoch_df = epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    epoch_df.sort_values("duration", ascending=False, inplace=True)
    return int(epoch_df[~epoch_df.environment.str.contains("sleep")].index[0])

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
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    cross_regions=(("CA1", "PFC"), ("CA1", "EC1|EC2|EC3|EC4|EC5|MEC")),
    allowed_pairs=[("CA1", "MEC"), ("CA1", "PFC")],
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=60,
    shuffles=500
):

    epoch_df = loading.load_epoch(basepath)
    epochs = nel.EpochArray(np.array([epoch_df.startTime, epoch_df.stopTime]).T)

    # load in spike data
    st, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=regions
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)
    # simplify brain region names
    cm.loc[cm.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"), "brainRegion"] = "MEC"
    cm.loc[cm.brainRegion.str.contains("CA1"), "brainRegion"] = "CA1"

    if st.isempty:
        return None
    if st.n_active < 2:
        return None
        
    # find longest non-sleep task
    task_idx = locate_task_epoch(epoch_df)

    # locate pairs of brain regions to include in the analysis
    pairs = find_valid_pairs(st, cm, allowed_pairs)
    pairs = remove_middle(cm, pairs)

    # get cross region assemblies
    assem_labels = []
    assembly_act = []
    abscissa_vals = []
    m1 = {}
    for cross_region in cross_regions:
        for deepSuperficial in ["Deep", "Superficial"]:
            m1_ = get_cross_region_assemblies(
                basepath,
                regions=cross_region,
                deepSuperficial=deepSuperficial,
                putativeCellType=putativeCellType,
                weight_dt=weight_dt,
                z_mat_dt=z_mat_dt,
                epoch=epochs[task_idx],
            )
            if m1 is None:
                continue

            assem_labels.append(
                [deepSuperficial + "_" + cross_region[1]] * assembly_act_.n_signals
            )
            m1[deepSuperficial + "_" + cross_region[1]] = m1_

            patterns_keep, is_member_keep, keep_assembly, is_member = find_sig_assembly.main(m1.patterns)

    if len(assembly_act) == 0:
        return None


    def load_results(save_path,verbose=False):
    """
    load_results: load results from a pickle file
    """
    print("Loading results...")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    label_df = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        label_df = pd.concat([label_df, results["label_df"]], ignore_index=True)
    return label_df