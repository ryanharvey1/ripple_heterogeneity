from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs
)
from ripple_heterogeneity.assembly import assembly_reactivation
import pandas as pd
import numpy as np
import nelpy as nel
import os
from sklearn.linear_model import LinearRegression
import glob
import pickle
import logging
import random

logging.getLogger().setLevel(logging.ERROR)

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

    assembly_act = m1.get_assembly_act(epoch=epoch)
    return m1, assembly_act

def locate_task_epoch(epoch_df):
    epoch_df = epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    epoch_df.sort_values("duration", ascending=False, inplace=True)
    return int(epoch_df[~epoch_df.environment.str.contains("sleep")].index[0])

def find_slope_over_time(assembly_act):
    X = assembly_act.abscissa_vals.reshape(-1, 1)
    y = assembly_act.data.T
    reg = LinearRegression().fit(X, y)
    return reg.coef_.flatten()

def find_slope_over_time_shuffle(assembly_act):
    random_order = random.sample(range(len(assembly_act.abscissa_vals)), len(assembly_act.abscissa_vals))
    X = assembly_act.abscissa_vals[random_order].reshape(-1, 1)
    y = assembly_act.data.T
    reg = LinearRegression().fit(X, y)
    return reg.coef_.flatten()

def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    cross_regions=(("CA1", "PFC"), ("CA1", "EC1|EC2|EC3|EC4|EC5|MEC")),
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=60,
    shuffles=500
):

    epoch_df = loading.load_epoch(basepath)
    epochs = nel.EpochArray(np.array([epoch_df.startTime, epoch_df.stopTime]).T)

    # find longest non-sleep task
    task_idx = locate_task_epoch(epoch_df)

    # get cross region assemblies
    assem_labels = []
    assembly_act = []
    abscissa_vals = []
    m1 = {}
    for cross_region in cross_regions:
        for deepSuperficial in ["Deep", "Superficial"]:
            m1_, assembly_act_ = get_cross_region_assemblies(
                basepath,
                regions=cross_region,
                deepSuperficial=deepSuperficial,
                putativeCellType=putativeCellType,
                weight_dt=weight_dt,
                z_mat_dt=z_mat_dt,
                epoch=epochs[task_idx],
            )
            if assembly_act_ is None:
                continue

            assembly_act.append(assembly_act_.data)
            abscissa_vals.append(assembly_act_.abscissa_vals)
            assem_labels.append(
                [deepSuperficial + "_" + cross_region[1]] * assembly_act_.n_signals
            )
            m1[deepSuperficial + "_" + cross_region[1]] = m1_

    if len(assembly_act) == 0:
        return None

    assembly_act = nel.AnalogSignalArray(
        data=np.vstack(assembly_act),
        timestamps=abscissa_vals[0],
    )
    assem_labels = np.hstack(assem_labels)

    slopes = find_slope_over_time(assembly_act)

    slopes_shuff = []
    for _ in range(shuffles):
        slopes_shuff.append(find_slope_over_time_shuffle(assembly_act))

    sig_event_idx, pvals = functions.get_significant_events(slopes, np.vstack(slopes_shuff), q=95, tail="both")


    label_df = pd.DataFrame()
    label_df["assembly_n"] = np.arange(assembly_act.n_signals).astype(int)
    label_df["slopes"] = slopes
    label_df["pvals"] = pvals
    label_df["cross_region_label"] = assem_labels
    label_df["basepath"] = basepath

    results = {
        "label_df": label_df,
        "assembly_act_task": assembly_act,
        }

    return results

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