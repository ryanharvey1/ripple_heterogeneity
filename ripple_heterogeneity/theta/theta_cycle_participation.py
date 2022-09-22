import glob
import pandas as pd
import numpy as np
import os
import pickle
import nelpy as nel
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import logging

logging.getLogger().setLevel(logging.ERROR)

def run(basepath):
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion="CA1", putativeCellType="Pyramidal Cell"
    )
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)

    if cell_metrics.shape[0] < 1:
        return None

    try:
        theta_cycles = loading.load_theta_cycles(basepath)
    except:
        return None
    theta_cycles = nel.EpochArray(np.array([theta_cycles.start, theta_cycles.stop]).T)

    # get brain states
    state_dict = loading.load_SleepState_states(basepath)
    theta_epochs = nel.EpochArray(state_dict["THETA"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    wake_theta_epochs = theta_cycles[wake_epochs]

    unit_mat = functions.get_participation(
        st.data, wake_theta_epochs.starts, wake_theta_epochs.stops, par_type="counts"
    )
    partici_prob = (unit_mat > 0).mean(axis=1)
    firing_rate = unit_mat.sum(axis=1) / wake_theta_epochs.duration

    results = pd.DataFrame()
    results["partici_prob"] = partici_prob
    results["firing_rate"] = firing_rate
    results["UID"] = cell_metrics.UID.values
    results["deepSuperficial"] = cell_metrics.deepSuperficial.values
    results["deepSuperficialDistance"] = cell_metrics.deepSuperficialDistance.values
    results["basepath"] = basepath

def load_results(save_path,verbose=False):
    """
    Loads results from save_path.
    Inputs:
        save_path: string, path to save results
    Outputs:
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")

    results_df = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        results_df = pd.concat([results_df, results], ignore_index=True)
        
    return results_df