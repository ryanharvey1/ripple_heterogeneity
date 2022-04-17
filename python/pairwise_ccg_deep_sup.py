import glob
import sys
import pandas as pd
import numpy as np
import os
import pickle
import nelpy as nel
import multiprocessing
from joblib import Parallel, delayed

sys.path.append(r"D:\github\ripple_heterogeneity\python")
import functions, loading, add_new_deep_sup


def load_data(basepath):
    """
    Loads data from a given basepath.
    """
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion="CA1", putativeCellType="Pyramidal Cell"
    )
    # need at least 2 cells to do pairwise CCG
    if cell_metrics.shape[0] < 2:
        return None, None, None, None, None

    cell_metrics = add_new_deep_sup.add_new_deep_sup_class(cell_metrics)
    ripples = loading.load_ripples_events(basepath)
    ripples = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)

    # get brain states                                                
    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict['NREMstate'])
    wake_epochs = nel.EpochArray(state_dict['WAKEstate'])

    return st, cell_metrics, ripples, nrem_epochs, wake_epochs


def main(basepath,states=None):

    # load data
    st, cell_metrics, ripples, nrem_epochs, wake_epochs = load_data(basepath)

    if st is None:
        return None

    # restrict to state
    if states is not None:
        if states == "nrem":
            st = st[nrem_epochs]
            ripples = ripples[nrem_epochs]
        elif states == "wake":
            st = st[wake_epochs]
            ripples = ripples[wake_epochs]
        else:
            raise ValueError("states must be 'nrem' or 'wake'")
        if len(st.data) == 0 | len(ripples.starts ) == 0:
            return None

    unit_mat = functions.get_participation(st.data, ripples.starts, ripples.stops, par_type="counts")
    rho, pval, corr_c = functions.pairwise_corr(unit_mat)

    # get ccg
    ccgs, c = functions.pairwise_cross_corr(
        st[ripples].data, nbins=100, binsize=0.004, return_index=True
    )

    # add id
    ccg_id_df = pd.DataFrame()
    ccg_id_df["ref"] = c[:, 0]
    ccg_id_df["target"] = c[:, 1]
    ccg_id_df["UID_ref"] = cell_metrics.UID.iloc[ccg_id_df["ref"]].values
    ccg_id_df["UID_target"] = cell_metrics.UID.iloc[ccg_id_df["target"]].values
    ccg_id_df["deepSuperficial_ref"] = cell_metrics.deepSuperficial.iloc[
        ccg_id_df["ref"]
    ].values
    ccg_id_df["deepSuperficial_target"] = cell_metrics.deepSuperficial.iloc[
        ccg_id_df["target"]
    ].values
    ccg_id_df["basepath"] = basepath

    results = {
        "ccgs": ccgs,
        "ccg_id_df": ccg_id_df,
        "rho": rho,
        "pval": pval,
        "corr_c": corr_c,
    }

    return results


def session_loop(basepath, df, save_path, states):

    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file):
        return

    results = main(basepath, states)
    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def run(df, save_path, parallel=True, states=None):
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(session_loop)(basepath, df, save_path, states) for basepath in basepaths
        )
    else:
        for basepath in basepaths:
            print(basepath)
            session_loop(basepath, df, save_path, states)


def load_results(save_path):
    sessions = glob.glob(save_path + os.sep + "*.pkl")

    ccgs = pd.DataFrame()
    ccg_id_df = pd.DataFrame()

    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        # horizontally concatenate pandas
        ccgs = pd.concat([ccgs, results["ccgs"]], axis=1)
        # add rho and pval
        results["ccg_id_df"]["rho"] = results["rho"]
        results["ccg_id_df"]["pval"] = results["pval"]
        # vertically concatenate pandas
        ccg_id_df = pd.concat([ccg_id_df, results["ccg_id_df"]], ignore_index=True)
    return ccgs, ccg_id_df
