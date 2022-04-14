import nelpy as nel
import functions, loading
import sys
import glob

sys.path.append(r"D:\github\neurocode\reactivation\assemblies")
import assembly
import numpy as np
import pickle
import pandas as pd
from scipy import stats
import os
import multiprocessing
from joblib import Parallel, delayed


def get_z_t(st, ds=0.001):
    """
    To increase the temporal resolution beyond the bin-size used to identify the assembly patterns,
    z(t) was obtained by convolving the spike-train of each neuron with a kernel-function
    """
    # bin to 1ms
    z_t = st.bin(ds=ds)
    # make binary
    z_t.data[z_t.data > 1] = 1
    # gaussian kernel to match the bin-size used to identify the assembly patterns
    z_t.smooth(sigma=0.025 / np.sqrt(12), inplace=True)
    # zscore
    return stats.zscore(z_t.data, axis=1), z_t.bin_centers


def main_analysis(
    st, ripple_epochs, behavioral_epochs, epoch_df, nrem_epochs, wake_epochs, dt=0.010
):
    """
    Input:
        st: spike train
        ripple_epochs: ripple epochs
        behavioral_epochs: behavioral epochs
        epoch_df: epoch dataframe
        nrem_epochs: nrem epochs
        wake_epochs: wake epochs
        dt: time bin size
    Output:
        results: dictionary of results
    """
    # spike times within ripples
    st_rip = st[ripple_epochs]

    # iter through each behavioral epoch
    patterns = []
    significance = []
    zactmat = []
    env = []
    for i, ep in enumerate(behavioral_epochs):
        if epoch_df.environment.iloc[i] == "sleep":
            temp_st = st_rip[nrem_epochs][ep]
        else:
            temp_st = st_rip[wake_epochs][ep]
        # check if temp_st is empty
        if temp_st.data is None:
            patterns.append(None)
            significance.append(None)
            zactmat.append(None)
            env.append(epoch_df.environment.iloc[i])
            continue
        # extract assembly patterns
        (patterns_, significance_, zactmat_) = assembly.runPatterns(
            temp_st.bin(ds=dt).data
        )

        # store results per epoch
        patterns.append(patterns_)
        significance.append(significance_)
        zactmat.append(zactmat_)
        env.append(epoch_df.environment.iloc[i])

    # package all results in dict
    results = {}
    results["patterns"] = patterns
    results["significance"] = significance
    results["zactmat"] = zactmat
    results["env"] = env

    return results


def session_loop(basepath, save_path, rip_window=0.050):
    """
    Loop through each session
    Input:
        basepath: path to the session
        save_path: path to save the results
        rip_window: time window to define ripples
    """
    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file):
        return

    _, _, ripples, fs_dat = loading.load_basic_data(basepath)

    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion="CA1", putativeCellType="Pyramidal"
    )

    if cell_metrics.shape[0] == 0:
        return

    # get ripple epochs
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])
    ripple_epochs = ripple_epochs.expand(rip_window, direction="both")

    # behavioral epochs
    epoch_df = loading.load_epoch(basepath)

    idx = functions.find_epoch_pattern(
        epoch_df.environment, ["sleep", "linear", "sleep"]
    )
    epoch_df = epoch_df[idx[0]]

    behavioral_epochs = nel.EpochArray(
        [np.array([epoch_df.startTime, epoch_df.stopTime]).T]
    )

    # get brain states
    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict["NREMstate"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    results = main_analysis(
        st, ripple_epochs, behavioral_epochs, epoch_df, nrem_epochs, wake_epochs
    )

    results["UID"] = cell_metrics.UID
    results["basepath"] = basepath
    results["deepSuperficial"] = cell_metrics.deepSuperficial
    results["deepSuperficialDistance"] = cell_metrics.deepSuperficialDistance

    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def assembly_run(df, save_path, parallel=True):
    """
    Run the assembly analysis on all sessions in the dataframe
    """
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(session_loop)(basepath, save_path) for basepath in basepaths
        )
    else:
        for basepath in basepaths:
            print(basepath)
            session_loop(basepath, save_path)


def load_assem_epoch_data(save_path):
    """
    Loads the assembly epoch data from the save_path
    """

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    assem_epoch_df = pd.DataFrame()

    for session in sessions:
        assem_epoch_df_temp = pd.DataFrame()
        prob_sig_member = []
        n_members = []
        n_assemblies = []
        n_cells = []
        epoch = []
        with open(session, "rb") as f:
            results = pickle.load(f)

            for i, pattern_ep in enumerate(results["patterns"]):
                (
                    patterns_keep,
                    is_member_keep,
                    keep_assembly,
                    is_member,
                ) = functions.find_sig_assemblies(pattern_ep)
                prob_sig_member.append(np.mean(is_member_keep))
                n_members.append(is_member_keep.sum())
                n_assemblies.append(patterns_keep.shape[0])
                n_cells.append(patterns_keep.shape[1])
                epoch.append(results['env'][i])

            assem_epoch_df_temp["prob_sig_member"] = prob_sig_member
            assem_epoch_df_temp["n_members"] = n_members
            assem_epoch_df_temp["n_assemblies"] = n_assemblies
            assem_epoch_df_temp["n_cells"] = n_cells
            assem_epoch_df_temp["epoch"] = epoch
            assem_epoch_df_temp["basepath"] = results["basepath"]

        assem_epoch_df = pd.concat(
            [assem_epoch_df, assem_epoch_df_temp], ignore_index=True
        )

    return assem_epoch_df
