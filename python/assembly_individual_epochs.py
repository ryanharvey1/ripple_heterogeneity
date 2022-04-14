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
            patterns.append(np.tile(np.nan, len(st.data)))
            significance.append(np.nan)
            zactmat.append(np.nan)
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

    save_path = r"Z:\home\ryanh\projects\ripple_heterogeneity\cell_assembly_epochs_10ms"

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    sessions_df = pd.DataFrame()
    sessions_df["sessions"] = sessions

    df = pd.DataFrame()

    UID = []
    deepSuperficial = []
    deepSuperficialDistance = []
    weights = []
    membership = []
    assembly_n = 0
    assembly_ = []
    basepath = []
    assembly_path = []
    epoch = []
    epoch_n = []
    for session in sessions_df.sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)

        # iterate over epochs
        for i_epoch, patterns in enumerate(results["patterns"]):
            # locate sig assemblies and assembly members
            (
                patterns,
                is_member_keep,
                keep_assembly,
                is_member,
            ) = functions.find_sig_assemblies(patterns)

            # iterate over assemblies
            for i_assemblies, pattern in enumerate(patterns):
                UID.append(results["UID"])
                deepSuperficial.append(results["deepSuperficial"])
                deepSuperficialDistance.append(results["deepSuperficialDistance"])
                weights.append(pattern)
                # thres = np.mean(pattern) + np.std(pattern) * 2
                # membership.append(pattern > thres)
                membership.append(is_member_keep[i_assemblies])
                assembly_.append([assembly_n] * len(pattern))
                assembly_n += 1
                basepath.append([results["basepath"]] * len(pattern))
                assembly_path.append([session] * len(pattern))
                epoch.append([results["env"][i_epoch]] * len(pattern))
                epoch_n.append(np.tile(i_epoch, len(pattern)))


    df["UID"] = np.hstack(UID)
    df["deepSuperficial"] = np.hstack(deepSuperficial)
    df["deepSuperficialDistance"] = np.hstack(deepSuperficialDistance)
    df["epoch"] = np.hstack(epoch)
    df["epoch_n"] = np.hstack(epoch_n)
    df["basepath"] = np.hstack(basepath)
    df["weights"] = np.hstack(weights)
    df["membership"] = np.hstack(membership)
    df["assembly_n"] = np.hstack(assembly_)
    df["assembly_path"] = np.hstack(assembly_path)

    return df
