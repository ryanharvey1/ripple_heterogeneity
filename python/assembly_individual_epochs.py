import nelpy as nel
import functions, loading
import sys

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

    # spike times within ripples
    st_rip = st[ripple_epochs]

    # iter through each behavioral epoch
    patterns = []
    significance = []
    zactmat = []
    for i, ep in enumerate(behavioral_epochs):
        if epoch_df.environment.iloc[i] == "sleep":
            temp_st = st_rip[nrem_epochs][ep]
        else:
            temp_st = st_rip[wake_epochs][ep]

        # extract assembly patterns
        (patterns_, significance_, zactmat_) = assembly.runPatterns(
            temp_st.bin(ds=dt).data
        )

        # store results per epoch
        patterns.append(patterns_)
        significance.append(significance_)
        zactmat.append(zactmat_)

    # package all results in dict
    results = {}
    results["patterns"] = patterns
    results["significance"] = significance
    results["zactmat"] = zactmat

    return results


def session_loop(basepath, save_path, rip_window=0.050):

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
    ripple_epochs = ripple_epochs.expand(rip_window,direction="both")

    # behavioral epochs
    epoch_df = loading.load_epoch(basepath)

    idx = functions.find_epoch_pattern(epoch_df.environment,["sleep","linear","sleep"])
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

    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def assembly_run(df, save_path, parallel=True):
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
