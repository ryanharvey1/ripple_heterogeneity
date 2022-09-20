import nelpy as nel
from ripple_heterogeneity.utils import loading
from ripple_heterogeneity.assembly import assembly, assembly_reactivation
import numpy as np
import pickle
import pandas as pd
from scipy import stats
import os
import multiprocessing
from joblib import Parallel, delayed


def load_basic_data(basepath):
    nChannels, fs, fs_dat, shank_to_channel = loading.loadXML(basepath)
    ripples = loading.load_ripples_events(basepath)
    cell_metrics, data = loading.load_cell_metrics(basepath)
    return cell_metrics, data, ripples, fs_dat


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


def main_analysis(st, ripple_epochs, dt=0.025):
    (
        patterns_outside_ripples,
        significance_outside_ripples,
        zactmat_outside_ripples,
    ) = assembly.runPatterns(st[~ripple_epochs].bin(ds=dt).data)

    (
        patterns_inside_ripples,
        significance_inside_ripples,
        zactmat_inside_ripples,
    ) = assembly.runPatterns(st[ripple_epochs].bin(ds=dt).data)

    results = {}
    results["patterns_outside_ripples"] = patterns_outside_ripples
    results["significance_outside_ripples"] = significance_outside_ripples
    results["zactmat_outside_ripples"] = zactmat_outside_ripples

    results["patterns_inside_ripples"] = patterns_inside_ripples
    results["significance_inside_ripples"] = significance_inside_ripples
    results["zactmat_inside_ripples"] = zactmat_inside_ripples

    return results


def session_loop(basepath, save_path):

    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file):
        return

    cell_metrics, data, ripples, fs_dat = load_basic_data(basepath)

    restrict_idx = (
        (cell_metrics.putativeCellType == "Pyramidal Cell")
        & (cell_metrics.brainRegion.str.contains("CA1"))
        & (cell_metrics.bad_unit == False)
    )

    # restrict cell metrics
    cell_metrics = cell_metrics[restrict_idx]

    if cell_metrics.shape[0] == 0:
        return

    # get ripple epochs
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])
    try:
        st = nel.SpikeTrainArray(
            timestamps=np.array(data["spikes"], dtype=object)[restrict_idx], fs=fs_dat
        )
    except:
        st = nel.SpikeTrainArray(
            timestamps=np.array(data["spikes"], dtype=object)[restrict_idx][0],
            fs=fs_dat,
        )

    results = main_analysis(st, ripple_epochs)

    results["UID"] = cell_metrics.UID
    results["deepSuperficial"] = cell_metrics.deepSuperficial
    results["deepSuperficialDistance"] = cell_metrics.deepSuperficialDistance
    results["basepath"] = basepath

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


def run(basepath, putativeCellType="Pyr", brainRegion="CA1"):

    state_dict = loading.load_SleepState_states(basepath)
    theta_epochs = nel.EpochArray(state_dict["THETA"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    m1 = assembly_reactivation.AssemblyReact(basepath, weight_dt=0.025)
    m1.load_data()
    if m1.st.isempty:
        return None
    m1.get_weights(theta_epochs[wake_epochs])

    results = {}
    results["patterns"] = m1.patterns

    # match outputs from above 
    results["UID"] = m1.cell_metrics.UID
    results["deepSuperficial"] = m1.cell_metrics.deepSuperficial
    results["deepSuperficialDistance"] = m1.cell_metrics.deepSuperficialDistance
    results["brainRegion"] = m1.cell_metrics.brainRegion
    results["basepath"] = basepath

    return results