import glob
import multiprocessing
from joblib import Parallel, delayed
import os
import pandas as pd
import numpy as np
import nelpy as nel
import pickle
from scipy import stats
from ripple_heterogeneity.utils import functions, loading, compress_repeated_epochs
from ripple_heterogeneity.assembly import assembly


class AssemblyReact(object):
    """
    Class for running assembly reactivation analysis

    Parameters:
    -----------
    basepath: str
        Path to the session folder
    brainRegion: str
        Brain region to restrict to
    putativeCellType: str
        Cell type to restrict to
    weight_dt: float
        Time resolution of the weight matrix
    z_mat_dt: float
        Time resolution of the z matrix

    attributes:
    -----------
    st: spike train
    cell_metrics: cell metrics
    ripples: ripples
    zactmat: z matrix
    ts: timestamps
    patterns: patterns
    assembly_act: assembly activity

    methods:
    --------
    load_data: load data
    restrict_to_epoch: restrict to epoch
    restrict_to_events: restrict to events
    get_z_mat: get z matrix
    get_weights: get weights
    get_assembly_act: get assembly activity
    session_loop_activation: session loop activation

    """

    def __init__(
        self,
        basepath,
        brainRegion="CA1",
        putativeCellType="Pyramidal Cell",
        weight_dt=0.01,
        z_mat_dt=0.002,
    ):
        self.basepath = basepath
        self.brainRegion = brainRegion
        self.putativeCellType = putativeCellType
        self.weight_dt = weight_dt
        self.z_mat_dt = z_mat_dt

    def add_st(self, st):
        self.st = st

    def add_ripples(self, ripples):
        self.ripples = ripples

    def add_epoch_df(self, epoch_df):
        self.epoch_df = epoch_df

    def load_spikes(self):
        """
        loads spikes from the session folder
        """
        self.st, self.cell_metrics = loading.load_spikes(
            self.basepath,
            brainRegion=self.brainRegion,
            putativeCellType=self.putativeCellType,
        )

    def load_ripples(self):
        """
        loads ripples from the session folder
        """
        ripples = loading.load_ripples_events(self.basepath)
        self.ripples = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])

    def load_epoch(self):
        """
        loads epochs from the session folder
        """
        epoch_df = loading.load_epoch(self.basepath)
        self.epochs = nel.EpochArray(
            [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
            label=epoch_df.environment.values,
        )

    def load_data(self):
        """
        loads data (spikes,ripples,epochs) from the session folder
        """
        self.load_spikes()
        self.load_ripples()
        self.load_epoch()

    def restrict_epochs_to_pre_task_post(self):
        """
        Restricts the epochs to the specified epochs
        """
        # fetch data
        epoch_df = loading.load_epoch(self.basepath)
        # compress back to back sleep epochs (an issue further up the pipeline)
        epoch_df = compress_repeated_epochs.main(epoch_df)
        # restrict to pre task post epochs
        idx = functions.find_pre_task_post(epoch_df.environment)
        epoch_df = epoch_df[idx[0]]
        # convert to epoch array and add to object
        self.epochs = nel.EpochArray(
            [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
            label=epoch_df.environment.values,
        )

    def restrict_to_epoch(self, epoch):
        """
        Restricts the spike data to a specific epoch
        """
        self.st_resticted = self.st[epoch]

    def get_z_mat(self, st):
        """
        To increase the temporal resolution beyond the
        bin-size used to identify the assembly patterns,
        z(t) was obtained by convolving the spike-train
        of each neuron with a kernel-function
        """
        # binning the spike train
        z_t = st.bin(ds=self.z_mat_dt)
        # gaussian kernel to match the bin-size used to identify the assembly patterns
        sigma = self.weight_dt / np.sqrt(int(1000 * self.weight_dt / 2))
        z_t.smooth(sigma=sigma, inplace=True)
        # zscore the z matrix
        z_scored_bst = stats.zscore(z_t.data, axis=1)
        # make sure there are no nans, important as strengths will all be nan otherwise
        z_scored_bst[np.isnan(z_scored_bst).any(axis=1)] = 0

        return z_scored_bst, z_t.bin_centers

    def get_weights(self, epoch=None):
        """
        Gets the assembly weights
        """
        if epoch is not None:
            bst = self.st[epoch].bin(ds=self.weight_dt).data
        else:
            bst = self.st.bin(ds=self.weight_dt).data

        self.patterns, _, _ = assembly.runPatterns(bst)

    def get_assembly_act(self, epoch=None):
        if epoch is not None:
            zactmat, ts = self.get_z_mat(self.st[epoch])
        else:
            zactmat, ts = self.get_z_mat(self.st)

        assembly_act = nel.AnalogSignalArray(
            data=assembly.computeAssemblyActivity(self.patterns, zactmat),
            timestamps=ts,
            fs=1 / self.z_mat_dt,
        )
        return assembly_act


def get_peak_activity(assembly_act, epochs):
    """
    Gets the peak activity of the assembly activity
    """
    strengths = []
    assembly_id = []
    for assembly_act in assembly_act[epochs]:
        strengths.append(assembly_act.max())
        assembly_id.append(np.arange(assembly_act.n_signals))

    return np.hstack(assembly_id), np.hstack(strengths)


def get_pre_post_assembly_strengths(basepath):
    """
    Gets the pre and post assembly strengths
    """
    # initialize session
    m1 = AssemblyReact(basepath, weight_dt=0.02)
    # load data
    m1.load_data()
    # check if no cells were found
    if m1.cell_metrics.shape[0] == 0:
        return None
    # restrict to pre/task/post epochs
    m1.restrict_epochs_to_pre_task_post()
    # get weights for task outside ripples
    # % (TODO: use more robust method to locate epochs than index)
    m1.get_weights(m1.epochs[1][~m1.ripples])

    # get assembly activity
    assembly_act_pre = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[0]])
    assembly_act_task = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[1]])
    assembly_act_post = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[2]])
    results = {
        "assembly_act_pre": assembly_act_pre,
        "assembly_act_task": assembly_act_task,
        "assembly_act_post": assembly_act_post,
        "react": m1,
    }

    return results


def session_loop(basepath, save_path):

    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file):
        return
    results = get_pre_post_assembly_strengths(basepath)
    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def run(df, save_path, parallel=True):
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


def load_results(save_path):
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    all_results = {}
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
            if results is None:
                continue
        all_results[results["react"].basepath] = results
    return all_results
