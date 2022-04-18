import multiprocessing
from joblib import Parallel, delayed
import os
import pandas as pd
import numpy as np
import nelpy as nel
import pickle
import loading
import sys
from scipy import stats

sys.path.append("D:/github/neurocode/reactivation/assemblies")
import assembly
import functions


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
        epoch_df = loading.load_epoch(self.basepath)
        idx = functions.find_pre_task_post(epoch_df.environment)
        epoch_df = epoch_df[idx[0]]
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
        # zscore and return
        return stats.zscore(z_t.data, axis=1), z_t.bin_centers

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

        self.assembly_act = nel.AnalogSignalArray(
            data=assembly.computeAssemblyActivity(self.patterns, zactmat),
            timestamps=ts,
            fs=1 / self.z_mat_dt,
        )

    def session_loop_activation(self):
        pass

    def standard_task_react(self):
        """
        Runs the standard task reactivation analysis
        """
        self.load_data()
        self.get_weights(self.ripples[self.epochs[1]])
        self.get_assembly_act()


def get_peak_activity(m1, epochs):
    """
    Gets the peak activity of the assembly activity
    """
    strengths = []
    assembly_id = []
    for assembly_act in m1.assembly_act[epochs]:
        strengths.append(assembly_act.max())
        assembly_id.append(np.arange(assembly_act.n_signals))

    return np.hstack(assembly_id), np.hstack(strengths)


def get_pre_post_assembly_strengths(basepath):
    """
    Gets the pre and post assembly strengths
    """
    # initialize session
    m1 = AssemblyReact(basepath)
    # load data
    m1.load_data()
    # restrict to pre/task/post epochs
    m1.restrict_epochs_to_pre_task_post()
    # get weights for task
    m1.get_weights(m1.ripples[m1.epochs[1]])

    
    df = pd.DataFrame()
    for ep in [0, 2]:
        # get assembly activity
        m1.get_assembly_act(epoch=m1.ripples[m1.epochs[ep]])
        # get pre assembly activity
        assembly_id, strengths = get_peak_activity(m1, m1.ripples[m1.epochs[ep]])
        temp_df = pd.DataFrame()
        temp_df["assembly_id"] = assembly_id
        temp_df["strengths"] = strengths
        temp_df["epoch"] = ep
        temp_df["basepath"] = basepath
        temp_df["brainRegion"] = m1.brainRegion
        temp_df["putativeCellType"] = m1.putativeCellType

        df = pd.concat([df, temp_df], ignore_index=True)
    results = {"df": df, "react": m1}
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
