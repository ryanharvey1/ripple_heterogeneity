from logging import exception
import scipy.io as sio
import sys, os
import pandas as pd
import numpy as np
import glob
import nelpy as nel
import warnings
from ripple_heterogeneity.utils import functions
from warnings import simplefilter
from typing import Union
import multiprocessing
from joblib import Parallel, delayed
from xml.dom import minidom
from scipy import signal

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


def loadXML(basepath: str):
    """
    path should be the folder session containing the XML file
    Function returns :
        1. the number of channels
        2. the sampling frequency of the dat file or the eeg file depending of what is present in the folder
            eeg file first if both are present or both are absent
        3. the mappings shanks to channels as a dict
    Args:
        path : string
    Returns:
        int, int, dict

    by Guillaume Viejo
    """
    # check if saved file exists
    try:
        basename = os.path.basename(basepath)
        filename = glob.glob(os.path.join(basepath, basename + ".xml"))[0]
    except:
        warnings.warn("xml file does not exist")
        return

    xmldoc = minidom.parse(filename)
    nChannels = (
        xmldoc.getElementsByTagName("acquisitionSystem")[0]
        .getElementsByTagName("nChannels")[0]
        .firstChild.data
    )
    fs_dat = (
        xmldoc.getElementsByTagName("acquisitionSystem")[0]
        .getElementsByTagName("samplingRate")[0]
        .firstChild.data
    )
    fs = (
        xmldoc.getElementsByTagName("fieldPotentials")[0]
        .getElementsByTagName("lfpSamplingRate")[0]
        .firstChild.data
    )

    shank_to_channel = {}
    groups = (
        xmldoc.getElementsByTagName("anatomicalDescription")[0]
        .getElementsByTagName("channelGroups")[0]
        .getElementsByTagName("group")
    )
    for i in range(len(groups)):
        shank_to_channel[i] = [
            int(child.firstChild.data)
            for child in groups[i].getElementsByTagName("channel")
        ]
    return int(nChannels), int(fs), int(fs_dat), shank_to_channel


def loadLFP(
    basepath: str,
    n_channels: int = 90,
    channel: int = 64,
    frequency: float = 1250.0,
    precision: str = "int16",
    ext: str = "lfp",
):
    path = ""
    if ext == "lfp":
        try:
            path = glob.glob(
                os.path.join(basepath, os.path.basename(basepath) + "*.lfp")
            )[0]
        except:
            path = glob.glob(
                os.path.join(basepath, os.path.basename(basepath) + "*.eeg")
            )[0]
    if ext == "dat":
        path = glob.glob(os.path.join(basepath, os.path.basename(basepath) + "*.dat"))[
            0
        ]

    # check if saved file exists
    if not os.path.exists(path):
        warnings.warn("file does not exist")
        return
    if channel is None:
        n_channels = int(n_channels)

        f = open(path, "rb")
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2
        n_samples = int((endoffile - startoffile) / n_channels / bytes_size)
        duration = n_samples / frequency
        f.close()
        data = np.memmap(path, np.int16, "r", shape=(n_samples, n_channels))
        timestep = np.arange(0, n_samples) / frequency
        return data, timestep

    if type(channel) is not list:
        f = open(path, "rb")
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2
        n_samples = int((endoffile - startoffile) / n_channels / bytes_size)
        duration = n_samples / frequency
        interval = 1 / frequency
        f.close()
        with open(path, "rb") as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:, channel]
            timestep = np.arange(0, len(data)) / frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(
                os.path.dirname(os.path.abspath(path)), "lfp_ts.npy"
            )
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)

            return data, timestep

    elif type(channel) is list:
        f = open(path, "rb")
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2

        n_samples = int((endoffile - startoffile) / n_channels / bytes_size)
        duration = n_samples / frequency
        f.close()
        with open(path, "rb") as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:, channel]
            timestep = np.arange(0, len(data)) / frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(
                os.path.dirname(os.path.abspath(path)), "lfp_ts.npy"
            )
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)
            return data, timestep


class LoadLfp(object):
    """
    Simple class to load LFP or wideband data from a recording folder
    Args:
        basepath : string (path to the recording folder)
        channels : int or list of int or None (default None, load all channels memmap)
        ext : string (lfp or dat)
        epoch: nelpy EpochArray or ndarray (default None, load all data)
    Returns:
        : nelpy analogsignalarray of shape (n_channels, n_samples)

    Example:
        # load lfp file
        >>> basepath = r"X:/data/Barrage/NN10/day10"
        >>> lfp = loading.LoadLfp(basepath,ext="lfp")
        >>> lfp
        >>>    <AnalogSignalArray at 0x25ba1576640: 128 signals> for a total of 5:33:58:789 hours

        # Loading dat file
        >>> dat = loading.LoadLfp(basepath,ext="dat")
        >>> dat
        >>>    <AnalogSignalArray at 0x25ba4fedc40: 128 signals> for a total of 5:33:58:790 hours
        >>> dat.lfp.data.shape
        >>>    (128, 400775808)
        >>> type(dat.lfp.data)
        >>>    numpy.memmap

    Ryan Harvey 2023
    """

    def __init__(
        self,
        basepath: str,
        channels: Union[int, list, None] = None,
        ext: str = "lfp",
        epoch:Union[np.ndarray,nel.EpochArray,None]=None,
    ) -> None:
        self.basepath = basepath  # path to the recording folder
        self.channels = channels  # channel number or list of channel numbers
        self.ext = ext  # lfp or dat
        self.epoch = epoch

        # get xml data
        self.get_xml_data()

        # set sampling rate based on the extension of the file (lfp or dat)
        if self.ext == "dat":
            self.fs = self.fs_dat

        # load lfp
        self.load_lfp()

    def get_xml_data(self):
        nChannels, fs, fs_dat, shank_to_channel = loadXML(self.basepath)
        self.nChannels = nChannels
        self.fs = fs
        self.fs_dat = fs_dat
        self.shank_to_channel = shank_to_channel

    def load_lfp(self):
        lfp, timestep = loadLFP(
            self.basepath,
            n_channels=self.nChannels,
            channel=self.channels,
            frequency=self.fs,
            ext=self.ext,
        )

        if isinstance(self.epoch, nel.EpochArray):
            intervals = self.epoch.data
        elif isinstance(self.epoch, np.ndarray):
            intervals = self.epoch
            if intervals.ndim == 1:
                intervals = intervals[np.newaxis, :]
        else:
            intervals = np.array([[min(timestep), max(timestep)]])

        idx = functions.in_intervals(timestep, intervals)

        self.lfp = nel.AnalogSignalArray(
            data=lfp[idx, None].T,
            timestamps=timestep[idx],
            fs=self.fs,
            support=nel.EpochArray(np.array([min(timestep[idx]), max(timestep[idx])])),
        )

    def __repr__(self) -> None:
        return self.lfp.__repr__()

    def get_phase(self, band2filter: list = [6, 12], ford=3):
        band2filter = np.array(band2filter, dtype=float)
        b, a = signal.butter(ford, band2filter / (self.fs / 2), btype="bandpass")
        filt_sig = signal.filtfilt(b, a, self.lfp.data, padtype="odd")
        return np.angle(signal.hilbert(filt_sig))

    def get_freq_phase_amp(self, band2filter: list = [6, 12], ford=3):
        """
        Uses the Hilbert transform to calculate the instantaneous phase and
        amplitude of the time series in sig
        Parameters
        ----------
        sig: np.array
            The signal to be analysed
        ford: int
            The order for the Butterworth filter
        band2filter: list
            The two frequencies to be filtered for e.g. [6, 12]
        """

        band2filter = np.array(band2filter, dtype=float)

        b, a = signal.butter(ford, band2filter / (self.fs / 2), btype="bandpass")

        filt_sig = signal.filtfilt(b, a, self.lfp.data, padtype="odd")
        phase = np.angle(signal.hilbert(filt_sig))
        amplitude = np.abs(signal.hilbert(filt_sig))
        amplitude_filtered = signal.filtfilt(b, a, amplitude, padtype="odd")
        return filt_sig, phase, amplitude, amplitude_filtered


def load_position(basepath, fs=39.0625):
    if not os.path.exists(basepath):
        print("The path " + basepath + " doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(basepath)
    whlfiles = [f for f in listdir if f.endswith(".whl")]
    if not len(whlfiles):
        print("Folder contains no whl files; Exiting ...")
        sys.exit()
    new_path = os.path.join(basepath, whlfiles[0])
    df = pd.read_csv(new_path, delimiter="\t", header=0, names=["x1", "y1", "x2", "y2"])
    df[df == -1] = np.nan
    return df, fs


def load_all_cell_metrics(basepaths):
    """
    load cell metrics from multiple sessions

    Input:
            basepaths: list of basepths, can be pandas column
    Output:
            cell_metrics: concatenated pandas df with metrics

    Note: to get waveforms, spike times, etc. use load_cell_metrics
    """

    # to speed up, use parallel
    num_cores = multiprocessing.cpu_count()
    cell_metrics = Parallel(n_jobs=num_cores)(
        delayed(load_cell_metrics)(basepath, True) for basepath in basepaths
    )

    return pd.concat(cell_metrics, ignore_index=True)


def load_cell_metrics(basepath: str, only_metrics: bool = False) -> tuple:
    """
    loader of cell-explorer cell_metrics.cellinfo.mat

    Inputs: basepath: path to folder with cell_metrics.cellinfo.mat
    outputs: df: data frame of single unit features
    data_: dict with data that does not fit nicely into a dataframe (waveforms, acgs, epochs, etc.)

    See https://cellexplorer.org/datastructure/standard-cell-metrics/ for details

    TODO: extract all fields from cell_metrics.cellinfo. There are more items that can be extracted

    - Ryan H
    """

    def extract_epochs(data):
        startTime = [
            ep["startTime"][0][0][0][0]
            for ep in data["cell_metrics"]["general"][0][0]["epochs"][0][0][0]
        ]
        stopTime = [
            ep["stopTime"][0][0][0][0]
            for ep in data["cell_metrics"]["general"][0][0]["epochs"][0][0][0]
        ]
        name = [
            ep["name"][0][0][0]
            for ep in data["cell_metrics"]["general"][0][0]["epochs"][0][0][0]
        ]

        epochs = pd.DataFrame()
        epochs["name"] = name
        epochs["startTime"] = startTime
        epochs["stopTime"] = stopTime
        return epochs

    def extract_general(data):
        # extract fr per unit with lag zero to ripple
        try:
            ripple_fr = [
                ev.T[0]
                for ev in data["cell_metrics"]["events"][0][0]["ripples"][0][0][0]
            ]
        except:
            ripple_fr = []
        # extract spikes times
        spikes = [
            spk.T[0] for spk in data["cell_metrics"]["spikes"][0][0]["times"][0][0][0]
        ]
        # extract epochs
        try:
            epochs = extract_epochs(data)
        except:
            epochs = []
        # extract avg waveforms
        try:
            waveforms = np.vstack(
                data["cell_metrics"]["waveforms"][0][0]["filt"][0][0][0]
            )
        except:
            try:
                waveforms = [
                    w.T for w in data["cell_metrics"]["waveforms"][0][0][0][0][0][0]
                ]
            except:
                waveforms = [w.T for w in data["cell_metrics"]["waveforms"][0][0][0]]
        # extract chanCoords
        try:
            chanCoords_x = data["cell_metrics"]["general"][0][0]["chanCoords"][0][0][0][
                0
            ]["x"].T[0]
            chanCoords_y = data["cell_metrics"]["general"][0][0]["chanCoords"][0][0][0][
                0
            ]["y"].T[0]
        except:
            chanCoords_x = []
            chanCoords_y = []

        # add to dictionary
        data_ = {
            "acg_wide": data["cell_metrics"]["acg"][0][0]["wide"][0][0],
            "acg_narrow": data["cell_metrics"]["acg"][0][0]["narrow"][0][0],
            "acg_log10": data["cell_metrics"]["acg"][0][0]["log10"][0][0],
            "ripple_fr": ripple_fr,
            "chanCoords_x": chanCoords_x,
            "chanCoords_y": chanCoords_y,
            "epochs": epochs,
            "spikes": spikes,
            "waveforms": waveforms,
        }
        return data_

    def un_nest_df(df):
        # Un-nest some strings are nested within brackets (a better solution exists...)
        # locate and iterate objects in df
        for item in df.keys()[df.dtypes == "object"]:
            # if you can get the size of the first item with [0], it is nested
            # otherwise it fails and is not nested
            try:
                df[item][0][0].size
                # the below line is from: https://www.py4u.net/discuss/140913
                df[item] = df[item].str.get(0)
            except:
                continue
        return df

    filename = os.path.join(basepath, os.path.basename(basepath) + ".cell_metrics.cellinfo.mat")
    # filename = glob.glob(os.path.join(basepath, "*.cell_metrics.cellinfo.mat"))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        if only_metrics:
            return None
        return None, None

    # load cell_metrics file
    data = sio.loadmat(filename)

    # construct data frame with features per neuron
    df = pd.DataFrame()
    # count units
    n_cells = data["cell_metrics"]["UID"][0][0][0].size
    dt = data["cell_metrics"].dtype
    for dn in dt.names:
        # check if var has the right n of units and is a vector
        try:
            if (data["cell_metrics"][dn][0][0][0][0].size == 1) & (
                data["cell_metrics"][dn][0][0][0].size == n_cells
            ):

                df[dn] = data["cell_metrics"][dn][0][0][0]
        except:
            continue

    # load in tag
    # check if tags exist within cell_metrics
    if "tags" in data.get("cell_metrics").dtype.names:
        # get names of each tag
        dt = data["cell_metrics"]["tags"][0][0].dtype
        if len(dt) > 0:
            # iter through each tag
            for dn in dt.names:
                # set up column for tag
                df["tags_" + dn] = [False] * df.shape[0]
                # iter through uid
                for uid in data["cell_metrics"]["tags"][0][0][dn][0][0].flatten():
                    df.loc[df.UID == uid, "tags_" + dn] = True

    # add bad unit tag for legacy
    df["bad_unit"] = [False] * df.shape[0]
    if "tags_Bad" in df.keys():
        df.bad_unit = df.tags_Bad
        df.bad_unit = df.bad_unit.replace({np.nan: False})

    # add data from general metrics
    df["basename"] = data["cell_metrics"]["general"][0][0]["basename"][0][0][0]
    df["basepath"] = basepath
    df["sex"] = data["cell_metrics"]["general"][0][0]["animal"][0][0]["sex"][0][0][0]
    df["species"] = data["cell_metrics"]["general"][0][0]["animal"][0][0]["species"][0][
        0
    ][0]
    df["strain"] = data["cell_metrics"]["general"][0][0]["animal"][0][0]["strain"][0][
        0
    ][0]
    try:
        df["geneticLine"] = data["cell_metrics"]["general"][0][0]["animal"][0][0][
            "geneticLine"
        ][0][0][0]
    except:
        pass
    df["cellCount"] = data["cell_metrics"]["general"][0][0]["cellCount"][0][0][0][0]

    # fix nesting issue for strings
    df = un_nest_df(df)

    # convert nans within tags columns to false
    cols = df.filter(regex="tags_").columns
    df[cols] = df[cols].replace({np.nan: False})

    if only_metrics:
        return df

    # extract other general data and put into dict
    data_ = extract_general(data)

    return df, data_


def load_SWRunitMetrics(basepath):
    """
    load_SWRunitMetrics loads SWRunitMetrics.mat into pandas dataframe

    returns pandas dataframe with the following fields
        particip: the probability of participation into ripples for each unit
        FRall: mean firing rate during ripples
        FRparticip: mean firing rate for ripples with at least 1 spike
        nSpkAll: mean number of spikes in all ripples
        nSpkParticip: mean number of spikes in ripples with at least 1 spike
        epoch: behavioral epoch label
    """

    def extract_swr_epoch_data(data, epoch):
        # get var names
        dt = data["SWRunitMetrics"][epoch][0][0].dtype

        df2 = pd.DataFrame()

        # get n units
        # there might be other fields within here like the epoch timestamps
        # skip those by returning empty df
        try:
            n_cells = data["SWRunitMetrics"][epoch][0][0][0]["particip"][0].shape[0]
        except:
            return df2

        for dn in dt.names:
            if (data["SWRunitMetrics"][epoch][0][0][0][dn][0].shape[1] == 1) & (
                data["SWRunitMetrics"][epoch][0][0][0][dn][0].shape[0] == n_cells
            ):
                df2[dn] = data["SWRunitMetrics"][epoch][0][0][0][dn][0].T[0]
        df2["epoch"] = epoch
        return df2

    try:
        filename = glob.glob(os.path.join(basepath, "*.SWRunitMetrics.mat"))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)

    df2 = pd.DataFrame()
    # loop through each available epoch and pull out contents
    for epoch in data["SWRunitMetrics"].dtype.names:
        if data["SWRunitMetrics"][epoch][0][0].size > 0:  # not empty

            # call content extractor
            df_ = extract_swr_epoch_data(data, epoch)

            # append conents to overall data frame
            if df_.size > 0:
                df2 = pd.concat([df2, df_], ignore_index=True)

    return df2


def add_manual_events(df, added_ts):
    """
    Add new rows to a dataframe representing manual events (from Neuroscope2)
    with durations equal to the mean duration of the existing events.

    Parameters:
    df (pandas DataFrame): The input dataframe, with at least two columns called
        'start' and 'stop', representing the start and stop times of the events

    added_ts (list): A list of timestamps representing the peaks of the new
        events to be added to the dataframe.

    Returns:
    pandas DataFrame: The modified dataframe with the new rows added and sorted by the 'peaks' column.
    """
    # Calculate the mean duration of the existing events
    mean_duration = (df["stop"] - df["start"]).mean()

    # Create a new dataframe with a 'peaks' column equal to the added_ts values
    df_added = pd.DataFrame()
    df_added["peaks"] = added_ts

    # Calculate the start and stop times of the new events based on the mean duration
    df_added["start"] = added_ts - mean_duration / 2
    df_added["stop"] = added_ts + mean_duration / 2

    # Calculate the duration of the new events as the mean duration
    df_added["duration"] = df_added.stop.values - df_added.start.values

    # Append the new events to the original dataframe
    df = pd.concat([df, df_added], ignore_index=True)

    # Sort the dataframe by the 'peaks' column
    df.sort_values(by=["peaks"], ignore_index=True, inplace=True)

    return df


def load_ripples_events(
    basepath: str, return_epoch_array: bool = False, manual_events: bool = True
):
    """
    load info from ripples.events.mat and store within df

    args:
        basepath: path to your session where ripples.events.mat is
        return_epoch_array: if you want the output in an EpochArray
        manual_events: add manually added events from Neuroscope2
            (interval will be calculated from mean event duration)

    returns pandas dataframe with the following fields
        start: start time of ripple
        stop: end time of ripple
        peaks: peak time of ripple
        amplitude: envlope value at peak time
        duration: ripple duration
        frequency: insta frequency at peak
        detectorName: the name of ripple detector used
        event_spk_thres: 1 or 0 for if a mua thres was used
        basepath: path name
        basename: session id
        animal: animal id *

        * Note that basepath/basename/animal relies on specific folder
        structure and may be incorrect for some data structures
    """

    # locate .mat file
    try:
        filename = glob.glob(basepath + os.sep + "*ripples.events.mat")[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename)

    # make data frame of known fields
    df = pd.DataFrame()
    try:
        df["start"] = data["ripples"]["timestamps"][0][0][:, 0]
        df["stop"] = data["ripples"]["timestamps"][0][0][:, 1]
    except:
        df["start"] = data["ripples"]["times"][0][0][:, 0]
        df["stop"] = data["ripples"]["times"][0][0][:, 1]

    for name in ["peaks", "amplitude", "duration", "frequency", "peakNormedPower"]:
        try:
            df[name] = data["ripples"][name][0][0]
        except:
            df[name] = np.nan

    if df.duration.isna().all():
        df["duration"] = df.stop - df.start

    try:
        df["detectorName"] = data["ripples"]["detectorinfo"][0][0]["detectorname"][0][
            0
        ][0]
    except:
        df["detectorName"] = data["ripples"]["detectorName"][0][0][0]

    # find ripple channel (this can be in several places depending on the file)
    try:
        df["ripple_channel"] = data["ripples"]["detectorinfo"][0][0]["detectionparms"][
            0
        ][0]["Channels"][0][0][0][0]
    except:
        try:
            df["ripple_channel"] = data["ripples"]["detectorParams"][0][0]["channel"][
                0
            ][0][0][0]
        except:
            try:
                df["ripple_channel"] = data["ripples"]["detectorinfo"][0][0][
                    "detectionparms"
                ][0][0]["channel"][0][0][0][0]
            except:
                try:
                    df["ripple_channel"] = data["ripples"]["detectorinfo"][0][0][
                        "detectionparms"
                    ][0][0]["ripple_channel"][0][0][0][0]
                except:
                    df["ripple_channel"] = data["ripples"]["detectorinfo"][0][0][
                        "detectionchannel1"
                    ][0][0][0][0]

    # remove flagged ripples, if exist
    try:
        df.drop(
            labels=np.array(data["ripples"]["flagged"][0][0]).T[0] - 1,
            axis=0,
            inplace=True,
        )
        df.reset_index(inplace=True)
    except:
        pass

    # adding manual events
    if manual_events:
        try:
            df = add_manual_events(df, data["ripples"]["added"][0][0].T[0])
        except:
            pass

    # adding if ripples were restricted by spikes
    dt = data["ripples"].dtype
    if "eventSpikingParameters" in dt.names:
        df["event_spk_thres"] = 1
    else:
        df["event_spk_thres"] = 0

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    if return_epoch_array:
        return nel.EpochArray([np.array([df.start, df.stop]).T], label="ripples")

    return df


def load_theta_cycles(basepath, return_epoch_array=False):
    """
    load theta cycles calculated from auto_theta_cycles.m
    """
    filename = glob.glob(os.path.join(basepath, "*.thetacycles.events.mat"))[0]
    data = sio.loadmat(filename, simplify_cells=True)
    df = pd.DataFrame()
    df["start"] = data["thetacycles"]["timestamps"][:, 0]
    df["stop"] = data["thetacycles"]["timestamps"][:, 1]
    df["duration"] = data["thetacycles"]["duration"]
    df["center"] = data["thetacycles"]["center"]
    df["trough"] = data["thetacycles"]["peaks"]
    df["theta_channel"] = data["thetacycles"]["detectorinfo"]["theta_channel"]

    if return_epoch_array:
        return nel.EpochArray([np.array([df.start, df.stop]).T], label="theta_cycles")
    return df


def load_barrage_events(
    basepath, return_epoch_array=False, restrict_to_nrem=True, file_name=None
):
    """
    load info from barrage.events.mat and store within df

    basepath: path to your session where ripples.events.mat is

    returns pandas dataframe with the following fields
        start: start time of ripple
        stop: end time of ripple
        peaks: peak time of ripple
        amplitude: envlope value at peak time
        duration: ripple duration
        basepath: path name
        basename: session id
        animal: animal id *

        * Note that basepath/basename/animal relies on specific folder
        structure and may be incorrect for some data structures
    """

    # locate .mat file
    try:
        if file_name is not None:
            filename = glob.glob(
                basepath + os.sep + "Barrage_Files" + os.sep + file_name + "*.mat"
            )[0]
        else:
            filename = glob.glob(
                os.path.join(
                    basepath, "Barrage_Files", os.path.basename(basepath) + ".HSE.mat"
                )
            )[0]
    except:
        warnings.warn("file does not exist")
        if return_epoch_array:
            return nel.EpochArray()
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename, simplify_cells=True)

    # make data frame of known fields
    df = pd.DataFrame()
    df["start"] = data["HSE"]["timestamps"][:, 0]
    df["stop"] = data["HSE"]["timestamps"][:, 1]
    df["peaks"] = data["HSE"]["peaks"]
    try:
        df["amplitude"] = data["HSE"]["amplitudes"]
    except:
        df["amplitude"] = np.nan
    try:
        df["duration"] = data["HSE"]["duration"]
    except:
        df["duration"] = df["stop"] - df["start"]

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    df = df.loc[np.array(data["HSE"]["keep"]).T - 1].reset_index(drop=True)

    if restrict_to_nrem:
        df = df.loc[np.array(data["HSE"]["NREM"]).T - 1].reset_index(drop=True)

    if return_epoch_array:
        return nel.EpochArray([np.array([df.start, df.stop]).T], label="barrage")

    return df


def load_ied_events(basepath: str, return_epoch_array: bool = False):
    """
    load info from ripples.events.mat and store within df

    args:
        basepath: path to your session where ripples.events.mat is
        return_epoch_array: if you want the output in an EpochArray
        manual_events: add manually added events from Neuroscope2
            (interval will be calculated from mean event duration)

    returns pandas dataframe

        * Note that basepath/basename/animal relies on specific folder
        structure and may be incorrect for some data structures
    """

    # locate .mat file
    try:
        filename = glob.glob(basepath + os.sep + "*IED.events.mat")[0]
    except:
        # warnings.warn("file does not exist")
        return pd.DataFrame()

    df = pd.DataFrame()

    data = sio.loadmat(filename, simplify_cells=True)
    struct_name = list(data.keys())[-1]
    df["start"] = data[struct_name]["timestamps"][:, 0]
    df["stop"] = data[struct_name]["timestamps"][:, 1]
    df["center"] = data[struct_name]["peaks"]

    if return_epoch_array:
        return nel.EpochArray([np.array([df.start, df.stop]).T], label="ied")

    return df


def load_dentate_spike(basepath):
    """
    load info from DS*.events.mat and store within df
    basepath: path to your session where DS*.events.mat is

    returns pandas dataframe with the following fields
        start: start time of DS
        stop: end time of DS
        peaks: peak time of DS
        amplitude: envlope value at peak time
        duration: DS duration
        detectorName: the name of DS detector used
        basepath: path name
        basename: session id
        animal: animal id *
        * Note that basepath/basename/animal relies on specific folder
        structure and may be incorrect for some data structures
    """

    def extract_data(s_type, data):
        # make data frame of known fields
        df = pd.DataFrame()
        df["start"] = data[s_type]["timestamps"][:, 0]
        df["stop"] = data[s_type]["timestamps"][:, 1]
        df["peaks"] = data[s_type]["peaks"]
        df["event_label"] = s_type
        df["amplitude"] = data[s_type]["amplitudes"]
        df["duration"] = data[s_type]["duration"]
        df["amplitudeUnits"] = data[s_type]["amplitudeUnits"]
        df["detectorName"] = data[s_type]["detectorinfo"]["detectorname"]
        df["ml_channel"] = data[s_type]["detectorinfo"]["ml_channel"]
        df["h_channel"] = data[s_type]["detectorinfo"]["h_channel"]
        return df

    # locate .mat file
    df = pd.DataFrame()
    for s_type in ["DS1", "DS2"]:
        filename = glob.glob(basepath + os.sep + "*" + s_type + ".events.mat")
        if len(filename) == 0:
            continue
        # load matfile
        filename = filename[0]
        data = sio.loadmat(filename, simplify_cells=True)
        # pull out data
        df = pd.concat([df, extract_data(s_type, data)], ignore_index=True)

    if df.shape[0] == 0:
        return df

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    return df


def load_theta_rem_shift(basepath):
    """
    load_theta_rem_shift: loads matlab structure from get_rem_shift.m
    """
    try:
        filename = glob.glob(basepath + os.sep + "*theta_rem_shift.mat")[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame(), np.nan

    data = sio.loadmat(filename)

    df = pd.DataFrame()

    df["UID"] = data["rem_shift_data"]["UID"][0][0][0]
    df["circ_dist"] = data["rem_shift_data"]["circ_dist"][0][0][0]
    df["rem_shift"] = data["rem_shift_data"]["rem_shift"][0][0][0]
    df["non_rem_shift"] = data["rem_shift_data"]["non_rem_shift"][0][0][0]

    # rem metrics
    df["m_rem"] = data["rem_shift_data"]["PhaseLockingData_rem"][0][0]["phasestats"][0][
        0
    ]["m"][0][0][0]
    df["r_rem"] = data["rem_shift_data"]["PhaseLockingData_rem"][0][0]["phasestats"][0][
        0
    ]["r"][0][0][0]
    df["k_rem"] = data["rem_shift_data"]["PhaseLockingData_rem"][0][0]["phasestats"][0][
        0
    ]["k"][0][0][0]
    df["p_rem"] = data["rem_shift_data"]["PhaseLockingData_rem"][0][0]["phasestats"][0][
        0
    ]["p"][0][0][0]
    df["mode_rem"] = data["rem_shift_data"]["PhaseLockingData_rem"][0][0]["phasestats"][
        0
    ][0]["mode"][0][0][0]

    # wake metrics
    df["m_wake"] = data["rem_shift_data"]["PhaseLockingData_wake"][0][0]["phasestats"][
        0
    ][0]["m"][0][0][0]
    df["r_wake"] = data["rem_shift_data"]["PhaseLockingData_wake"][0][0]["phasestats"][
        0
    ][0]["r"][0][0][0]
    df["k_wake"] = data["rem_shift_data"]["PhaseLockingData_wake"][0][0]["phasestats"][
        0
    ][0]["k"][0][0][0]
    df["p_wake"] = data["rem_shift_data"]["PhaseLockingData_wake"][0][0]["phasestats"][
        0
    ][0]["p"][0][0][0]
    df["mode_wake"] = data["rem_shift_data"]["PhaseLockingData_wake"][0][0][
        "phasestats"
    ][0][0]["mode"][0][0][0]

    def get_distros(data, state):
        return np.vstack(data["rem_shift_data"][state][0][0]["phasedistros"][0][0].T)

    def get_spikephases(data, state):
        return data["rem_shift_data"][state][0][0]["spkphases"][0][0][0]

    # add to dictionary
    data_dict = {
        "rem": {
            "phasedistros": get_distros(data, "PhaseLockingData_rem"),
            "spkphases": get_spikephases(data, "PhaseLockingData_rem"),
        },
        "wake": {
            "phasedistros": get_distros(data, "PhaseLockingData_wake"),
            "spkphases": get_spikephases(data, "PhaseLockingData_wake"),
        },
    }

    return df, data_dict


def load_SleepState_states(basepath):
    """
    loader of SleepState.states.mat

    returns dict of structures contents.

    TODO: extract more from file, this extracts the basics for now.

    """
    try:
        filename = glob.glob(os.path.join(basepath, "*.SleepState.states.mat"))[0]
    except:
        warnings.warn("file does not exist")
        return

    # load cell_metrics file
    data = sio.loadmat(filename)

    # get epoch id
    wake_id = (
        np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "WAKE")[0][0]
        + 1
    )
    rem_id = (
        np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "REM")[0][0]
        + 1
    )
    nrem_id = (
        np.where(data["SleepState"]["idx"][0][0]["statenames"][0][0][0] == "NREM")[0][0]
        + 1
    )

    # get states and timestamps vectors
    states = data["SleepState"]["idx"][0][0]["states"][0][0]
    timestamps = data["SleepState"]["idx"][0][0]["timestamps"][0][0]

    # set up dict
    dict_ = {
        "wake_id": wake_id,
        "rem_id": rem_id,
        "nrem_id": nrem_id,
        "states": states,
        "timestamps": timestamps,
    }

    # iter through states and add to dict
    dt = data["SleepState"]["ints"][0][0].dtype
    for dn in dt.names:
        dict_[dn] = data["SleepState"]["ints"][0][0][dn][0][0]

    return dict_


def load_animal_behavior(
    basepath: str, alternative_file: Union[str, None] = None
) -> pd.DataFrame:
    """
    load_animal_behavior loads basename.animal.behavior.mat files created by general_behavior_file.m
    The output is a pandas data frame with [time,x,y,z,linerized,speed,acceleration,trials,epochs]

    Ryan H 2021
    """

    if alternative_file is None:
        try:
            filename = glob.glob(os.path.join(basepath, "*.animal.behavior.mat"))[0]
        except:
            warnings.warn("file does not exist")
            return pd.DataFrame()
    else:
        try:
            filename = glob.glob(
                os.path.join(basepath, "*" + alternative_file + ".mat")
            )[0]
        except:
            warnings.warn("file does not exist")
            return pd.DataFrame()
    data = []
    data = sio.loadmat(filename, simplify_cells=True)

    df = pd.DataFrame()
    # add timestamps first which provide the correct shape of df
    # here, I'm naming them time, but this should be depreciated
    df["time"] = data["behavior"]["timestamps"]

    # add all other position coordinates to df (will add everything it can within position)
    for key in data["behavior"]["position"].keys():
        try:
            df[key] = data["behavior"]["position"][key]
        except:
            pass
    # add other fields from behavior to df (acceleration,speed,states)
    for key in data["behavior"].keys():
        try:
            df[key] = data["behavior"][key]
        except:
            pass
    # add speed and acceleration
    if "speed" not in df.columns:
        df["speed"] = functions.get_speed(df[["x", "y"]].values, df.time.values)
    if "acceleration" not in df.columns:
        df.loc[1:, "acceleration"] = np.diff(df["speed"])

    trials = data["behavior"]["trials"]
    try:
        for t in range(trials.shape[0]):
            idx = (df.time >= trials[t, 0]) & (df.time <= trials[t, 1])
            df.loc[idx, "trials"] = t
    except:
        pass

    epochs = load_epoch(basepath)
    for t in range(epochs.shape[0]):
        idx = (df.time >= epochs.startTime.iloc[t]) & (
            df.time <= epochs.stopTime.iloc[t]
        )
        df.loc[idx, "epochs"] = epochs.name.iloc[t]
        df.loc[idx, "environment"] = epochs.environment.iloc[t]
    return df


def load_epoch(basepath: str) -> pd.DataFrame:
    """
    Loads epoch info from cell explorer basename.session and stores in df
    """
    try:
        filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename, simplify_cells=True)

    try:
        epoch_df = pd.DataFrame(data["session"]["epochs"])
        epoch_df["basepath"] = basepath
        return epoch_df
    except:
        epoch_df = pd.DataFrame([data["session"]["epochs"]])
        epoch_df["basepath"] = basepath
        return epoch_df


def load_trials(basepath):
    """
    Loads trials from cell explorer basename.session.behavioralTracking and stores in df
    """
    try:
        filename = glob.glob(os.path.join(basepath, "*.animal.behavior.mat"))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename, simplify_cells=True)

    try:
        df = pd.DataFrame(data=data["behavior"]["trials"])
        df.columns = ["startTime", "stopTime"]
        df["trialsID"] = data["behavior"]["trialsID"]
        return df
    except:
        df = pd.DataFrame(data=[data["behavior"]["trials"]])
        df.columns = ["startTime", "stopTime"]
        df["trialsID"] = data["behavior"]["trialsID"]
        return df


def load_brain_regions(basepath):
    """
    Loads brain region info from cell explorer basename.session and stores in dict

    Example:
        Input:
            brainRegions = load_epoch("Z:\\Data\\GirardeauG\\Rat09\\Rat09-20140327")
            print(brainRegions.keys())
            print(brainRegions['CA1'].keys())
            print(brainRegions['CA1']['channels'])
            print(brainRegions['CA1']['electrodeGroups'])
        output:
            dict_keys(['CA1', 'Unknown', 'blv', 'bmp', 'ven'])
            dict_keys(['channels', 'electrodeGroups'])
            [145 146 147 148 149 153 155 157 150 151 154 159 156 152 158 160 137 140
            129 136 138 134 130 132 142 143 144 141 131 139 133 135]
            [17 18 19 20]
    """
    filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]

    # load file
    data = sio.loadmat(filename)
    data = data["session"]
    brainRegions = {}
    for dn in data["brainRegions"][0][0].dtype.names:
        brainRegions[dn] = {
            "channels": data["brainRegions"][0][0][dn][0][0][0][0][0][0],
            "electrodeGroups": data["brainRegions"][0][0][dn][0][0][0][0][1][0],
        }

    return brainRegions


def get_animal_id(basepath):
    """return animal ID from basepath using basename.session.mat"""
    try:
        filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)
    return data["session"][0][0]["animal"][0][0]["name"][0]


def load_basic_data(basepath):

    try:
        nChannels, fs, fs_dat, shank_to_channel = loadXML(basepath)
    except:
        fs_dat = load_extracellular_metadata(basepath).get("sr")

    ripples = load_ripples_events(basepath)
    cell_metrics, data = load_cell_metrics(basepath)

    return cell_metrics, data, ripples, fs_dat


def load_spikes(
    basepath,
    putativeCellType=[],  # restrict spikes to putativeCellType
    brainRegion=[],  # restrict spikes to brainRegion
    remove_bad_unit=True,  # true for not loading bad cells (tagged in CE)
    brain_state=[],  # restrict spikes to brainstate
    other_metric=None,  # restrict spikes to other_metric
    other_metric_value=None,  # restrict spikes to other_metric_value
    support=None,  # provide time support
):
    """
    Load specific cells' spike times
    """
    if not isinstance(putativeCellType, list):
        putativeCellType = [putativeCellType]
    if not isinstance(brainRegion, list):
        brainRegion = [brainRegion]

    # get sample rate from xml or session
    try:
        _, _, fs_dat, _ = loadXML(basepath)
    except:
        fs_dat = load_extracellular_metadata(basepath).get("sr")

    # load cell metrics and spike data
    cell_metrics, data = load_cell_metrics(basepath)

    if cell_metrics is None or data is None:
        return None, None
    
    # put spike data into array st
    st = np.array(data["spikes"], dtype=object)

    # restrict cell metrics
    if len(putativeCellType) > 0:
        restrict_idx = []
        for cell_type in putativeCellType:
            restrict_idx.append(
                cell_metrics.putativeCellType.str.contains(cell_type).values
            )
        restrict_idx = np.any(restrict_idx, axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    if len(brainRegion) > 0:
        restrict_idx = []
        for brain_region in brainRegion:
            restrict_idx.append(
                cell_metrics.brainRegion.str.contains(brain_region).values
            )
        restrict_idx = np.any(restrict_idx, axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    # restrict cell metrics by arbitrary metric
    if other_metric is not None:
        # make other_metric_value a list if not already
        if not isinstance(other_metric, list):
            other_metric = [other_metric]
        if not isinstance(other_metric_value, list):
            other_metric_value = [other_metric_value]
        # check that other_metric_value is the same length as other_metric
        if len(other_metric) != len(other_metric_value):
            raise ValueError(
                "other_metric and other_metric_value must be of same length"
            )

        restrict_idx = []
        for metric, value in zip(other_metric, other_metric_value):
            restrict_idx.append(cell_metrics[metric].str.contains(value).values)
        restrict_idx = np.any(restrict_idx, axis=0)
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    if remove_bad_unit:
        # bad units will be tagged true, so only keep false values
        restrict_idx = cell_metrics.bad_unit.values == False
        cell_metrics = cell_metrics[restrict_idx]
        st = st[restrict_idx]

    # get spike train array
    try:
        if support is not None:
            st = nel.SpikeTrainArray(timestamps=st, fs=fs_dat, support=support)
        else:
            st = nel.SpikeTrainArray(timestamps=st, fs=fs_dat)
    except:  # if only single cell... should prob just skip session
        if support is not None:
            st = nel.SpikeTrainArray(timestamps=st[0], fs=fs_dat, support=support)
        else:
            st = nel.SpikeTrainArray(timestamps=st[0], fs=fs_dat)

    if len(brain_state) > 0:
        # get brain states
        brain_states = ["WAKEstate", "NREMstate", "REMstate", "THETA", "nonTHETA"]
        if brain_state not in brain_states:
            assert print("not correct brain state. Pick one", brain_states)
        else:
            state_dict = load_SleepState_states(basepath)
            state_epoch = nel.EpochArray(state_dict[brain_state])
            st = st[state_epoch]

    return st, cell_metrics


def load_deepSuperficialfromRipple(basepath, bypass_mismatch_exception=False):
    """
    Load deepSuperficialfromRipple file created by classification_DeepSuperficial.m

    """
    # locate .mat file
    file_type = "*.deepSuperficialfromRipple.channelinfo.mat"
    filename = glob.glob(basepath + os.sep + file_type)[0]

    # load matfile
    data = sio.loadmat(filename)

    channel_df = pd.DataFrame()
    name = "deepSuperficialfromRipple"

    # sometimes more channels positons will be in deepSuperficialfromRipple than in xml
    #   this is because they used channel id as an index.
    channel_df = pd.DataFrame()
    channels = np.hstack(data[name]["channel"][0][0]) * np.nan
    shanks = np.hstack(data[name]["channel"][0][0]) * np.nan

    channels_, shanks_ = zip(
        *[
            (values[0], np.tile(shank, len(values[0])))
            for shank, values in enumerate(data[name]["ripple_channels"][0][0][0])
        ]
    )
    channel_sort_idx = np.hstack(channels_) - 1
    channels[channel_sort_idx] = np.hstack(channels_)
    shanks[channel_sort_idx] = np.hstack(shanks_) + 1

    channel_df["channel"] = channels
    channel_df.loc[
        np.arange(len(channel_sort_idx)), "channel_sort_idx"
    ] = channel_sort_idx
    channel_df["shank"] = shanks

    # add distance from pyr layer (will only be accurate if polarity rev)
    channel_df["channelDistance"] = data[name]["channelDistance"][0][0].T[0]

    # add channel class (deep or superficial)
    channelClass = []
    for item in data[name]["channelClass"][0][0]:
        try:
            channelClass.append(item[0][0])
        except:
            channelClass.append("unknown")
    channel_df["channelClass"] = channelClass

    # add if shank has polarity reversal
    for shank in channel_df.shank.unique():
        if channel_df[channel_df.shank == shank].channelClass.unique().shape[0] == 2:
            channel_df.loc[channel_df.shank == shank, "polarity_reversal"] = True
        else:
            channel_df.loc[channel_df.shank == shank, "polarity_reversal"] = False

    # add ripple and sharp wave features
    labels = ["ripple_power", "ripple_amplitude", "SWR_diff", "SWR_amplitude"]
    for label in labels:
        try:
            channel_df.loc[channel_sort_idx, label] = np.hstack(
                data[name][label][0][0][0]
            )[0]
        except:
            x = np.arange(len(channel_sort_idx)) * np.nan
            x[0 : len(np.hstack(data[name][label][0][0][0])[0])] = np.hstack(
                data[name][label][0][0][0]
            )[0]
            channel_df.loc[channel_sort_idx, label] = x

    # pull put avg ripple traces and ts
    ripple_time_axis = data[name]["ripple_time_axis"][0][0][0]
    ripple_average = np.ones([channel_df.shape[0], len(ripple_time_axis)]) * np.nan

    rip_map = []
    for ch, values in zip(channels_, data[name]["ripple_average"][0][0][0]):
        if values.shape[1] > 0:
            rip_map.append(values)
        else:
            rip_map.append(np.zeros([len(ripple_time_axis), len(ch)]) * np.nan)

    ripple_average[channel_sort_idx] = np.hstack(rip_map).T

    brainRegions = load_brain_regions(basepath)
    for key, value in brainRegions.items():
        if ("ca1" in key.lower()) | ("ca2" in key.lower()):
            for shank in value["electrodeGroups"]:
                channel_df.loc[channel_df.shank == shank, "ca1_shank"] = True

    if (ripple_average.shape[0] != channel_df.shape[0]) & (~bypass_mismatch_exception):
        raise Exception(
            "size mismatch "
            + str(np.hstack(ripple_average).shape[1])
            + " and "
            + str(channel_df.shape[0])
        )

    channel_df["basepath"] = basepath

    return channel_df, ripple_average, ripple_time_axis


def load_mua_events(basepath):
    """
    Loads the MUA data from the basepath.
    Meant to load .mat file created by find_HSE.m

    input:
        basepath: str
            The path to the folder containing the MUA data.
    output:
        mua_data: pandas.DataFrame
            The pandas.DataFrame containing the MUA data

    TODO: if none exist in basepath, create one
    """

    # locate .mat file
    try:
        filename = glob.glob(basepath + os.sep + "*mua_ca1_pyr.events.mat")[0]
    except:
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename)

    # pull out and package data
    df = pd.DataFrame()
    df["start"] = data["HSE"]["timestamps"][0][0][:, 0]
    df["stop"] = data["HSE"]["timestamps"][0][0][:, 1]
    df["peaks"] = data["HSE"]["peaks"][0][0]
    df["center"] = data["HSE"]["center"][0][0]
    df["duration"] = data["HSE"]["duration"][0][0]
    df["amplitude"] = data["HSE"]["amplitudes"][0][0]
    df["amplitudeUnits"] = data["HSE"]["amplitudeUnits"][0][0][0]
    df["detectorName"] = data["HSE"]["detectorinfo"][0][0]["detectorname"][0][0][0]

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df["basepath"] = basepath
    df["basename"] = path_components[-2]
    df["animal"] = path_components[-3]

    return df


def load_manipulation(
    basepath, struct_name=None, return_epoch_array=True, merge_gap=None
):
    """
    Loads the data from the basename.eventName.manipulations.mat file and returns a pandas dataframe

    file structure defined here:
        https://cellexplorer.org/datastructure/data-structure-and-format/#manipulations

    inputs:
        basepath: string, path to the basename.eventName.manipulations.mat file
        struct_name: string, name of the structure in the mat file to load. If None, loads all the manipulation files.
        return_epoch_array: bool, if True, returns only the epoch array
        merge_gap: int, if not None, merges the epochs that are separated by less than merge_gap (sec). return_epoch_array must be True.
    outputs:
        df: pandas dataframe, with the following columns:
            - start (float): start time of the manipulation in frames
            - stop (float): stop time of the manipulation in frames
            - peaks (float): list of the peak times of the manipulation in frames
            - center (float): center time of the manipulation in frames
            - duration (float): duration of the manipulation in frames
            - amplitude (float): amplitude of the manipulation
            - amplitudeUnit (string): unit of the amplitude
    Example:
        >> basepath = r"Z:\Data\Can\OML22\day8"
        >> df_manipulation = load_manipulation(basepath,struct_name="optoStim",return_epoch_array=False)
        >> df_manipulation.head(2)

                start	    stop	    peaks	    center	    duration amplitude amplitudeUnits
        0	8426.83650	8426.84845	8426.842475	8426.842475	0.01195	19651	pulse_respect_baseline
        1	8426.85245	8426.86745	8426.859950	8426.859950	0.01500	17516	pulse_respect_baseline

        >> basepath = r"Z:\Data\Can\OML22\day8"
        >> df_manipulation = load_manipulation(basepath,struct_name="optoStim",return_epoch_array=True)
        >> df_manipulation

        <EpochArray at 0x1faba577520: 5,774 epochs> of length 1:25:656 minutes
    """
    try:
        if struct_name is None:
            filename = glob.glob(basepath + os.sep + "*manipulation.mat")
            print(filename)
            if len(filename) > 1:
                raise ValueError(
                    "multi-file not implemented yet...than one manipulation file found"
                )
            filename = filename[0]
        else:
            filename = glob.glob(
                basepath + os.sep + "*" + struct_name + ".manipulation.mat"
            )[0]
    except:
        return None
    # load matfile
    data = sio.loadmat(filename)

    if struct_name is None:
        struct_name = list(data.keys())[-1]

    df = pd.DataFrame()
    df["start"] = data[struct_name]["timestamps"][0][0][:, 0]
    df["stop"] = data[struct_name]["timestamps"][0][0][:, 1]
    df["peaks"] = data[struct_name]["peaks"][0][0]
    df["center"] = data[struct_name]["center"][0][0]
    df["duration"] = data[struct_name]["duration"][0][0]
    df["amplitude"] = data[struct_name]["amplitude"][0][0]
    df["amplitudeUnits"] = data[struct_name]["amplitudeUnits"][0][0][0]

    # extract event label names
    eventIDlabels = []
    for name in data[struct_name]["eventIDlabels"][0][0][0]:
        eventIDlabels.append(name[0])

    # extract numeric category labels associated with label names
    eventID = np.array(data[struct_name]["eventID"][0][0]).ravel()

    # add eventIDlabels and eventID to df
    for ev_label, ev_num in zip(eventIDlabels, np.unique(eventID)):
        df.loc[eventID == ev_num, "ev_label"] = ev_label

    if return_epoch_array:
        # get session epochs to add support for epochs
        epoch_df = load_epoch(basepath)
        # get session bounds to provide support
        session_bounds = nel.EpochArray(
            [epoch_df.startTime.iloc[0], epoch_df.stopTime.iloc[-1]]
        )
        # if many types of manipulations, add them to dictinary
        if df.ev_label.unique().size > 1:
            manipulation_epoch = {}
            for label in df.ev_label.unique():
                manipulation_epoch_ = nel.EpochArray(
                    np.array(
                        [
                            df[df.ev_label == label]["start"],
                            df[df.ev_label == label]["stop"],
                        ]
                    ).T,
                    domain=session_bounds,
                )
                if merge_gap is not None:
                    manipulation_epoch_ = manipulation_epoch_.merge(gap=merge_gap)

                manipulation_epoch[label] = manipulation_epoch_
        else:
            manipulation_epoch = nel.EpochArray(
                np.array([df["start"], df["stop"]]).T, domain=session_bounds
            )
            if merge_gap is not None:
                manipulation_epoch = manipulation_epoch.merge(gap=merge_gap)

        return manipulation_epoch
    else:
        return df


def load_channel_tags(basepath):
    """
    load_channel_tags returns dictionary of tags located in basename.session.channelTags
    """
    filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]
    data = sio.loadmat(filename, simplify_cells=True)
    return data["session"]["channelTags"]


def load_extracellular_metadata(basepath):
    """
    load_extracellular returns dictionary of metadata located
        in basename.session.extracellular
    """
    filename = glob.glob(os.path.join(basepath, "*.session.mat"))[0]
    data = sio.loadmat(filename, simplify_cells=True)
    return data["session"]["extracellular"]
