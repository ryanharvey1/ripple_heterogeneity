import scipy.io as sio
import sys,os
import pandas as pd
import numpy as np
import glob
import warnings


def loadXML(path):
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
    if not os.path.exists(path):
        print("The path "+path+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(path)
    xmlfiles = [f for f in listdir if f.endswith('.xml')]
    if not len(xmlfiles):
        print("Folder contains no xml files; Exiting ...")
        sys.exit()
    new_path = os.path.join(path, xmlfiles[0])

    from xml.dom import minidom
    xmldoc = minidom.parse(new_path)
    nChannels = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('nChannels')[0].firstChild.data
    fs_dat = xmldoc.getElementsByTagName('acquisitionSystem')[0].getElementsByTagName('samplingRate')[0].firstChild.data
    fs = xmldoc.getElementsByTagName('fieldPotentials')[0].getElementsByTagName('lfpSamplingRate')[0].firstChild.data

    shank_to_channel = {}
    groups = xmldoc.getElementsByTagName('anatomicalDescription')[0].getElementsByTagName('channelGroups')[0].getElementsByTagName('group')
    for i in range(len(groups)):
        shank_to_channel[i] = np.sort([int(child.firstChild.data) for child in groups[i].getElementsByTagName('channel')])
    return int(nChannels), int(fs), int(fs_dat), shank_to_channel

def loadLFP(path, n_channels=90, channel=64, frequency=1250.0, precision='int16'):
    if type(channel) is not list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2
        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        interval = 1/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(os.path.dirname(os.path.abspath(path)),'lfp_ts.npy')
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)

            return data, timestep # nts.Tsd(timestep, data, time_units = 's')
        
    elif type(channel) is list:
        f = open(path, 'rb')
        startoffile = f.seek(0, 0)
        endoffile = f.seek(0, 2)
        bytes_size = 2

        n_samples = int((endoffile-startoffile)/n_channels/bytes_size)
        duration = n_samples/frequency
        f.close()
        with open(path, 'rb') as f:
            data = np.fromfile(f, np.int16).reshape((n_samples, n_channels))[:,channel]
            timestep = np.arange(0, len(data))/frequency
            # check if lfp time stamps exist
            lfp_ts_path = os.path.join(os.path.dirname(os.path.abspath(path)),'lfp_ts.npy')
            if os.path.exists(lfp_ts_path):
                timestep = np.load(lfp_ts_path).reshape(-1)
            return data,timestep # nts.TsdFrame(timestep, data, time_units = 's')


def load_position(path,fs=39.0625):
    if not os.path.exists(path):
        print("The path "+path+" doesn't exist; Exiting ...")
        sys.exit()
    listdir = os.listdir(path)
    whlfiles = [f for f in listdir if f.endswith('.whl')]
    if not len(whlfiles):
        print("Folder contains no whl files; Exiting ...")
        sys.exit()
    new_path = os.path.join(path, whlfiles[0])
    df = pd.read_csv(new_path,delimiter="\t",header=0,names=['x1','y1','x2','y2'])
    df[df==-1] = np.nan
    return df,fs


def load_cell_metrics(basepath):
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
        startTime = [ep['startTime'][0][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]
        stopTime = [ep['stopTime'][0][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]
        name = [ep['name'][0][0][0] for ep in data['cell_metrics']['general'][0][0]['epochs'][0][0][0]]

        epochs = pd.DataFrame()
        epochs['name'] = name
        epochs['startTime'] = startTime
        epochs['stopTime'] = stopTime
        return epochs

    def extract_general(data):
        # extract fr per unit with lag zero to ripple
        try:
            ripple_fr = [ev.T[0] for ev in data['cell_metrics']['events'][0][0]['ripples'][0][0][0]]
        except:
            ripple_fr = []
        # extract spikes times
        spikes = [spk.T[0] for spk in data['cell_metrics']['spikes'][0][0]['times'][0][0][0]]
        # extract epochs
        try:
            epochs = extract_epochs(data)
        except:
            epochs = []
        # extract avg waveforms (one wavefrom per channel on shank)
        try:
            waveforms = [w.T for w in data['cell_metrics']['waveforms'][0][0][0][0][0][0]]
        except:
            waveforms = [w.T for w in data['cell_metrics']['waveforms'][0][0][0]]
        # extract chanCoords
        try:
            chanCoords_x = data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['x'].T[0]
            chanCoords_y = data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['y'].T[0]
        except:
            chanCoords_x = []
            chanCoords_y = []

        # add to dictionary 
        data_ = {
            "acg_wide": data['cell_metrics']['acg'][0][0]['wide'][0][0],
            "acg_narrow": data['cell_metrics']['acg'][0][0]['narrow'][0][0],
            "acg_log10": data['cell_metrics']['acg'][0][0]['log10'][0][0],
            "ripple_fr": ripple_fr,
            "chanCoords_x": chanCoords_x,
            "chanCoords_y": chanCoords_y,
            "epochs": epochs,
            "spikes": spikes,
            "waveforms": waveforms
            }
        return data_ 

    def un_nest_df(df):
        # Un-nest some strings are nested within brackets (a better solution exists...)
        # locate and iterate objects in df
        for item in df.keys()[df.dtypes =="object"]:
            # if you can get the size of the first item with [0], it is nested
            # otherwise it fails and is not nested
            try:
                df[item][0][0].size
                # the below line is from: https://www.py4u.net/discuss/140913
                df[item] = df[item].str.get(0)
            except:
                continue
        return df

    filename = glob.glob(os.path.join(basepath,'*.cell_metrics.cellinfo.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return 

    # load cell_metrics file
    data = sio.loadmat(filename)

    # construct data frame with features per neuron
    df = pd.DataFrame()
    # count units
    n_cells = data['cell_metrics']['UID'][0][0][0].size
    dt = data['cell_metrics'].dtype
    for dn in dt.names:
        # check if var has the right n of units and is a vector
        try:
            if ((data['cell_metrics'][dn][0][0][0][0].size == 1) &
                    (data['cell_metrics'][dn][0][0][0].size == n_cells)):
                    
                df[dn] = data['cell_metrics'][dn][0][0][0]
        except:
            continue

    # add data from general metrics        
    df['basename'] = data['cell_metrics']['general'][0][0]['basename'][0][0][0]
    df['basepath'] = data['cell_metrics']['general'][0][0]['basepath'][0][0][0]
    df['sex'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['sex'][0][0][0]
    df['species'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['species'][0][0][0]
    df['strain'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['strain'][0][0][0]
    df['geneticLine'] = data['cell_metrics']['general'][0][0]['animal'][0][0]['geneticLine'][0][0][0]
    df['cellCount'] = data['cell_metrics']['general'][0][0]['cellCount'][0][0][0][0]

    # fix nesting issue for strings
    df = un_nest_df(df)

    # extract other general data and put into dict    
    data_ = extract_general(data)

    return df,data_

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

    def extract_swr_epoch_data(data,epoch):
        # get var names
        dt = data['SWRunitMetrics'][epoch][0][0].dtype

        df2 = pd.DataFrame()

        # get n units
        # there might be other fields within here like the epoch timestamps
        # skip those by returning empty df
        try:
            n_cells = data['SWRunitMetrics'][epoch][0][0][0]['particip'][0].shape[0]
        except:
            return df2

        for dn in dt.names:
            if (
                (data['SWRunitMetrics'][epoch][0][0][0][dn][0].shape[1] == 1) &
                (data['SWRunitMetrics'][epoch][0][0][0][dn][0].shape[0] == n_cells)
                ):
                df2[dn] = data['SWRunitMetrics'][epoch][0][0][0][dn][0].T[0]
        df2['epoch'] = epoch
        return df2

    filename = glob.glob(os.path.join(basepath,'*.SWRunitMetrics.mat'))[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load file
    data = sio.loadmat(filename)

    df2 = pd.DataFrame()
    # loop through each available epoch and pull out contents
    for epoch in data['SWRunitMetrics'].dtype.names:
        if data['SWRunitMetrics'][epoch][0][0].size>0: # not empty
            
            # call content extractor 
            df_ = extract_swr_epoch_data(data,epoch)

            # append conents to overall data frame
            if df_.size>0:
                df2 = df2.append(df_,ignore_index=True)

    return df2

def load_ripples_events(basepath):
    """
    load info from ripples.events.mat and store within df

    basepath: path to your session where ripples.events.mat is
    
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
    filename = glob.glob(basepath+os.sep+'*ripples.events.mat')[0]

    # check if saved file exists
    if not os.path.exists(filename):
        warnings.warn("file does not exist")
        return pd.DataFrame()

    # load matfile
    data = sio.loadmat(filename)

    # make data frame of known fields 
    df = pd.DataFrame()
    df['start'] = data['ripples']['timestamps'][0][0][:,0]
    df['stop'] = data['ripples']['timestamps'][0][0][:,1]
    df['peaks'] = data['ripples']['peaks'][0][0]
    df['amplitude'] = data['ripples']['amplitude'][0][0]
    df['duration'] = data['ripples']['duration'][0][0]
    df['frequency'] = data['ripples']['frequency'][0][0]
    
    try:
        df['detectorName'] = data['ripples']['detectorinfo'][0][0]['detectorname'][0][0][0]
    except:
        df['detectorName'] = data['ripples']['detectorName'][0][0][0]

    dt = data['ripples'].dtype
    if "eventSpikingParameters" in dt.names:
        df['event_spk_thres'] = 1
    else:
        df['event_spk_thres'] = 0

    # get basename and animal
    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    df['basepath'] = basepath  
    df['basename'] = path_components[-2]
    df['animal'] = path_components[-3]

    return df