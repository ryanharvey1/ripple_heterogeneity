import numpy as np
import pandas as pd
# import hdf5storage
# import h5py
import scipy.io as sio

from scipy.signal import find_peaks
import os
import sys

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.
    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

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
    return int(nChannels), int(fs), shank_to_channel    

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

def writeNeuroscopeEvents(path, ep, name):
    f = open(path, 'w')
    for i in range(len(ep)):
        f.writelines(str(ep.as_units('ms').iloc[i]['start']) + " "+name+" start "+ str(1)+"\n")
        #f.writelines(str(ep.as_units('ms').iloc[i]['peak']) + " "+name+" start "+ str(1)+"\n")
        f.writelines(str(ep.as_units('ms').iloc[i]['end']) + " "+name+" end "+ str(1)+"\n")
    f.close()
    return


def load_cell_metrics(filename):
    """ 
    loader of cell-explorer cell_metrics.cellinfo.mat

    Inputs: filename: path to cell_metrics.cellinfo.mat
    outputs: df: data frame of single unit features
    data_: dict with data that does not fit nicely into a dataframe (waveforms, acgs, epochs, etc.)

    TODO: extract all fields from cell_metrics.cellinfo. 

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
        ripple_fr = [ev.T[0] for ev in data['cell_metrics']['events'][0][0]['ripples'][0][0][0]]
        # extract spikes times
        spikes = [spk.T[0] for spk in data['cell_metrics']['spikes'][0][0]['times'][0][0][0]]
        # extract epochs
        epochs = extract_epochs(data)
        # extract avg waveforms (one wavefrom per channel on shank)
        waveforms = [w.T for w in data['cell_metrics']['waveforms'][0][0][0][0][0][0]]
        # add to dictionary 
        data_ = {
            "acg_wide": data['cell_metrics']['acg'][0][0]['wide'][0][0],
            "acg_narrow": data['cell_metrics']['acg'][0][0]['narrow'][0][0],
            "acg_log10": data['cell_metrics']['acg'][0][0]['log10'][0][0],
            "ripple_fr": ripple_fr,
            "chanCoords_x": data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['x'].T[0],
            "chanCoords_y": data['cell_metrics']['general'][0][0]['chanCoords'][0][0][0][0]['y'].T[0],
            "epochs": epochs,
            "spikes": spikes,
            "waveforms": waveforms
            }
        return data_ 

    # load cell_metrics file
    data = sio.loadmat(filename)
    dt = data['cell_metrics'].dtype

    # construct data frame with features per neuron
    df = pd.DataFrame()
    for dn in dt.names:
        try:
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

    # extract other general data and put into dict    
    data_ = extract_general(data)

    return df,data_