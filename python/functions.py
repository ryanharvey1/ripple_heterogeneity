import numpy as np
import pandas as pd
# import hdf5storage
# import h5py
import scipy.io as sio
import glob

from scipy.signal import find_peaks
import sys,os
from sklearn.decomposition import PCA


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

def writeNeuroscopeEvents(path, ep, name):
    f = open(path, 'w')
    for i in range(len(ep)):
        f.writelines(str(ep.as_units('ms').iloc[i]['start']) + " "+name+" start "+ str(1)+"\n")
        #f.writelines(str(ep.as_units('ms').iloc[i]['peak']) + " "+name+" start "+ str(1)+"\n")
        f.writelines(str(ep.as_units('ms').iloc[i]['end']) + " "+name+" end "+ str(1)+"\n")
    f.close()
    return

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

    # load cell_metrics file
    data = sio.loadmat(filename)
    dt = data['cell_metrics'].dtype

    # construct data frame with features per neuron
    df = pd.DataFrame()
    # count units
    n_cells = data['cell_metrics']['UID'][0][0][0].size
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

    # extract other general data and put into dict    
    data_ = extract_general(data)

    return df,data_

def load_SWRunitMetrics(basepath):
    filename = glob.glob(os.path.join(basepath,'*.SWRunitMetrics.mat'))
    data = sio.loadmat(filename[0])
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

    df2 = pd.DataFrame()
    for epoch in data['SWRunitMetrics'].dtype.names:
        if data['SWRunitMetrics'][epoch][0][0].size>0: # not empty
            df_ = extract_swr_epoch_data(data,epoch)
            if df_.size>0:
                df2 = df2.append(df_,ignore_index=True)

    return df2

def load_ripples_events(filename):
    """
    load info from ripples.events.mat and store within df
    """
    df = pd.DataFrame()

    data = sio.loadmat(filename)
    dt = data['ripples'].dtype

    df['start'] = data['ripples']['timestamps'][0][0][:,0]
    df['stop'] = data['ripples']['timestamps'][0][0][:,1]
    df['peaks'] = data['ripples']['peaks'][0][0]
    df['amplitude'] = data['ripples']['amplitude'][0][0]
    df['duration'] = data['ripples']['duration'][0][0]
    df['frequency'] = data['ripples']['frequency'][0][0]

    if "detectorName" not in dt.names:
        df['detectorName'] = data['ripples']['detectorinfo'][0][0]['detectorname'][0][0][0]
    else:
        df['detectorName'] = data['ripples']['detectorName'][0][0][0]

    if "eventSpikingParameters" in dt.names:
        df['event_spk_thres'] = 1
    else:
        df['event_spk_thres'] = 0

    normalized_path = os.path.normpath(filename)
    path_components = normalized_path.split(os.sep)
    head_tail = os.path.split(normalized_path)
    df['basename'] = path_components[-2]
    df['basepath'] = head_tail[0]   
    df['animal'] = path_components[-3]

    return df

def linearize_position(x,y):
    """
    use PCA (a dimensionality reduction technique) to find
     the direction of maximal variance in our position data, 
     and we use this as our new 1D linear track axis.
    -Ryan H
    """
    # locate and remove nans (sklearn pca does not like nans)
    badidx = (np.isnan(x)) | (np.isnan(y)) 
    badidx_pos = np.where(badidx)
    goodidx_pos = np.where(~badidx)
    n = len(x)

    x = x[~badidx]
    y = y[~badidx]

    # perform pca and return the first 2 components
    pca = PCA(n_components=2)
    # transform our coords 
    linear = pca.fit_transform(np.array([x,y]).T)

    # add back nans
    x = np.zeros([n])
    x[badidx_pos] = np.nan
    x[goodidx_pos] = linear[:,0]

    y = np.zeros([n])
    y[badidx_pos] = np.nan
    y[goodidx_pos] = linear[:,1]

    # pca will center data at 0,0... adjust for this here
    x = x + np.abs(np.nanmin(x))
    y = y + np.abs(np.nanmin(y))

    return x,y


def crossCorr(t1, t2, binsize, nbins):
    ''' 
        Fast crossCorr 
        # crossCorr functions from Guillaume Viejo of Peyrache Lab
        # https://github.com/PeyracheLab/StarterPack/blob/master/python/main6_autocorr.py
    '''
    nt1 = len(t1)
    nt2 = len(t2)
    if np.floor(nbins/2)*2 == nbins:
        nbins = nbins+1

    m = -binsize*((nbins+1)/2)
    B = np.zeros(nbins)
    for j in range(nbins):
        B[j] = m+j*binsize

    w = ((nbins/2) * binsize)
    C = np.zeros(nbins)
    i2 = 1

    for i1 in range(nt1):
        lbound = t1[i1] - w
        while i2 < nt2 and t2[i2] < lbound:
            i2 = i2+1
        while i2 > 1 and t2[i2-1] > lbound:
            i2 = i2-1

        rbound = lbound
        l = i2
        for j in range(nbins):
            k = 0
            rbound = rbound+binsize
            while l < nt2 and t2[l] < rbound:
                l = l+1
                k = k+1

            C[j] += k

    # for j in range(nbins):
    # C[j] = C[j] / (nt1 * binsize)
    C = C/(nt1 * binsize)

    return C

def compute_AutoCorrs(spks, binsize = 0.001, nbins = 100):
    # First let's prepare a pandas dataframe to receive the data
    times = np.arange(0, binsize*(nbins+1), binsize) - (nbins*binsize)/2	
    autocorrs = pd.DataFrame(index = times, columns = np.arange(len(spks)))

    # Now we can iterate over the dictionnary of spikes
    for i,s in enumerate(spks):		
        # Calling the crossCorr function
        autocorrs[i] = crossCorr(s, s, binsize, nbins)

    # And don't forget to replace the 0 ms for 0
    autocorrs.loc[0] = 0.0
    return autocorrs

def BurstIndex_Royer_2012(autocorrs):
    # calc burst index from royer 2012
    # burst_idx will range from -1 to 1
    # -1 being non-bursty and 1 being bursty
    
    # peak range 2 - 9 ms
    peak = autocorrs.loc[0.002:0.009].max()
    # baseline idx 40 - 50 ms
    baseline = autocorrs.loc[0.04:0.05].mean()

    burst_idx = []
    for p, b in zip(peak, baseline):
        if p > b:
            burst_idx.append((p-b)/p)
        elif p < b:
            burst_idx.append((p-b)/b)
        else:
            burst_idx.append(np.nan)
    return burst_idx

def compress_repeated_epochs(epoch_df):
    """
    compress_repeated_epochs: Compresses epoch_df loaded by loading.load_epoch()
    If there are back to back epochs of the same name, it will combine them

    Input: epoch_df (uses: loading.load_epoch(basepath))
    Output: Compressed epoch_df

    Ryan H
    """
    match = np.zeros([epoch_df.environment.shape[0]])
    match[match == 0] = np.nan
    for i,ep in enumerate(epoch_df.environment[:-1]):
        if np.isnan(match[i]):
            # find match in current and next epoch
            if ep == epoch_df.environment.iloc[i+1]:
                match[i:i+2] = i
                # given match, see if there are more matches
                for match_i in np.arange(1,epoch_df.environment[:-1].shape[0]):
                    if i+1+match_i == epoch_df.environment.shape[0]:
                        break
                    if ep == epoch_df.environment.iloc[i+1+match_i]:

                        match[i:i+1+match_i+1] = i
                    else:
                        break

    for i in range(len(match)):
        if np.isnan(match[i]):
            match[i] = (i+1)*2000 # make nans large numbers that are unlikely to be real epoch

    # iter through each epoch indicator to get start and stop
    results = pd.DataFrame()
    no_nan_match = match[~np.isnan(match)]
    for m in pd.unique(no_nan_match):
        temp_dict = {'name': epoch_df[match==m].name.iloc[0],
                    'startTime':epoch_df[match==m].startTime.iloc[0],
                    'stopTime':epoch_df[match==m].stopTime.iloc[-1],
                    'environment':epoch_df[match==m].environment.iloc[0],
                    'behavioralParadigm':epoch_df[match==m].behavioralParadigm.iloc[0],
                    }
        temp_df = pd.DataFrame.from_dict(temp_dict,orient='index').T

        results = results.append(temp_df,ignore_index=True)
    return results    

def spatial_information(ratemap_,occupancy_):
    """
    spatial_information
    input:
        ratemap: 1 or 2d binned firing rate
        occupancy: binned occupancy same dim as ratemap
    output: 
        spatial information of a cell (bits/spike)

    See Adrien Peyrache 2008 Methods

    Ryan H
    """ 

    ratemap = ratemap_.copy()
    occupancy = occupancy_.copy()

    ratemap = ratemap.ravel()
    occupancy = occupancy.ravel()

    #  normalize to probability
    occupancy = occupancy/np.nansum(occupancy)
    f = np.nansum(occupancy*ratemap)
    ratemap = ratemap/f
    ix = ratemap != 0
    SB = (occupancy[ix]*ratemap[ix])*np.log2(ratemap[ix])
    return np.nansum(SB)
