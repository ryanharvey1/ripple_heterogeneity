import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from numba import jit
from numba import int8
from numba.typed import List
from numba import typeof
import sys
import nelpy as nel

def set_plotting_defaults():
    tex_fonts = {
    #     # Use LaTeX to write all text
        "font.family": "serif",
        # Use 10pt font in plots
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "svg.fonttype": 'none'
    }
    plt.style.use('seaborn-paper')
    plt.rcParams.update(tex_fonts)


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

def writeNeuroscopeEvents(path, ep, name):
    f = open(path, 'w')
    for i in range(len(ep)):
        f.writelines(str(ep.as_units('ms').iloc[i]['start']) + " "+name+" start "+ str(1)+"\n")
        #f.writelines(str(ep.as_units('ms').iloc[i]['peak']) + " "+name+" start "+ str(1)+"\n")
        f.writelines(str(ep.as_units('ms').iloc[i]['end']) + " "+name+" end "+ str(1)+"\n")
    f.close()
    return

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

@jit(nopython=True)
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

def compute_psth(spikes,event,bin_width=0.002,n_bins=100):

    times = np.arange(0, bin_width*(n_bins+1), bin_width) - (n_bins*bin_width)/2	
    ccg = pd.DataFrame(index = times, columns = np.arange(len(spikes)))
    # Now we can iterate over spikes
    for i,s in enumerate(spikes):		
        ccg[i] = crossCorr(event, s, bin_width, n_bins)
    return ccg

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

def Otsu(vector):
    """
    The Otsu method for splitting data into two groups.
    This is somewhat equivalent to kmeans(vector,2), but while the kmeans implementation
    finds a local minimum and may therefore produce different results each time,
    the Otsu implementation is guaranteed to find the best division every time.

    input: 
        vector: arbitrary vector
    output:
        group: binary class
        threshold: threshold used for classification
        em: effectiveness metric

    From Raly
    """
    sorted = np.sort(vector)
    n = len(vector)
    intraClassVariance = [np.nan]*n
    for i in np.arange(n):
        p = (i+1)/n
        p0 = 1 - p
        intraClassVariance[i] = p*np.var(sorted[0:i+1]) + p0*np.var(sorted[i+1:])
    
    minIntraVariance = np.nanmin(intraClassVariance)
    idx = np.nanargmin(intraClassVariance)
    threshold = sorted[idx]
    group = (vector > threshold)
    
    em = 1 - (minIntraVariance/np.var(vector))

    return group,threshold,em

def find_sig_assemblies(patterns):
    is_member = []
    keep_assembly = []
    for pat in patterns:
        isMember,_,_ = Otsu(np.abs(pat))
        is_member.append(isMember)

        if any(pat[isMember] < 0) & any(pat[isMember] > 0):
            keep_assembly.append(False)
        elif sum(isMember) == 0:
            keep_assembly.append(False)
        else:
            keep_assembly.append(True)

    is_member = np.array(is_member)

    return patterns[keep_assembly],is_member[keep_assembly],keep_assembly,is_member

@jit(nopython=True)
def get_participation_(st,event_starts,event_stops):
    unit_mat = np.zeros((len(st),len(event_starts)),dtype=int8)
    for rip in range(len(event_starts)):
        for i,s in enumerate(st):
            idx = (s >= event_starts[rip]) & (s <= event_stops[rip])
            unit_mat[i,rip] = sum(idx) > 0
    return unit_mat 

def fix_array_or_list(l, dtype=None):
    if isinstance(l, np.ndarray):
        return l
    if not l:
        raise ValueError("List must contain at least one element")
    if dtype is None:
        dt = typeof(l[0])
    else:
        dt = dtype
    new_list = List.empty_list(dt)
    for x in l:
        new_list.append(x)
    return new_list

def get_participation(st,event_starts,event_stops):
    """
    get participation prob.
    make matrix n rows (units) by n cols (ripple epochs)
    Input:
        st: spike train nelpy object that is epoched by ripples
        ripple_epochs: ripple events in nelpy epoch object 
    """
    return get_participation_(fix_array_or_list(list(st)),event_starts,event_stops)

def get_significant_events(scores, shuffled_scores, q=95):
    """Return the significant events based on percentiles.
    NOTE: The score is compared to the distribution of scores obtained
    using the randomized data and a Monte Carlo p-value can be computed
    according to: p = (r+1)/(n+1), where r is the number of
    randomizations resulting in a score higher than (ETIENNE EDIT: OR EQUAL TO?)
    the real score and n is the total number of randomizations performed.
    Parameters
    ----------
    scores : array of shape (n_events,)
    shuffled_scores : array of shape (n_shuffles, n_events)
    q : float in range of [0,100]
        Percentile to compute, which must be between 0 and 100 inclusive.
    Returns
    -------
    sig_event_idx : array of shape (n_sig_events,)
        Indices (from 0 to n_events-1) of significant events.
    pvalues :
    """
    # check shape and correct if needed
    if shuffled_scores.shape[1] != len(scores):
        shuffled_scores = shuffled_scores.T

    n, _ = shuffled_scores.shape
    r = np.sum(abs(shuffled_scores) >= abs(scores), axis=0)
    pvalues = (r+1)/(n+1)

    # set nan scores to 1
    pvalues[np.isnan(scores)] = 1
    
    sig_event_idx = np.argwhere(scores > np.percentile(
        shuffled_scores,
        axis=0,
        q=q)).squeeze()

    return np.atleast_1d(sig_event_idx), np.atleast_1d(pvalues)

def find_laps(Vts,Vdata,newLapThreshold=15,good_laps=True,
                edgethresh=0.1,completeprop=0.2,posbins=50):
    """
    Find Laps in linear track 

    INPUT:
    Vts: timestamps
    Vdata: x coords 

    newLapThreshold: endpoint proximity threshold in percent of track length (default = 15%);
                    whenever rat enters the proximity zone of e.g. 15% of tracklength near a end, a new lap
                    is started and the maximum (or minimum) is searched
                    for a Lap-Top  or Lap-Bottom (around 0 end).

    good_laps: run find_good_laps to remove laps with excess nans and 
                parts of laps where rat turns around in middle of track

    OUTPUT:
    laps  .... 1*nLaps struct array with fields
    laps(i).start_ts  ... start timestamp of i-th lap
    laps(i).pos       ... the value of input position V at lap start point
    laps(i).start_idx ... the index of the new lap start frame in input V 
    laps(i).direction ... +1/-1 for up/down laps 

    From NSMA toolbox
    Author: PL
    Version: 0.9  05/12/2005
    edited by Ryan Harvey to work with standard linear track
    edited for use in python by Ryan h 2022
    """
    
    TL = np.abs(np.nanmax(Vdata) - np.nanmin(Vdata))   #% track length
    th1= np.nanmin(Vdata) + TL*newLapThreshold / 100 #        % lower threshold for lower end
    th2 = np.nanmax(Vdata) - TL*newLapThreshold / 100 #      % upper threshold for upper end

    # % loop over all frames
    laps = pd.DataFrame()
    laps.loc[0,'start_ts'] = Vts[0]
    laps.loc[0,'pos'] = Vdata[0]
    laps.loc[0,'start_idx'] = 1
    laps.loc[0,'direction'] = 0
    iLap = 0

    newUpThCross = 1     #% flag for new lap top search
    newDownThCross = 1    # % flag for new lap top search
    for i in range(len(Vdata)):
        if Vdata[i] < th1:   # % search for min
            if newUpThCross == 1:  #     % start a new lap
                newUpThCross = 0
                newDownThCross = 1
                iLap = iLap + 1
                laps.loc[iLap,'start_ts'] = Vts[i]
                laps.loc[iLap,'pos'] = Vdata[i]
                laps.loc[iLap,'start_idx'] = i
                laps.loc[iLap,'direction'] = 1
            
            if Vdata[i] < laps.iloc[iLap].pos:     # % record new min if any
                laps.loc[iLap,'start_ts'] = Vts[i]
                laps.loc[iLap,'pos'] = Vdata[i]
                laps.loc[iLap,'start_idx'] = i
        
        if Vdata[i] > th2:  # % search for max
            if newDownThCross: #      % start a new lap
                newUpThCross = 1
                newDownThCross = 0
                iLap = iLap + 1
                laps.loc[iLap,'start_ts'] = Vts[i]
                laps.loc[iLap,'pos'] = Vdata[i]
                laps.loc[iLap,'start_idx'] = i
                laps.loc[iLap,'direction'] = -1
            
            if Vdata[i] > laps.iloc[iLap].pos:     #  % record new min if any
                laps.loc[iLap,'start_ts'] = Vts[i]
                laps.loc[iLap,'pos'] = Vdata[i]
                laps.loc[iLap,'start_idx'] = i

    # % fix direction of first lap which was unknown above
    # % make first lap direction opposite of second lap's direction (laps alternate!)
    laps.iloc[0].direction = -laps.iloc[1].direction

    # % make sure laps cross the halfway point
    middle = np.nanmedian(np.arange(np.nanmin(Vdata),np.nanmax(Vdata)))
    i = 0
    while True:
        try:
            positions = np.arange(laps.iloc[i].pos,laps.iloc[i+1].pos)
        except:
            positions = [np.nan,np.nan]
        if (any(positions > middle) == True) & (any(positions < middle) == False):
            laps = laps.drop(laps.index[i+1])
        i = i+1
        if i+1 >= len(laps.pos):
            if len(laps.pos) < iLap:
                laps.iloc[0].direction = -laps.iloc[1].direction
            break

    if good_laps:
        laps = find_good_laps(Vts,
                            Vdata,
                            laps,
                            edgethresh=edgethresh,
                            completeprop=completeprop,
                            posbins=posbins)

    return laps

def peakdetz(v, delta, lookformax=1, backwards=0):
    """
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDETZ(V, DELTA, lookformax, backwards) finds 
    %        the local maxima and minima ("peaks") in the vector V.
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA. MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    %
    % ZN edit 04/2010: added option to specify looking for troughs or peaks
    % first (lookformax variable: if 1, will look for peaks first, if 0 will
    % look for troughs; default is look for peaks); and option to go backwards
    % (so that find last instance of a peak/trough value instead of the first
    % instance: backwards variable: if 1 will go backwards, if 0 or absent, 
    % will go forwards); and changed it so that last min/max value will be 
    % assigned

    edited for use in python by Ryan H 2022
    """

    maxtab = []
    mintab = []

    v = np.asarray(v)

    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    if backwards == 0:
        inc = 1
        first = 0
        last = len(v)
    elif backwards:
        inc = -1
        first = len(v)
        last = 0

    mn = np.inf
    mx = -np.inf
    mnpos = np.nan
    mxpos = np.nan

    for ii in np.arange(first,last,inc):
        this = v[ii]
        if this > mx:
            mx = this
            mxpos = ii
        if this < mn:
            mn = this
            mnpos = ii
    
        if lookformax:
            try:
                idx = mx-delta>mintab[-1]
            except:
                idx = mx-delta>mintab

            if (this < mx-delta) | ((ii==last) & (len(mintab)>0) & idx):
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = ii
                lookformax = 0
        else:
            try:
                idx = mx-delta>maxtab[-1]
            except:
                idx = mx-delta>maxtab
            if (this > mn+delta) | ((ii==last) & (len(maxtab)>0) & idx):
                mintab.append((mnpos, mn))
                mx = this
                mxpos = ii
                lookformax = 1

    if (len(maxtab)==0) & (len(mintab)==0):
        if lookformax:
            if mx-mn>delta:
                maxtab = [mxpos, mx]
        else:
            if mx-mn>delta:
                mintab = [mnpos, mn]
    return maxtab, mintab

def find_good_laps(ts,V_rest,laps,edgethresh=0.1,completeprop=0.2,posbins=50):
    """
    % [startgoodlaps, stopgoodlaps, laps] = 
    %        find_good_laps(V_rest,laps,edgethresh,completeprop,posbins)
    %
    % find and eliminate laps which have too many NaNs (because rat was off 
    % track), and parts of laps where rat turns around in middle of track
    % 
    % inputs: V_rest: V coordinates of rat with off track periods masked out 
    %                 (as NaNs)
    %         laps: struct with lap start and end times (generated by
    %               find_laps)
    %         edgethresh: threshold for detection of a turn around point 
    %                     (proportion of length of track) (default = 0.1)
    %         completeprop: the amount of lap that can be missing (NaNs) to
    %                       still be considered a lap (default = 0.2).
    %         plotlaps: flag for making plots of each lap, and pause for user
    %                   to hit key to continue (default = 1)
    %         posbins: number of bins to divide the track into to determine
    %                  position coverage percentage; at 60frames/s want at 
    %                  least 2cm/bin (default = 50bins; this works for 100+ cm 
    %                  track, as long as V_rest is in cm)
    % outputs: 
    %          laps: a new laps struct, with the bad laps removed
    %
    % ZN 04/2011
    Edited for use in python by Ryan H 2022
    """

    if edgethresh > 1:    # % in case edgethresh is input as a percentage instead of a proportion
        edgethresh = edgethresh / 100

    if completeprop > 1:     # % in case completeprop is input as a percentage instead of a proportion
        completeprop = completeprop / 100

    bottomend = np.nanmin(V_rest)
    topend = np.nanmax(V_rest)
    bins = np.arange(bottomend,topend+(topend-bottomend)/posbins,(topend-bottomend)/posbins)
    delta = (topend - bottomend)*edgethresh   #  % threshold for peak/trough detection
    startgoodlaps = []
    stopgoodlaps = []

    l = 0
    while l < len(laps)-1:
        #% select out just this lap
        if l == len(laps):
            endoflap = ts[-1]
        else:
            endoflap = laps.iloc[l+1].start_ts
        
        v = V_rest[np.where(ts == laps.iloc[l].start_ts)[0][0]:np.where(ts == endoflap)[0][0]]
        t = ts[np.where(ts==laps.iloc[l].start_ts)[0][0]:np.where(ts==endoflap)[0][0]]

        #% find turn around points during this lap
        lookformax = laps.iloc[l].direction == 1
        peak,trough = peakdetz(v, delta, lookformax, 0)
        
        if lookformax:
        #% find the direct path from bottomend to topend (or mark lap for
        #% deleting if the turn around points are not in those ranges)
            if len(trough) > 0:
                #% find the last trough in range of bottomend (start of lap)
                gt = len(trough)
                while ((gt>0) & (trough(gt,2) >= 2*delta + bottomend)):
                    gt=gt-1
                
                #% assign the next peak after that trough as the end of the lap
                #% (or mark lap for deleting, if that peak is not at topend)
                if gt == 0:
                    if (peak[1,2] > topend-2*delta):
                        t = t[0:peak[0]]
                        v = v[0:peak[0]]
                    else:
                        #% this marks the lap for deleting
                        t = t[0:5]
                        v = v[0:5]
                else:
                    et = len(peak)
                    if (gt+1 > et):
                        gt = 0
                        t=t[0:2]
                        v=v[0:2]
                    else:
                        t=t[trough[gt,1]:peak[gt+1,1]]
                        v=v[trough[gt,1]:peak[gt+1,1]]
                
            else:
            #% make sure peak exists and is in range of topend
                if len(peak) == 0:
                    if len(t) > 2:
                        t=t[0:2]
                        v=v[0:2]
                elif peak[1] < topend-2*delta:
                    #% this marks the lap for deleting
                    if len(t) > 5:
                        t=t[0:5]
                        v=v[0:5]
                 
        else: # % if lookformax
        #% find the direct path from topend to bottomend (or mark lap for
        #% deleting if the turn around points are not in those ranges)
            if len(peak) > 0:
            #% find the last peak in range of topend (start of lap)
                gt = len(peak)
                while (gt>0) & (peak[gt,2]<=topend-2*delta):
                    gt=gt-1
                # % assign the next trough after that peak as the end of the lap
                # % (or mark lap for deleting, if that trough is not at bottomend)
                if gt==0:
                    if trough(1,2)<bottomend+2*delta:
                        t=t[1:trough[0]]
                        v=v[1:trough[0]]
                    else:
                        # % this marks the lap for deleting
                        t=t[0:5]
                        v=v[0:5]
                else:
                    et = len(trough)
                    if gt+1 > et:
                        t=t[0:2]
                        v=v[0:2]
                        gt=0
                    else:
                        t=t[peak[gt,1]:trough[gt+1,1]]
                        v=v[peak[gt,1]:trough[gt+1,1]]
             
            else: #% if ~isempty(peak)
                #% make sure trough exists and is in range of bottomend
                if len(trough) == 0:
                 
                    if len(t)>2:
                        t = t[0:2]
                        v = v[0:2]
                    
                elif trough[1] > bottomend+2*delta:
                    #% this marks the lap for deleting
                    if len(t) > 5:
                        t = t[0:5]
                        v = v[0:5]
                     
        vcovered,_ = np.histogram(v,bins=bins)

        if len(v) < 3:
        #% eliminate the lap if it is non-existent (as is sometimes the case for lap 1)
            laps = laps.drop(laps.index[l])
         
        #% eliminate the lap if >completeprop of it is NaNs or if it has been marked for
        #% deleting above
        elif (len(v) < 6) | (sum(vcovered==0)>completeprop*posbins):
            laps.drop(laps.index[l])
            #% remove the other lap from the lap pair
            if l % 2 == 0:
                #% delete previous lap from laps
                laps = laps.drop(laps.index[l-1])
             
                #% change goodlaps markers to delete previous lap from laps
                if len(stopgoodlaps) > 0:
                    if 'lastlapend' not in locals() | (startgoodlaps[-1] > lastlapend):
                        startgoodlaps[-1] = []
                        stopgoodlaps[-1] = []
                    else:
                        stopgoodlaps[-1] = lastlapend
                    
                l=l-1
            elif l<=len(laps) & l > 1:
                #% delete next lap from laps
                laps = laps.drop(laps.index[l])
              
        else: #% if lap is good
            #% store last lap end just in case have to delete this lap with next lap
            if len(stopgoodlaps) > 0:
                lastlapend = stopgoodlaps[-1]
            
            # % add this lap to goodlaps
            try:
                idx = stopgoodlaps[-1]==t[0]
            except:
                idx = stopgoodlaps == t[0]
            if (len(stopgoodlaps)>0) & (idx):
                stopgoodlaps[-1] = t[-1]
            else:
                startgoodlaps.append(t[0])
                stopgoodlaps.append(t[-1])
            
            l = l+1
  
    return laps

def get_linear_track_lap_epochs(ts,x,newLapThreshold=15,
                                good_laps=True,edgethresh=0.1,
                                completeprop=0.2,posbins=50):

    """
    get_linear_track_lap_epochs: def that calls find_laps and outputs nelpy epochs
        for out and inbound running directions
    """
    laps = find_laps(np.array(ts),
                    np.array(x),
                    newLapThreshold=newLapThreshold,
                    good_laps=good_laps,
                    edgethresh=edgethresh,
                    completeprop=completeprop,
                    posbins=posbins)    

    outbound_start = []
    outbound_stop = []
    inbound_start = []
    inbound_stop = []

    for i in range(len(laps)-1):
        if laps.iloc[i].direction == 1:
            outbound_start.append(laps.iloc[i].start_ts)
            outbound_stop.append(laps.iloc[i+1].start_ts)

        if laps.iloc[i].direction == -1:
            inbound_start.append(laps.iloc[i].start_ts)
            inbound_stop.append(laps.iloc[i+1].start_ts)

    outbound_epochs = nel.EpochArray([np.array([outbound_start,outbound_stop]).T])
    inbound_epochs = nel.EpochArray([np.array([inbound_start,inbound_stop]).T])
    
    return outbound_epochs,inbound_epochs