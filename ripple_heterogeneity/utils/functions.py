import itertools
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
import warnings
from scipy import stats
from ripple_heterogeneity.assembly import find_sig_assembly
from itertools import combinations
from scipy import signal

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
        "svg.fonttype": "none",
    }
    plt.style.use("seaborn-paper")
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
    if width == "thesis":
        width_pt = 426.79135
    elif width == "beamer":
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)


def writeNeuroscopeEvents(path, ep, name):
    f = open(path, "w")
    for i in range(len(ep)):
        f.writelines(
            str(ep.as_units("ms").iloc[i]["start"])
            + " "
            + name
            + " start "
            + str(1)
            + "\n"
        )
        # f.writelines(str(ep.as_units('ms').iloc[i]['peak']) + " "+name+" start "+ str(1)+"\n")
        f.writelines(
            str(ep.as_units("ms").iloc[i]["end"]) + " " + name + " end " + str(1) + "\n"
        )
    f.close()
    return


def linearize_position(x, y):
    """
    use PCA (a dimensionality reduction technique) to find
    the direction of maximal variance in our position data,
    and we use this as our new 1D linear track axis.
    
    Input:
        x: numpy array of shape (n,1)
        y: numpy array of shape (n,1)
    Output:
        x_lin: numpy array of shape (n,1)
        y_lin: numpy array of shape (n,1)

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
    linear = pca.fit_transform(np.array([x, y]).T)

    # add back nans
    x = np.zeros([n])
    x[badidx_pos] = np.nan
    x[goodidx_pos] = linear[:, 0]

    y = np.zeros([n])
    y[badidx_pos] = np.nan
    y[goodidx_pos] = linear[:, 1]

    # pca will center data at 0,0... adjust for this here
    x = x + np.abs(np.nanmin(x))
    y = y + np.abs(np.nanmin(y))

    return x, y


@jit(nopython=True)
def crossCorr(t1, t2, binsize, nbins):
    """
    Fast crossCorr
    # crossCorr functions from Guillaume Viejo of Peyrache Lab
    # https://github.com/PeyracheLab/StarterPack/blob/master/python/main6_autocorr.py
    """
    nt1 = len(t1)
    nt2 = len(t2)
    if np.floor(nbins / 2) * 2 == nbins:
        nbins = nbins + 1

    m = -binsize * ((nbins + 1) / 2)
    B = np.zeros(nbins)
    for j in range(nbins):
        B[j] = m + j * binsize

    w = (nbins / 2) * binsize
    C = np.zeros(nbins)
    i2 = 1

    for i1 in range(nt1):
        lbound = t1[i1] - w
        while i2 < nt2 and t2[i2] < lbound:
            i2 = i2 + 1
        while i2 > 1 and t2[i2 - 1] > lbound:
            i2 = i2 - 1

        rbound = lbound
        l = i2
        for j in range(nbins):
            k = 0
            rbound = rbound + binsize
            while l < nt2 and t2[l] < rbound:
                l = l + 1
                k = k + 1

            C[j] += k

    # for j in range(nbins):
    # C[j] = C[j] / (nt1 * binsize)
    C = C / (nt1 * binsize)

    return C


def compute_AutoCorrs(spks, binsize=0.001, nbins=100):
    # First let's prepare a pandas dataframe to receive the data
    times = np.arange(0, binsize * (nbins + 1), binsize) - (nbins * binsize) / 2
    autocorrs = pd.DataFrame(index=times, columns=np.arange(len(spks)))

    # Now we can iterate over the dictionnary of spikes
    for i, s in enumerate(spks):
        # Calling the crossCorr function
        autocorrs[i] = crossCorr(s, s, binsize, nbins)

    # And don't forget to replace the 0 ms for 0
    autocorrs.loc[0] = 0.0
    return autocorrs

def pairwise_corr(X,method='pearson',pairs=None):
    """
    Compute pairwise correlations between all rows of matrix
    Input: 
        X: numpy array of shape (n,p)
    Output:
        corr: numpy array rho
        pval: numpy array pval
        c: numpy array ref and target from which the correlation was computed
    """
    if pairs is None:
        x = np.arange(0, X.shape[0])
        pairs = np.array(list(itertools.combinations(x, 2)))

    rho = []
    pval = []
    for i, s in enumerate(pairs):
        if method == 'pearson':
            rho_, pval_ = stats.pearsonr(X[s[0], :], X[s[1], :])
        elif method == 'spearman':
            rho_, pval_ = stats.spearmanr(X[s[0], :], X[s[1], :])
        elif method == 'kendall':
            rho_, pval_ = stats.kendalltau(X[s[0], :], X[s[1], :])
        else:
            raise ValueError('method must be pearson, spearman or kendall')
        rho.append(rho_)
        pval.append(pval_)
    return rho, pval, pairs

def pairwise_cross_corr(spks, binsize=0.001, nbins=100, return_index=False, pairs=None):
    """
    Compute pairwise time-lagged correlations between cells
    Input:
        spks: list of numpy arrays of shape (n,)
        binsize: float, size of bins in seconds
        nbins: int, number of bins
        return_index: bool, return the index of the cells used for the correlation
        pairs: list of pairs of cells to compute the correlation
    Output:
        crosscorrs: pandas dataframe of shape (t,n pairs)
    
    """
    # Get unique combo without repeats
    if pairs is None:
        x = np.arange(0, spks.shape[0])
        pairs = np.array(list(itertools.combinations(x, 2)))

    # prepare a pandas dataframe to receive the data
    times = np.linspace(-(nbins*binsize)/2,(nbins*binsize)/2,nbins+1)

    crosscorrs = pd.DataFrame(index=times, columns=np.arange(len(pairs)))

    # Now we can iterate over spikes
    for i, s in enumerate(pairs):
        # Calling the crossCorr function
        crosscorrs[i] = crossCorr(spks[s[0]], spks[s[1]], binsize, nbins)

    if return_index:
        return crosscorrs, pairs
    else:
        return crosscorrs
        
def pairwise_spatial_corr(X, return_index=False, pairs=None):
    """
    Compute pairwise spatial correlations between cells
    Input:
        X: numpy array of shape (n_cells, n_space, n_space)
        return_index: bool, return the index of the cells used for the correlation
        pairs: list of pairs of cells to compute the correlation
    Output:
        spatial_corr: the pearson correlation between the cells in pairs
        pairs: list of pairs of cells used for the correlation
    
    """
    # Get unique combo without repeats
    if pairs is None:
        x = np.arange(0, X.shape[0])
        pairs = np.array(list(combinations(x, 2)))

    X[np.isnan(X)] = 0

    spatial_corr = []
    # Now we can iterate over spikes
    for i, s in enumerate(pairs):
        # Calling the crossCorr function
        spatial_corr.append(np.corrcoef(X[s[0],:,:].flatten(),X[s[1],:,:].flatten())[0,1])

    if return_index:
        return np.array(spatial_corr), pairs
    else:
        return np.array(spatial_corr)

def compute_psth(spikes, event, bin_width=0.002, n_bins=100):

    # times = np.arange(0, bin_width * (n_bins + 1), bin_width) - (n_bins * bin_width) / 2
    times = np.linspace(-(n_bins*bin_width)/2,(n_bins*bin_width)/2,n_bins+1)
    ccg = pd.DataFrame(index=times, columns=np.arange(len(spikes)))
    # Now we can iterate over spikes
    for i, s in enumerate(spikes):
        ccg[i] = crossCorr(event, s, bin_width, n_bins)
    return ccg

def compute_cross_correlogram(X, dt=1, window=0.5):
    """
    Cross-correlate two N-dimensional arrays (pairwise).
    Input:
        X: N-dimensional array of shape  (n_signals, n_timepoints)
        dt: time step in seconds, default 1 is nlags
        window: window size in seconds, output will be +- window
    Output:
        cross_correlogram: pandas dataframe with pairwise cross-correlogram
    """

    crosscorrs = {}
    pairs = list(combinations(np.arange(X.shape[0]), 2))
    for i, j in pairs:
        auc = signal.correlate(X[i], X[j])
        times = signal.correlation_lags(len(X[i]), len(X[j])) * dt
        # normalize by coeff
        normalizer = np.sqrt((X[i] ** 2).sum(axis=0) * (X[j] ** 2).sum(axis=0))
        auc /= normalizer

        crosscorrs[(i, j)] = pd.Series(index=times, data=auc, dtype="float32")
    crosscorrs = pd.DataFrame.from_dict(crosscorrs)

    if window is None:
        return crosscorrs
    else:
        return crosscorrs[(crosscorrs.index >= -window) & (crosscorrs.index <= window)]

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
            burst_idx.append((p - b) / p)
        elif p < b:
            burst_idx.append((p - b) / b)
        else:
            burst_idx.append(np.nan)
    return burst_idx


def select_burst_spikes(spikes, mode="bursts", isiBursts=0.006, isiSpikes=0.020):
    """
    select_burst_spikes - Discriminate bursts vs single spikes.
    adpated from: http://fmatoolbox.sourceforge.net/Contents/FMAToolbox/Analyses/SelectSpikes.html

    Input:
        spikes: list of spike times
        mode: either 'bursts' (default) or 'single'
        isiBursts: max inter-spike interval for bursts (default = 0.006)
        isiSpikes: min for single spikes (default = 0.020)
    Output:
        selected: a logical vector indicating for each spike whether it
                    matches the criterion
    """

    dt = np.diff(spikes)

    if mode == "bursts":
        b = dt < isiBursts
        # either next or previous isi < threshold
        selected = np.insert(b, 0, False, axis=0) | np.append(b, False)
    else:
        s = dt > isiSpikes
        # either next or previous isi > threshold
        selected = np.insert(s, 0, False, axis=0) & np.append(s, False)

    return selected


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
    for i, ep in enumerate(epoch_df.environment[:-1]):
        if np.isnan(match[i]):
            # find match in current and next epoch
            if ep == epoch_df.environment.iloc[i + 1]:
                match[i : i + 2] = i
                # given match, see if there are more matches
                for match_i in np.arange(1, epoch_df.environment[:-1].shape[0]):
                    if i + 1 + match_i == epoch_df.environment.shape[0]:
                        break
                    if ep == epoch_df.environment.iloc[i + 1 + match_i]:

                        match[i : i + 1 + match_i + 1] = i
                    else:
                        break

    for i in range(len(match)):
        if np.isnan(match[i]):
            # make nans large numbers that are unlikely to be real epoch
            match[i] = (i + 1) * 2000

    # iter through each epoch indicator to get start and stop
    results = pd.DataFrame()
    no_nan_match = match[~np.isnan(match)]
    for m in pd.unique(no_nan_match):
        temp_dict = {}
        for item in epoch_df.keys():
            temp_dict[item] = epoch_df[match == m][item].iloc[0]
            
        temp_dict["startTime"] = epoch_df[match == m].startTime.min()
        temp_dict["stopTime"] = epoch_df[match == m].stopTime.max()

        temp_df = pd.DataFrame.from_dict(temp_dict, orient="index").T

        results = pd.concat([results,temp_df], ignore_index=True)
    return results


def spatial_information(ratemap_, occupancy_):
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
    occupancy = occupancy / np.nansum(occupancy)
    f = np.nansum(occupancy * ratemap)
    ratemap = ratemap / f
    ix = ratemap != 0
    SB = (occupancy[ix] * ratemap[ix]) * np.log2(ratemap[ix])
    return np.nansum(SB)


def find_sig_assemblies(patterns):
    """
    Input:
        patterns: a list of patterns
    Output:
        patterns[keep_assembly]: a list of patterns that are significant
        is_member[keep_assembly]: a list of booleans indicating whether each pattern is a significant assembly
        keep_assembly: a list of indices of the significant assemblies
        is_member: a list of booleans indicating whether each unit is a significant member of an assembly
    """

    print("functions.find_sig_assemblies: is deprecated, use find_sig_assembly.main instead")
    patterns, is_member_sig, keep_assembly, is_member = find_sig_assembly.main(patterns)
    return patterns, is_member_sig, keep_assembly, is_member


@jit(nopython=True)
def get_participation_(st, event_starts, event_stops):
    """
    get_participation_: Calculates the participation of each unit to each event
    Input:
        st: spike times
        event_starts: start times of events
        event_stops: stop times of events
    Output:
        unit_mat: a matrix of participation of each unit to each event
    """
    unit_mat = np.zeros((len(st), len(event_starts)), dtype=int8)
    for rip in range(len(event_starts)):
        for i, s in enumerate(st):
            idx = (s >= event_starts[rip]) & (s <= event_stops[rip])
            unit_mat[i, rip] = sum(idx)
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


def get_participation(st, event_starts, event_stops, par_type="binary"):
    """
    get participation prob.
    make matrix n rows (units) by n cols (ripple epochs)
    Input:
        st: spike train list
        event_starts: event starts
        event_stops: event stops
        par_type: participation type (counts, binary, firing_rate)

        quick binning solution using searchsorted from:
        https://stackoverflow.com/questions/57631469/extending-histogram-function-to-overlapping-bins-and-bins-with-arbitrary-gap
    """
    # convert to numpy array
    event_starts, event_stops = np.array(event_starts), np.array(event_stops)

    # initialize matrix
    unit_mat = np.zeros((len(st), (len(event_starts))))

    # loop over units and bin spikes into epochs
    for i, s in enumerate(st):
        idx1 = np.searchsorted(s, event_starts, 'right')
        idx2 = np.searchsorted(s, event_stops, 'left')
        unit_mat[i, :] = idx2 - idx1

    if par_type == "counts":
        pass
    elif par_type == "binary":
        unit_mat = (unit_mat > 0) * 1
    elif par_type == "firing_rate":
        unit_mat = unit_mat / (event_stops - event_starts)

    return unit_mat


def get_significant_events(scores, shuffled_scores, q=95, tail="both"):
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
    if tail == "both":
        r = np.sum(np.abs(shuffled_scores) >= np.abs(scores), axis=0)
    elif tail == "right":
        r = np.sum(shuffled_scores >= scores, axis=0)
    elif tail == "left":
        r = np.sum(shuffled_scores <= scores, axis=0)
    else:
        raise ValueError("tail must be 'left', 'right', or 'both'")
    pvalues = (r + 1) / (n + 1)

    # set nan scores to 1
    pvalues[np.isnan(scores)] = 1

    sig_event_idx = np.argwhere(
        scores > np.percentile(shuffled_scores, axis=0, q=q)
    ).squeeze()

    return np.atleast_1d(sig_event_idx), np.atleast_1d(pvalues)


def find_laps(
    Vts,
    Vdata,
    newLapThreshold=15,
    good_laps=True,
    edgethresh=0.1,
    completeprop=0.2,
    posbins=50,
):
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

    TL = np.abs(np.nanmax(Vdata) - np.nanmin(Vdata))  # % track length
    th1 = (
        np.nanmin(Vdata) + TL * newLapThreshold / 100
    )  # % lower threshold for lower end
    th2 = (
        np.nanmax(Vdata) - TL * newLapThreshold / 100
    )  # % upper threshold for upper end

    # % loop over all frames
    laps = pd.DataFrame()
    laps.loc[0, "start_ts"] = Vts[0]
    laps.loc[0, "pos"] = Vdata[0]
    laps.loc[0, "start_idx"] = 1
    laps.loc[0, "direction"] = 0
    iLap = 0

    newUpThCross = 1  # % flag for new lap top search
    newDownThCross = 1  # % flag for new lap top search
    for i in range(len(Vdata)):
        if Vdata[i] < th1:  # % search for min
            if newUpThCross == 1:  # % start a new lap
                newUpThCross = 0
                newDownThCross = 1
                iLap = iLap + 1
                laps.loc[iLap, "start_ts"] = Vts[i]
                laps.loc[iLap, "pos"] = Vdata[i]
                laps.loc[iLap, "start_idx"] = i
                laps.loc[iLap, "direction"] = 1

            if Vdata[i] < laps.iloc[iLap].pos:  # % record new min if any
                laps.loc[iLap, "start_ts"] = Vts[i]
                laps.loc[iLap, "pos"] = Vdata[i]
                laps.loc[iLap, "start_idx"] = i

        if Vdata[i] > th2:  # % search for max
            if newDownThCross:  # % start a new lap
                newUpThCross = 1
                newDownThCross = 0
                iLap = iLap + 1
                laps.loc[iLap, "start_ts"] = Vts[i]
                laps.loc[iLap, "pos"] = Vdata[i]
                laps.loc[iLap, "start_idx"] = i
                laps.loc[iLap, "direction"] = -1

            if Vdata[i] > laps.iloc[iLap].pos:  # % record new min if any
                laps.loc[iLap, "start_ts"] = Vts[i]
                laps.loc[iLap, "pos"] = Vdata[i]
                laps.loc[iLap, "start_idx"] = i

    # % fix direction of first lap which was unknown above
    # % make first lap direction opposite of second lap's direction (laps alternate!)
    laps.iloc[0].direction = -laps.iloc[1].direction

    # % make sure laps cross the halfway point
    middle = np.nanmedian(np.arange(np.nanmin(Vdata), np.nanmax(Vdata)))
    i = 0
    while True:
        try:
            positions = np.arange(laps.iloc[i].pos, laps.iloc[i + 1].pos)
        except:
            positions = [np.nan, np.nan]
        if (np.any(positions > middle) == True) & (np.any(positions < middle) == False):
            laps = laps.drop(laps.index[i + 1])
        i = i + 1
        if i + 1 >= len(laps.pos):
            if len(laps.pos) < iLap:
                laps.iloc[0].direction = -laps.iloc[1].direction
            break

    if good_laps:
        laps = find_good_laps(
            Vts,
            Vdata,
            laps,
            edgethresh=edgethresh,
            completeprop=completeprop,
            posbins=posbins,
        )

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
        sys.exit("Input argument delta must be a scalar")

    if delta <= 0:
        sys.exit("Input argument delta must be positive")

    if backwards == 0:
        inc = 1
        first = 0
        last = len(v)
        iter_ = np.arange(first, last, inc)
    elif backwards:
        inc = -1
        first = len(v)
        last = 0
        iter_ = np.arange(first, last, inc)

    mn = np.inf
    mx = -np.inf
    mnpos = np.nan
    mxpos = np.nan

    for ii in iter_:
        this = v[ii]
        if this > mx:
            mx = this
            mxpos = ii
        if this < mn:
            mn = this
            mnpos = ii

        if lookformax:
            try:
                idx = mx - delta > mintab[-1]
            except:
                idx = mx - delta > mintab

            if (this < mx - delta) | ((ii == last - 1) & (len(mintab) > 0) & idx):
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = ii
                lookformax = 0
        else:
            try:
                idx = mx - delta < maxtab[-1]
            except:
                idx = mx - delta < maxtab
            if (this > mn + delta) | ((ii == last - 1) & (len(maxtab) > 0) & idx):
                mintab.append((mnpos, mn))
                mx = this
                mxpos = ii
                lookformax = 1

    if (len(maxtab) == 0) & (len(mintab) == 0):
        if lookformax:
            if mx - mn > delta:
                maxtab = [mxpos, mx]
        else:
            if mx - mn > delta:
                mintab = [mnpos, mn]
    return maxtab, mintab


def find_good_laps(ts, V_rest, laps, edgethresh=0.1, completeprop=0.2, posbins=50):
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

    if (
        edgethresh > 1
    ):  # % in case edgethresh is input as a percentage instead of a proportion
        edgethresh = edgethresh / 100

    if (
        completeprop > 1
    ):  # % in case completeprop is input as a percentage instead of a proportion
        completeprop = completeprop / 100

    bottomend = np.nanmin(V_rest)
    topend = np.nanmax(V_rest)
    bins = np.arange(
        bottomend,
        topend + (topend - bottomend) / posbins,
        (topend - bottomend) / posbins,
    )
    # % threshold for peak/trough detection
    delta = (topend - bottomend) * edgethresh
    startgoodlaps = []
    stopgoodlaps = []

    l = 0
    while l < len(laps) - 1:
        # % select out just this lap
        if l == len(laps):
            endoflap = ts[-1]
        else:
            endoflap = laps.iloc[l + 1].start_ts

        v = V_rest[
            np.where(ts == laps.iloc[l].start_ts)[0][0] : np.where(ts == endoflap)[0][0]
        ]
        t = ts[
            np.where(ts == laps.iloc[l].start_ts)[0][0] : np.where(ts == endoflap)[0][0]
        ]

        # % find turn around points during this lap
        lookformax = laps.iloc[l].direction == 1
        peak, trough = peakdetz(v, delta, lookformax, 0)

        if lookformax:
            # % find the direct path from bottomend to topend (or mark lap for
            # % deleting if the turn around points are not in those ranges)
            if len(trough) > 0:
                # % find the last trough in range of bottomend (start of lap)
                gt = len(trough)
                while (gt > 0) & (trough(gt, 2) >= 2 * delta + bottomend):
                    gt = gt - 1

                # % assign the next peak after that trough as the end of the lap
                # % (or mark lap for deleting, if that peak is not at topend)
                if gt == 0:
                    if peak[1, 2] > topend - 2 * delta:
                        t = t[0 : peak[0]]
                        v = v[0 : peak[0]]
                    else:
                        # % this marks the lap for deleting
                        t = t[0:5]
                        v = v[0:5]
                else:
                    et = len(peak)
                    if gt + 1 > et:
                        gt = 0
                        t = t[0:2]
                        v = v[0:2]
                    else:
                        t = t[trough[gt, 1] : peak[gt + 1, 1]]
                        v = v[trough[gt, 1] : peak[gt + 1, 1]]

            else:
                # % make sure peak exists and is in range of topend
                if len(peak) == 0:
                    if len(t) > 2:
                        t = t[0:2]
                        v = v[0:2]
                elif peak[1] < topend - 2 * delta:
                    # % this marks the lap for deleting
                    if len(t) > 5:
                        t = t[0:5]
                        v = v[0:5]
        else:  # % if lookformax
            # % find the direct path from topend to bottomend (or mark lap for
            # % deleting if the turn around points are not in those ranges)
            if len(peak) > 0:
                # % find the last peak in range of topend (start of lap)
                gt = len(peak)
                while (gt > 0) & (peak[gt, 2] <= topend - 2 * delta):
                    gt = gt - 1
                # % assign the next trough after that peak as the end of the lap
                # % (or mark lap for deleting, if that trough is not at bottomend)
                if gt == 0:
                    if trough(1, 2) < bottomend + 2 * delta:
                        t = t[1 : trough[0]]
                        v = v[1 : trough[0]]
                    else:
                        # % this marks the lap for deleting
                        t = t[0:5]
                        v = v[0:5]
                else:
                    et = len(trough)
                    if gt + 1 > et:
                        t = t[0:2]
                        v = v[0:2]
                        gt = 0
                    else:
                        t = t[peak[gt, 1] : trough[gt + 1, 1]]
                        v = v[peak[gt, 1] : trough[gt + 1, 1]]
            else:  # % if ~isempty(peak)
                # % make sure trough exists and is in range of bottomend
                if len(trough) == 0:
                    if len(t) > 2:
                        t = t[0:2]
                        v = v[0:2]

                elif trough[1] > bottomend + 2 * delta:
                    # % this marks the lap for deleting
                    if len(t) > 5:
                        t = t[0:5]
                        v = v[0:5]
        vcovered, _ = np.histogram(v, bins=bins)

        if len(v) < 3:
            # % eliminate the lap if it is non-existent (as is sometimes the case for lap 1)
            laps = laps.drop(laps.index[l])
        # % eliminate the lap if >completeprop of it is NaNs or if it has been marked for
        # % deleting above
        elif (len(v) < 6) | (sum(vcovered == 0) > completeprop * posbins):
            laps.drop(laps.index[l])
            # % remove the other lap from the lap pair
            if l % 2 == 0:
                # % delete previous lap from laps
                laps = laps.drop(laps.index[l - 1])
                # % change goodlaps markers to delete previous lap from laps
                if len(stopgoodlaps) > 0:
                    if "lastlapend" not in locals() | (startgoodlaps[-1] > lastlapend):
                        startgoodlaps[-1] = []
                        stopgoodlaps[-1] = []
                    else:
                        stopgoodlaps[-1] = lastlapend

                l = l - 1
            elif l <= len(laps) & l > 1:
                # % delete next lap from laps
                laps = laps.drop(laps.index[l])
        else:  # % if lap is good
            # % store last lap end just in case have to delete this lap with next lap
            if len(stopgoodlaps) > 0:
                lastlapend = stopgoodlaps[-1]

            # % add this lap to goodlaps
            try:
                idx = stopgoodlaps[-1] == t[0]
            except:
                idx = stopgoodlaps == t[0]
            if (len(stopgoodlaps) > 0) & (idx):
                stopgoodlaps[-1] = t[-1]
            else:
                startgoodlaps.append(t[0])
                stopgoodlaps.append(t[-1])

            l = l + 1

    return laps


def get_linear_track_lap_epochs(
    ts,
    x,
    newLapThreshold=15,
    good_laps=False,
    edgethresh=0.1,
    completeprop=0.2,
    posbins=50,
):
    """
    get_linear_track_lap_epochs: def that calls find_laps and outputs nelpy epochs
        for out and inbound running directions
    """
    laps = find_laps(
        np.array(ts),
        np.array(x),
        newLapThreshold=newLapThreshold,
        good_laps=good_laps,
        edgethresh=edgethresh,
        completeprop=completeprop,
        posbins=posbins,
    )

    outbound_start = []
    outbound_stop = []
    inbound_start = []
    inbound_stop = []

    for i in range(len(laps) - 1):
        if laps.iloc[i].direction == 1:
            outbound_start.append(laps.iloc[i].start_ts)
            outbound_stop.append(laps.iloc[i + 1].start_ts)

        if laps.iloc[i].direction == -1:
            inbound_start.append(laps.iloc[i].start_ts)
            inbound_stop.append(laps.iloc[i + 1].start_ts)

    outbound_epochs = nel.EpochArray([np.array([outbound_start, outbound_stop]).T])
    inbound_epochs = nel.EpochArray([np.array([inbound_start, inbound_stop]).T])

    return outbound_epochs, inbound_epochs


def find_pre_task_post(env, pre_post_label="sleep"):
    """
    given list of environment, finds first contigous epochs that meet pre/task/post

    Input:
        environment list, can be pandas column
    Output:
        indices of where pre-sleep/task/post-sleep exist

    example:
    pre_task_post = find_pre_task_post(epoch_df.environment)
    epoch_df.loc[pre_task_post]

            name	                        startTime	stopTime	environment
        1	OR15day1_sleep1_180116_110120	2001.600	8087.29195	sleep
        2	OR15day1_2_180116_171020	    8087.292	9952.05995	wmaze
        3	OR15day1_sleep2_180116_181618	9952.060	10182.92795	sleep
    """
    if len(env) < 3:
        return None, None
    numeric_idx = (pre_post_label == env) * 1
    dummy = np.zeros_like(numeric_idx) == 1
    if all(numeric_idx[:3] == [1, 0, 1]):
        dummy[:3] = True
        return dummy, [0, 1, 2]
    else:
        for i in np.arange(len(numeric_idx) + 3):
            if 3 + i > len(numeric_idx):
                return None, None
            if all(numeric_idx[0 + i : 3 + i] == [1, 0, 1]):
                dummy[0 + i : 3 + i] = True
                return dummy, [0, 1, 2] + i


def find_epoch_pattern(env, pattern):
    """
    given list of environment, finds contigous epochs that meet pattern

    Limitation: stops on the first instance of finding the pattern

    Input:
        env: environment list, can be pandas column
        pattern: pattern you are searching for
    Output:
        indices of where pattern exist

    example:
    epoch_df = loading.load_epoch(basepath)
    pattern_idx,_ = find_epoch_pattern(epoch_df.environment,['sleep','linear','sleep'])
    epoch_df.loc[pattern_idx]

        name	                startTime	stopTime	environment	behavioralParadigm	notes
    0	preSleep_210411_064951	0.0000	    9544.56315	sleep	    NaN	                NaN
    1	maze_210411_095201	    9544.5632	11752.80635	linear	    novel	            novel
    2	postSleep_210411_103522	11752.8064	23817.68955	sleep	    novel	            novel
    """

    env = list(env)
    pattern = list(pattern)

    if len(env) < len(pattern):
        return None, None

    dummy = np.zeros(len(env))

    for i in range(len(env) - len(pattern) + 1):
        if pattern == env[i : i + len(pattern)]:
            dummy[i : i + len(pattern)] = 1
            dummy = dummy == 1
            return dummy, np.arange(i, i + len(pattern))
    return None, None


def get_rank_order(
    st,
    epochs,
    method="peak_fr",  # 'first_spike' or 'peak_fr'
    ref="cells",  # 'cells' or 'epochs'
    padding=0.05,
    dt=0.001,
    sigma=0.01,
    min_units=5,
):
    """
    get rank order of spike train within epoch
    Input:
        st: spike train nelpy array
        epochs: epoch array, windows in which to calculate rank order
        method: method of rank order 'first_spike' or 'peak_fr' (default: peak_fr)
        ref: frame of reference for rank order ('cells' or 'epoch') (default: cells)
        padding: +- padding for epochs
        dt: bin width (s) for finding relative time (epoch ref)
        sigma: smoothing sigma (s) (peak_fr method)
    Output:
        median_rank: median rank order over all epochs (0-1)
        rank_order: matrix (n cells X n epochs) each column shows rank per cell per epoch (0-1)

    Example:
        st,_ = loading.load_spikes(basepath,putativeCellType='Pyr')
        forward_replay = nel.EpochArray(np.array([starts,stops]).T)
        median_rank,rank_order = get_rank_order(st,forward_replay)
    """
    # filter out specific warnings
    warnings.filterwarnings(
        "ignore", message="ignoring events outside of eventarray support"
    )
    warnings.filterwarnings("ignore", message="Mean of empty slice")

    if method not in ["first_spike", "peak_fr"]:
        assert Exception("method " + method + " not implemented")
    if ref not in ["cells", "epoch"]:
        assert Exception("ref " + ref + " not implemented")

    def get_min_ts(st_temp):
        min_ts = []
        for ts in st_temp.data:
            # nan if no spikes
            if len(ts) == 0:
                min_ts.append(np.nan)
            else:
                min_ts.append(np.nanmin(ts))
        return min_ts

    def rank_order_first_spike(st_epoch, epochs, dt, min_units, ref):
        # set up empty matrix for rank order
        rank_order = np.ones([st_epoch.data.shape[0], st_epoch.n_intervals]) * np.nan

        unit_id = np.arange(st_epoch.data.shape[0])

        # iter over every event
        for event_i, st_temp in enumerate(st_epoch):

            if ref == "cells":
                # get firing order
                idx = np.array(st_temp.get_event_firing_order()) - 1
                # reorder unit ids by order and remove non-active
                units = unit_id[idx][st_temp.n_events[idx] > 0]
                # how many are left?
                nUnits = len(units)

                if nUnits < min_units:
                    rank_order[:, event_i] = np.nan
                else:
                    # arange 1 to n units in order of units
                    rank_order[units, event_i] = np.arange(nUnits)
                    # normalize by n units
                    rank_order[units, event_i] = rank_order[units, event_i] / nUnits
            elif ref == "epoch":
                # find first spike time for each cell
                min_ts = get_min_ts(st_temp)
                # make time stamps for interpolation
                epoch_ts = np.arange(epochs[event_i].start, epochs[event_i].stop, dt)
                # make normalized range 0-1
                norm_range = np.linspace(0, 1, len(epoch_ts))
                # get spike order relative to normalized range
                if len(min_ts) < min_units:
                    rank_order[:, event_i] = np.nan
                else:
                    rank_order[:, event_i] = np.interp(min_ts, epoch_ts, norm_range)
        return rank_order

    def rank_order_fr(st_epoch, dt, sigma, min_units, ref):
        # set up empty matrix for rank order
        rank_order = np.zeros([st_epoch.data.shape[0], st_epoch.n_intervals]) * np.nan

        unit_id = np.arange(st_epoch.data.shape[0])

        # bin spike train here (smooth later per epoch to not have edge issues)
        z_t = st_epoch.bin(ds=dt)
        # iter over epochs
        for event_i, z_t_temp in enumerate(z_t):
            # smooth spike train in order to estimate peak
            z_t_temp.smooth(sigma=sigma, inplace=True)

            if ref == "cells":

                # find loc of each peak and get sorted idx of active units
                idx = np.argsort(np.argmax(z_t_temp.data, axis=1))
                # reorder unit ids by order and remove non-active
                units = unit_id[idx][np.sum(z_t_temp.data[idx, :] > 0, axis=1) > 0]

                nUnits = len(units)

                if nUnits < min_units:
                    rank_order[:, event_i] = np.nan
                else:
                    # arange 1 to n units in order of units
                    rank_order[units, event_i] = np.arange(nUnits)
                    # normalize by n units
                    rank_order[units, event_i] = rank_order[units, event_i] / nUnits
            elif ref == "epoch":
                # iterate over each cell
                for cell_i, unit in enumerate(z_t_temp.data):
                    # if the cell is not active apply nan
                    if not np.any(unit > 0):
                        rank_order[cell_i, event_i] = np.nan
                    else:
                        # calculate normalized rank order (0-1)
                        rank_order[cell_i, event_i] = np.argmax(unit) / len(unit)
        return rank_order

    # create epoched spike array
    st_epoch = st[epochs.expand(padding)]

    # if no spikes in epoch, break out
    if st_epoch.n_active == 0:
        return np.tile(np.nan, st.data.shape), None

    # set up empty matrix for rank order
    if method == "peak_fr":
        rank_order = rank_order_fr(st_epoch, dt, sigma, min_units, ref)
    elif method == "first_spike":
        rank_order = rank_order_first_spike(st_epoch, epochs, dt, min_units, ref)
    else:
        raise Exception("other method, " + method + " is not implemented")

    return np.nanmedian(rank_order, axis=1), rank_order

def overlap_intersect(epoch,interval,return_indices=True):
    """
    Returns the epochs with overlap with interval
    Input:
        epoch: nelpy.EpochArray
            The epochs to check
        interval: nelpy.IntervalArray
            The interval to check for overlap
        return_indices: bool
            If True, returns the indices of the epochs (interval) that overlap    
    Output:
        epoch: nelpy.EpochArray
            The epochs with overlap with interval
    """
    new_intervals = []
    indices = []
    for epa in (epoch):
        if any((interval.starts < epa.stop) & (interval.stops > epa.start)):
            new_intervals.append([epa.start, epa.stop])
            cand_ep_idx = np.where((interval.starts < epa.stop) & (interval.stops > epa.start))
            indices.append(cand_ep_idx[0][0])
    out = type(epoch)(new_intervals)
    out._domain = epoch.domain
    if return_indices:
        return out, indices
    return out

def get_velocity(position, time=None):
    if time is None:
        time = np.arange(position.shape[0])
    return np.gradient(position, time, axis=0)

def get_speed(position, time=None):
    velocity = get_velocity(position, time=time)
    return np.sqrt(np.sum(velocity ** 2, axis=1))