import itertools
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from numba import jit,njit
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
from ripple_heterogeneity.utils import compress_repeated_epochs as comp_rep_ep
import random
from nelpy import core
from ripple_heterogeneity.utils import loading
from scipy.linalg import toeplitz


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


def pairwise_corr(X, method="pearson", pairs=None):
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
        if method == "pearson":
            rho_, pval_ = stats.pearsonr(X[s[0], :], X[s[1], :])
        elif method == "spearman":
            rho_, pval_ = stats.spearmanr(X[s[0], :], X[s[1], :])
        elif method == "kendall":
            rho_, pval_ = stats.kendalltau(X[s[0], :], X[s[1], :])
        else:
            raise ValueError("method must be pearson, spearman or kendall")
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
    times = np.linspace(-(nbins * binsize) / 2, (nbins * binsize) / 2, nbins + 1)

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

    spatial_corr = []
    # Now we can iterate over spikes
    for i, s in enumerate(pairs):
        # Calling the crossCorr function
        x1 = X[s[0], :, :].flatten()
        x2 = X[s[1], :, :].flatten()
        bad_idx = np.isnan(x1) | np.isnan(x1)
        spatial_corr.append(np.corrcoef(x1[~bad_idx], x2[~bad_idx])[0, 1])

    if return_index:
        return np.array(spatial_corr), pairs
    else:
        return np.array(spatial_corr)


def compute_psth(spikes, event, bin_width=0.002, n_bins=100):

    # times = np.arange(0, bin_width * (n_bins + 1), bin_width) - (n_bins * bin_width) / 2
    times = np.linspace(-(n_bins * bin_width) / 2, (n_bins * bin_width) / 2, n_bins + 1)
    ccg = pd.DataFrame(index=times, columns=np.arange(len(spikes)))
    # Now we can iterate over spikes
    for i, s in enumerate(spikes):
        ccg[i] = crossCorr(event, s, bin_width, n_bins)
    return ccg

def deconvolve_peth(signal,events,bin_width=0.002, n_bins=100):
    """
    This function performs deconvolution of a peri-event time histogram (PETH) signal.
    
    Parameters:
    signal (array): An array representing the discrete events.
    events (array): An array representing the discrete events.
    bin_width (float, optional): The width of a time bin in seconds (default value is 0.002 seconds).
    n_bins (int, optional): The number of bins to use in the PETH (default value is 100 bins).
    
    Returns:
    deconvolved (array): An array representing the deconvolved signal.
    times (array): An array representing the time points corresponding to the bins.
    
    Based on DeconvolvePETH.m from https://github.com/ayalab1/neurocode/blob/master/spikes/DeconvolvePETH.m
    """

    # calculate time lags for peth
    times = np.linspace(-(n_bins * bin_width) / 2, (n_bins * bin_width) / 2, n_bins + 1)

    # Calculate the autocorrelogram of the signal and the PETH of the events and the signal
    autocorrelogram = crossCorr(signal, signal, bin_width, n_bins*2)
    raw_peth = crossCorr(events, signal, bin_width, n_bins*2)

    # Subtract the mean value from the raw_peth
    const = np.mean(raw_peth)
    raw_peth = raw_peth - const

    # Calculate the Toeplitz matrix using the autocorrelogram and
    #   the cross-correlation of the autocorrelogram
    T0 = toeplitz(
        autocorrelogram,
        np.hstack([autocorrelogram[0], np.zeros(len(autocorrelogram)-1)])
    )
    T = T0[n_bins:, :n_bins+1]

    # Calculate the deconvolved signal by solving a linear equation
    deconvolved = np.linalg.solve(
        T, raw_peth[int(n_bins / 2) : int(n_bins / 2 * 3 + 1)].T + const / len(events)
    )   
    return deconvolved, times


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


def get_raster_points(data, time_ref, bin_width=0.002, n_bins=100, window=None):
    """
    Generate points for a raster plot centered around each reference time in the `time_ref` array.

    Parameters
    ----------
    data : ndarray
        A 1D array of time values.
    time_ref : ndarray
        A 1D array of reference times.
    bin_width : float, optional
        The width of each bin in the raster plot, in seconds. Default is 0.002 seconds.
    n_bins : int, optional
        The number of bins in the raster plot. Default is 100.
    window : tuple, optional
        A tuple containing the start and end times of the window to be plotted around each reference time.
        If not provided, the window will be centered around each reference time and have a width of `n_bins * bin_width` seconds.

    Returns
    -------
    x : ndarray
        A 1D array of x values representing the time offsets of each data point relative to the corresponding reference time.
    y : ndarray
        A 1D array of y values representing the reference times.
    times : ndarray
        A 1D array of time values corresponding to the bins in the raster plot.
    """

    if window is not None:
        times = np.arange(window[0], window[1] + bin_width/2, bin_width)
    else:
        times = np.linspace(-(n_bins * bin_width) / 2, (n_bins * bin_width) / 2, n_bins + 1)
    x = []
    y = []
    for i, r in enumerate(time_ref):
        idx = (data > r + times.min()) & (data < r + times.max())
        cur_data = data[idx]
        # if any(cur_data):
        x.append(cur_data - r)
        y.append(np.ones_like(cur_data)+i)       
    x = list(itertools.chain(*x))
    y = list(itertools.chain(*y))
    return x, y, times

def peth_matrix(data, time_ref, bin_width=0.002, n_bins=100, window=None):

    x, y, t = get_raster_points(
        data, time_ref, bin_width=bin_width, n_bins=n_bins, window=window
    )
    dt = np.diff(t)[0]
    x, y = np.array(x), np.array(y)
    H, xedges, yedges = np.histogram2d(
        x, y, bins=(
            np.arange(t.min(),t.max()+dt,dt),
            np.arange(.5,len(time_ref) + 1.5)
            )
    )
    return H, t[:-1]+dt/2

def event_triggered_average_irregular_sample(
    timestamps, data, time_ref, bin_width=0.002, n_bins=100, window=None
):
    """
    Compute the average and standard deviation of data values within a window around 
    each reference time.

    Specifically for irregularly sampled data

    Parameters
    ----------
    timestamps : ndarray
        A 1D array of times associated with data.
    data : ndarray
        A 1D array of data values.
    time_ref : ndarray
        A 1D array of reference times.
    bin_width : float, optional
        The width of each bin in the window, in seconds. Default is 0.002 seconds.
    n_bins : int, optional
        The number of bins in the window. Default is 100.
    window : tuple, optional
        A tuple containing the start and end times of the window to be plotted around each reference time.
        If not provided, the window will be centered around each reference time and have a 
        width of `n_bins * bin_width` seconds.

    Returns
    -------
    pd.DataFrame, pd.DataFrame
        two dataframes, the first containing the average values, the second the 
        standard deviation of data values within the window around each reference time.
    """

    if window is not None:
        times = np.arange(window[0], window[1] + bin_width, bin_width)
    else:
        times = np.linspace(
            -(n_bins * bin_width) / 2, (n_bins * bin_width) / 2, n_bins + 1
        )
    x = []
    y = []
    for i, r in enumerate(time_ref):
        idx = (timestamps > r + times.min()) & (timestamps < r + times.max())
        x.append((timestamps - r)[idx])
        y.append(data[idx])

    temp_df = pd.DataFrame()
    if len(x) == 0:
        return temp_df, temp_df
    temp_df["time"] = np.hstack(x)
    temp_df["data"] = np.hstack(y)
    temp_df = temp_df.sort_values(by="time",ascending=True)

    average_val = np.zeros(len(times)-1)
    std_val = np.zeros(len(times)-1)
    for i in range(len(times)-1):
        average_val[i] = temp_df[temp_df.time.between(times[i], times[i+1])].data.mean()
        std_val[i] = temp_df[temp_df.time.between(times[i], times[i+1])].data.std()

    avg = pd.DataFrame(index=times[:-1] + bin_width/2)
    avg[0] = average_val

    std = pd.DataFrame(index=times[:-1] + bin_width/2)
    std[0] = std_val

    return avg, std


def event_triggered_average(
        timestamps:np.ndarray, signal:np.ndarray, events:np.ndarray, window=[-0.5, 0.5]
) -> np.ndarray:
    """
    Calculates the spike-triggered averages of signals in a time window
    relative to the event times of a corresponding events for multiple
    signals each. The function receives n signals and either one or
    n events. In case it is one event this one is muliplied n-fold
    and used for each of the n signals.

    adapted from elephant.sta.spike_triggered_average to be used with ndarray

    Parameters
    ----------
    timestamps : ndarray (n samples)

    signal : ndarray (n samples x n signals)

    events : one numpy ndarray or a list of n of either of these.

    window : tuple of 2.
        'window' is the start time and the stop time, relative to a event, of
        the time interval for signal averaging.
        If the window size is not a multiple of the sampling interval of the
        signal the window will be extended to the next multiple.

    Returns
    -------
    result_sta : ndarray
        'result_sta' contains the event-triggered averages of each of the
        signals with respect to the event in the corresponding
        events. The length of 'result_sta' is calculated as the number
        of bins from the given start and stop time of the averaging interval
        and the sampling rate of the signal. If for an signal
        no event was either given or all given events had to be ignored
        because of a too large averaging interval, the corresponding returned
        signal has all entries as nan.


    Examples
    --------

    >>> m1 = assembly_reactivation.AssemblyReact(basepath=r"Z:\Data\HMC2\day5")

    >>> m1.load_data()
    >>> m1.get_weights(epoch=m1.epochs[1])
    >>> assembly_act = m1.get_assembly_act()

    >>> peth_avg, time_lags = event_triggered_average(
    ...    assembly_act.data.T, assembly_act.abscissa_vals, m1.ripples.starts, window=[-0.5, 0.5]
    ... )

    >>> plt.plot(time_lags,peth_avg)
    >>> plt.show()
    """
    window_starttime, window_stoptime = window

    if len(signal.shape) == 1:
        signal = signal[:, np.newaxis]

    num_signals = signal.shape[1]

    sampling_rate = 1 / np.diff(timestamps)[0]
    # window_bins: number of bins of the chosen averaging interval
    window_bins = int(np.ceil(((window_stoptime - window_starttime) * sampling_rate)))
    # result_sta: array containing finally the spike-triggered averaged signal
    result_sta = np.zeros((window_bins, num_signals))
    # setting of correct times of the spike-triggered average
    # relative to the spike
    time_lags = np.arange(window_starttime, window_stoptime, 1 / sampling_rate)

    used_events = np.zeros(num_signals, dtype=int)
    total_used_events = 0

    for i in range(num_signals):
        # summing over all respective signal intervals around spiketimes
        for event in events:
            # checks for sufficient signal data around spiketime
            if (
                event + window_starttime >= timestamps[0]
                and event + window_stoptime <= timestamps[-1]
            ):
                # calculating the startbin in the analog signal of the
                # averaging window for event
                startbin = int(
                    np.floor(
                        ((event + window_starttime - timestamps[0]) * sampling_rate)
                    )
                )
                # adds the signal in selected interval relative to the spike
                result_sta[:, i] += signal[startbin : startbin + window_bins, i]
                # counting of the used event
                used_events[i] += 1

        # normalization
        result_sta[:, i] = result_sta[:, i] / used_events[i]

        total_used_events += used_events[i]

    if total_used_events == 0:
        warnings.warn("No events at all was either found or used for averaging")

    return result_sta, time_lags


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

        results = pd.concat([results, temp_df], ignore_index=True)
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

    print(
        "functions.find_sig_assemblies: is deprecated, use find_sig_assembly.main instead"
    )
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
        idx1 = np.searchsorted(s, event_starts, "right")
        idx2 = np.searchsorted(s, event_stops, "left")
        unit_mat[i, :] = idx2 - idx1

    par_type_funcs = {
        "counts": lambda x: x,
        "binary": lambda x: (x > 0) * 1,
        "firing_rate": lambda x: x / (event_stops - event_starts),
    }
    calc_func = par_type_funcs[par_type]
    unit_mat = calc_func(unit_mat)

    return unit_mat


def get_significant_events(scores, shuffled_scores, q=95, tail="both"):
    """
    Return the significant events based on percentiles,
    the p-values and the standard deviation of the scores 
    in terms of the shuffled scores.
    Parameters
    ----------
    scores : array of shape (n_events,)
        The array of scores for which to calculate significant events
    shuffled_scores : array of shape (n_shuffles, n_events)
        The array of scores obtained from randomized data 
    q : float in range of [0,100]
        Percentile to compute, which must be between 0 and 100 inclusive.
    Returns
    -------
    sig_event_idx : array of shape (n_sig_events,)
        Indices (from 0 to n_events-1) of significant events.
    pvalues : array of shape (n_events,)
        The p-values 
    stddev : array of shape (n_events,)
        The standard deviation of the scores in terms of the shuffled scores
    """
    # check shape and correct if needed
    if isinstance(scores,list) | isinstance(scores,np.ndarray):
        if shuffled_scores.shape[1] != len(scores):
            shuffled_scores = shuffled_scores.T

    n = shuffled_scores.shape[0]
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
    if isinstance(np.isnan(scores),np.ndarray):
        pvalues[np.isnan(scores)] = 1

    sig_event_idx = np.argwhere(
        scores > np.percentile(shuffled_scores, axis=0, q=q)
    ).squeeze()

    # calculate how many standard deviations away from shuffle
    stddev = (np.abs(scores) - np.nanmean(np.abs(shuffled_scores), axis=0)) / np.nanstd(
        np.abs(shuffled_scores), axis=0
    )

    return np.atleast_1d(sig_event_idx), np.atleast_1d(pvalues), np.atleast_1d(stddev)


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


def find_good_lap_epochs(pos, dir_epoch, thres=0.5, binsize=6, min_laps=10):
    """
    find_good_laps: finds good laps in behavior data
        Made to find good laps in nelpy array for replay analysis
    input:
        pos: nelpy analog array with single dim
        dir_epoch: epoch to find good lap
        thres: occupancy threshold for good lap
        binsize: size of bins to calculate occupancy
    output:
        good_laps: epoch array of good laps
    """
    # make bin edges to calc occupancy
    x_edges = np.arange(np.nanmin(pos.data[0]), np.nanmax(pos.data[0]), binsize)
    # initialize occupancy matrix (position x time)
    occ = np.zeros([len(x_edges) - 1, dir_epoch.n_intervals])
    # iterate through laps
    for i, ep in enumerate(dir_epoch):
        # bin position per lap
        occ[:, i], _ = np.histogram(pos[ep].data[0], bins=x_edges)
    # calc percent occupancy over position bins per lap and find good laps
    good_laps = np.where(~((np.sum(occ == 0, axis=0) / occ.shape[0]) > thres))[0]
    # if no good laps, return empty epoch
    if (len(good_laps) == 0) | (len(good_laps) < min_laps):
        dir_epoch = nel.EpochArray()
    else:
        dir_epoch = dir_epoch[good_laps]
    return dir_epoch


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

def find_multitask_pre_post(env, task_tag = 'open_field|linear_track|box|tmaze|wmaze'):
    """
    Find the row index for pre_task/post_task sleep for a given enviornment from cell explorer session.epochs dataframe 
    Returns list of pre/task_post epochs for each task. 
    input: 
        df: data frame consisting of cell explorer session.epochs data 
        col_name: name of column to query task tag. Default is environemnt
        task_tag: string within col_name that indicates a task. 
    output: 
        list of epoch indicies [pre_task, task, post_task] of size n = # of task epochs

    LB/RH 1/5/2022    
    """
    # Find the row indices that contain the search string in the specified column
    task_bool = env.str.contains(task_tag, case=False)
    sleep_bool = env.str.contains('sleep', case=False)

    task_idx = np.where(task_bool)[0]
    task_idx = np.delete(task_idx,task_idx == 0,0)
    sleep_idx = np.where(sleep_bool)[0]
    
    pre_task_post = []
    for task in task_idx:
        temp = sleep_idx - task
        pre_task = sleep_idx[temp < 0]
        post_task = sleep_idx[temp > 0]

        if len(post_task) == 0:
            print('no post_task sleep for task epoch '+str(task))
        elif len(pre_task) == 0:
            print('no pre_task sleep for task epoch '+str(task))
        else:
            pre_task_post.append([pre_task[-1], task, post_task[0]])

    if len(pre_task_post) == 0:
        pre_task_post = None
        
    return pre_task_post

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


def find_env_paradigm_pre_task_post(epoch_df, env="sleep", paradigm="memory"):
    """
    find_env_paradigm_pre_task_post: use env and paradigm to find pre task post
    Made because: FujisawaS data has Spontaneous alternation task & Working memory task
        both flanked by sleep. We want to locate the working memory task pre/task/post
    ex.

    >> epoch_df
        name	startTime	stopTime	environment	behavioralParadigm	            notes
    0	EE.042	0.0	        995.9384	sleep	    NaN	                            NaN
    1	EE.045	995.9384	3336.3928	tmaze	    Spontaneous alternation task	NaN
    2	EE.046	3336.3928	5722.444	sleep	    NaN	                            NaN
    3	EE.049	5722.444	7511.244	tmaze	    Working memory task	            NaN
    4	EE.050	7511.244	9387.644	sleep	    NaN	                            NaN

    >> idx = find_env_paradigm_pre_task_post(epoch_df)
    >> idx
    array([False, False,  True,  True,  True])

    """
    # compress back to back sleep epochs
    epoch_df_ = comp_rep_ep.main(epoch_df, epoch_name="sleep")
    # make col with env and paradigm
    epoch_df_["sleep_ind"] = (
        epoch_df_.environment + "_" + epoch_df_.behavioralParadigm.astype(str)
    )
    # locate env and paradigm of choice with this col
    epoch_df_["sleep_ind"] = epoch_df_["sleep_ind"].str.contains(env + "|" + paradigm)
    # the pattern we are looking for is all True

    # https://stackoverflow.com/questions/48710783/pandas-find-and-index-rows-that-match-row-sequence-pattern
    pat = np.asarray([True, True, True])
    N = len(pat)
    idx = (
        epoch_df_["sleep_ind"]
        .rolling(window=N, min_periods=N)
        .apply(lambda x: (x == pat).all())
        .mask(lambda x: x == 0)
        .bfill(limit=N - 1)
        .fillna(0)
        .astype(bool)
    ).values
    return idx

def add_animal_id(df:pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    df["animal_id"] = df.basepath.map(
        dict(
            [
                (basepath, loading.get_animal_id(basepath))
                    for basepath in df.basepath.unique()
            ]
        )
    )
    return df

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


def randomize_epochs(epoch, randomize_each=True, start_stop=None):
    """Randomly shifts the epochs of a EpochArray object and wraps them around the original time boundaries.

    This method takes a EpochArray object as input, and can either randomly shift each epoch by a different amount
    (if `randomize_each` is True) or shift all the epochs by the same amount (if `randomize_each` is False).
    In either case, the method wraps the shifted epochs around the original time boundaries to make sure they remain
    within the original time range. It then returns the modified EpochArray object.

    Args:
        epoch (EpochArray): The EpochArray object whose epochs should be shifted and wrapped.
        randomize_each (bool, optional): If True, each epoch will be shifted by a different random amount.
            If False, all the epochs will be shifted by the same random amount. Defaults to True.
        start_stop (array, optional): If not None, time support will be taken from start_stop

    Returns:
        new_epochs: The modified EpochArray object with the shifted and wrapped epochs.
    """

    def wrap_intervals(intervals, start, stop):
        idx = np.any(intervals > stop, axis=1)
        intervals[idx] = intervals[idx] - stop + start

        idx = np.any(intervals < start, axis=1)
        intervals[idx] = intervals[idx] - start + stop
        return intervals

    new_epochs = epoch.copy()

    if start_stop is None:
        start = new_epochs.start
        stop = new_epochs.stop
    else:
        start, stop = start_stop

    ts_range = stop - start

    if randomize_each:
        # Randomly shift each epoch by a different amount
        random_order = random.sample(
            range(-int(ts_range), int(ts_range)), new_epochs.n_intervals
        )

        new_intervals = new_epochs.data + np.expand_dims(random_order, axis=1)
        new_epochs._data = wrap_intervals(new_intervals, start, stop)
    else:
        # Shift all the epochs by the same amount
        random_shift = random.randint(-int(ts_range), int(ts_range))
        new_epochs._data = wrap_intervals((new_epochs.data + random_shift), start, stop)

    if not new_epochs.isempty:
        if np.any(new_epochs.data[:, 1] - new_epochs.data[:, 0] < 0):
            raise ValueError("start must be less than or equal to stop")
            
    new_epochs._sort()

    return new_epochs


def overlap_intersect(epoch, interval, return_indices=True):
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
    for epa in epoch:
        if any((interval.starts < epa.stop) & (interval.stops > epa.start)):
            new_intervals.append([epa.start, epa.stop])
            cand_ep_idx = np.where(
                (interval.starts < epa.stop) & (interval.stops > epa.start)
            )
            indices.append(cand_ep_idx[0][0])
    out = type(epoch)(new_intervals)
    out._domain = epoch.domain
    if return_indices:
        return out, indices
    return out

@jit(nopython=True)
def find_intersecting_intervals_(set1, set2):

    intersecting_intervals = []
    for i, (start1, end1) in enumerate(set1):
        # Check if any of the intervals in set2 intersect with the current interval in set1
        intersects = False
        for start2, end2 in set2:
            if start2 <= end1 and end2 >= start1:
                intersects = True
                break
        intersecting_intervals.append(intersects)
    return intersecting_intervals

def find_intersecting_intervals(set1, set2):
    """
    Find whether each interval in set1 intersects with any interval in set2.

    Parameters
    ----------
    set1 : nelpy EpochArray
    set2 : nelpy EpochArray

    Returns
    -------
    list
        A list of bools, where each bool indicates whether the corresponding interval in set1 intersects with any interval in set2.

    Examples
    --------
    >>> set1 = nel.EpochArray([(1, 3), (5, 7), (9, 10)])
    >>> set2 = nel.EpochArray([(2, 4), (6, 8)])
    >>> find_intersecting_intervals(set1, set2)
    [True, True, False]
    """
    if not isinstance(set1, core.IntervalArray) & isinstance(set2, core.IntervalArray):
        raise ValueError("only EpochArrays are supported")

    return np.array(find_intersecting_intervals_(set1.data, set2.data))


def get_velocity(position, time=None):
    if time is None:
        time = np.arange(position.shape[0])
    return np.gradient(position, time, axis=0)


def get_speed(position, time=None):
    velocity = get_velocity(position, time=time)
    return np.sqrt(np.sum(velocity**2, axis=1))

def find_interval(logical):
    """
    Find consecutive intervals of True values in a list of boolean values.
    
    Parameters:
    logical (List[bool]): The list of boolean values.
    
    Returns:
    List[Tuple[int, int]]: A list of tuples representing the start and end indices of each consecutive interval of True values in the logical list.
    
    Example:
    find_interval([0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1]) -> [(2, 4), (6, 7), (10, 11)]
    find_interval([1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1]) -> [(0, 2), (4, 5), (9, 10)]
    """
    intervals = []
    start = None
    for i, value in enumerate(logical):
        if value and start is None:
            start = i
        elif not value and start is not None:
            intervals.append((start, i - 1))
            start = None
    if start is not None:
        intervals.append((start, len(logical) - 1))
    return intervals

@njit(parallel=True)
def in_intervals(timestamps, intervals):
    """
    Find which timestamps fall within the given intervals.

    Parameters
    ----------
    timestamps : ndarray
        An array of timestamp values. assumes sorted
    intervals : ndarray
        An array of time intervals, represented as pairs of start and end times.

    Returns
    -------
    ndarray
        A logical index indicating which timestamps fall within the intervals.

    Examples
    --------
    >>> timestamps = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    >>> intervals = np.array([[2, 4], [5, 7]])
    >>> in_intervals(timestamps, intervals)
    array([False,  True,  True,  True,  True,  True,  True, False])
    """
    in_interval = np.empty(timestamps.shape, dtype=np.bool_) * False
    for i, (start, end) in enumerate(intervals):
        in_interval[(timestamps >= start) & (timestamps <= end)] = True

    return in_interval


def count_events(events, time_ref, time_range):
    """
    Count the number of events that occur within a given time range after each reference event.
    Parameters
    ----------
    events : ndarray
        A 1D array of event times.
    time_ref : ndarray
        A 1D array of reference times.
    time_range : tuple
        A tuple containing the start and end times of the time range.
    Returns
    -------
    counts : ndarray
        A 1D array of event counts, one for each reference time (same len as time_ref).
    """
    # Initialize an array to store the event counts
    counts = np.zeros_like(time_ref)

    # Iterate over the reference times
    for i, r in enumerate(time_ref):
        # Check if any events occur within the time range
        idx = (events > r + time_range[0]) & (events < r + time_range[1])
        # Increment the event count if any events are found
        counts[i] = len(events[idx])

    return counts


def clean_lfp(lfp, thresholds=[5,10], artifact_time_expand=[0.25,0.1]):
    """
    Remove artefacts and noise from a local field potential (LFP) signal.

    Parameters
    ----------
    lfp : nelpy AnalogSignalArray
        The LFP signal to be cleaned. Single signal only
    thresholds : list, optional
        A list of two thresholds for detecting artefacts and noise. The first threshold is used to detect large global
        artefacts by finding values in the z-scored LFP signal that deviate by more than the threshold number of sigmas
        from the mean. The second threshold is used to detect noise by finding values in the derivative of the z-scored
        LFP signal that are greater than the threshold. Default is [5, 10].
    artifact_time_expand : list, optional
        A list of two time intervals around detected artefacts and noise. The first interval is used to expand the detected
        large global artefacts. The second interval is used to expand the detected noise. Default is [0.25, 0.1].

    Returns
    -------
    ndarray
        The cleaned LFP signal.

    Based on https://github.com/ayalab1/neurocode/blob/master/lfp/CleanLFP.m by Ralitsa Todorova

    Examples
    --------
    >>> lfp = nel.AnalogSignalArray(data=np.random.randn(1250),timestamps=np.arange(1250)/1250)
    >>> clean_lfp(lfp)
    array([-1.73104885,  1.08192036,  1.40332741, ..., -2.78671212,
        -1.63661574, -1.10868426])
    """
    threshold1 = thresholds[0]  # in sigmas deviating from the mean
    aroundArtefact1 = artifact_time_expand[0]  # interval to expand large global artefacts

    threshold2 = thresholds[1]  # for derivative of z-scored signal
    aroundArtefact2 = artifact_time_expand[1]  # interval to expand detected noise

    t = lfp.time  # time points of LFP signal
    values = lfp.copy().data.flatten()  # values of LFP signal
    z = lfp.zscore().data.flatten()  # z-scored values of LFP signal
    d = np.append(np.diff(z), 0)  # derivative of z-scored LFP signal

    # Detect large global artefacts [0]
    artefactInterval = t[np.array(find_interval(np.abs(z) > threshold1), dtype=int)]
    artefactInterval = nel.EpochArray(artefactInterval)
    if not artefactInterval.isempty:
        artefactInterval = artefactInterval.expand(aroundArtefact1)

    # Find noise using the derivative of the z-scored signal [1]
    noisyInterval = t[np.array(find_interval(np.abs(d) > threshold2), dtype=int)]
    noisyInterval = nel.EpochArray(noisyInterval)
    if not noisyInterval.isempty:
        noisyInterval = noisyInterval.expand(aroundArtefact2)

    # Combine intervals for artefacts and noise
    bad = (artefactInterval | noisyInterval).merge()

    if bad.isempty:
        return values

    # Find timestamps within intervals for artefacts and noise
    in_interval = in_intervals(t, bad.data)

    # Interpolate values for timestamps within intervals for artefacts and noise
    values[in_interval] = np.interp(t[in_interval], t[~in_interval], values[~in_interval])
    
    return values

