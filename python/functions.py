import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from numba import jit


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