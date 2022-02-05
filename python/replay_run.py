import numpy as np

import nelpy as nel
import nelpy.plotting as npl
from nelpy.analysis import replay

import os
import matplotlib.pyplot as plt

import pandas as pd
import statistics 
from scipy import stats

import multiprocessing
from joblib import Parallel, delayed

import statsmodels.api as sm
import pickle
import copy
import loading


def decode_and_score(bst, tc, pos):
    # access decoding accuracy on behavioral time scale 
    posteriors, lengths, mode_pth, mean_pth = nel.decoding.decode1D(bst,
                                                                    tc,
                                                                    xmin=np.nanmin(pos.data),
                                                                    xmax=np.nanmax(pos.data))
    actual_pos = pos(bst.bin_centers)
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(actual_pos, mode_pth)
    median_error = np.nanmedian(np.abs(actual_pos - mode_pth))
    
    return rvalue,median_error

def pooled_incoherent_shuffle_bst(bst):
    out = copy.deepcopy(bst)
    data = out._data

    for uu in range(bst.n_units):
        segment = np.atleast_1d(np.squeeze(data[uu, :]))
        segment = np.roll(segment, np.random.randint(len(segment)))
        data[uu, :] = segment
    return out

def decode_and_shuff(bst, tc, pos, n_shuffles=500):
    """
    """
    rvalue, median_error = decode_and_score(bst, tc, pos)

    scores = np.zeros(bst.n_epochs)
    if n_shuffles > 0:
        rvalue_time_swap = np.zeros((n_shuffles,1))
        median_error_time_swap = np.zeros((n_shuffles,1))
        
    for shflidx in range(n_shuffles):
        bst_shuff = pooled_incoherent_shuffle_bst(bst)
        rvalue_time_swap[shflidx], median_error_time_swap[shflidx] = decode_and_score(bst_shuff, tc, pos)
            
    return rvalue, median_error, rvalue_time_swap, median_error_time_swap

def _m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)


def _cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - _m(x, w)) * (y - _m(y, w))) / np.sum(w)


def _corr(x, y, w):
    """Weighted Correlation"""
    return _cov(x, y, w) / np.sqrt(_cov(x, x, w) * _cov(y, y, w))

def weighted_correlation(posterior, time, place_bin_centers):
    """ From Eric Denovellis """
    place_bin_centers = place_bin_centers.squeeze()
    posterior[np.isnan(posterior)] = 0.0

    return _corr(time[:, np.newaxis],
                 place_bin_centers[np.newaxis, :], posterior)

def score_array(posterior):
    """
    takes in posterior matrix (distance by time) and conducts
    weighted least squares
    """
    nan_loc = np.isnan(posterior).any(axis=0)

    rows, cols = posterior.shape

    x = np.arange(cols)
    y = posterior.argmax(axis=0)
    w = posterior.max(axis=0)

    x = x[~nan_loc]
    y = y[~nan_loc]
    w = w[~nan_loc]
    
    # if only one time bin is active
    if len(x)==1:
        return np.nan,np.nan,np.nan,np.nan
    
    X = sm.add_constant(x)
    wls_model = sm.WLS(y,X, weights=w)
    results = wls_model.fit()
    
    slope = results.params[1]
    intercept = results.params[0]
    log_like = wls_model.loglike(results.params)

    return results.rsquared,slope,intercept,log_like

def get_score_coef(bst,bdries,posterior):
    """
    runs score_array on each event epoch in bst (binned spike train)
    """
    scores = np.zeros(bst.n_epochs)
    slope = np.zeros(bst.n_epochs)
    intercept = np.zeros(bst.n_epochs)
    log_like = np.zeros(bst.n_epochs)

    for idx in range(bst.n_epochs):
        posterior_array = posterior[:, bdries[idx]:bdries[idx+1]]
        scores[idx],slope[idx],intercept[idx],log_like[idx] = score_array(posterior_array)
    return scores,slope,intercept,log_like


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

def shuff(posterior_array,time,place_bin_centers,dt,dp):
    
    posterior_ts = replay.time_swap_array(posterior_array)
    posterior_cs = replay.column_cycle_array(posterior_array)


    w_corr_time_swap = weighted_correlation(posterior_ts, time, place_bin_centers)
    w_corr_col_cycle = weighted_correlation(posterior_cs, time, place_bin_centers)
    
    return w_corr_time_swap,w_corr_col_cycle

def get_features(bst_placecells,
                 posteriors,
                 bdries,
                 mode_pth,
                 pos,
                 figs=False,
                 dp=3):
    """
    Using the posterior probability matrix, calculate several features on spatial trajectory
    and detects if the trajectory is foward or reverse depending on the rat's current position
    """
    
    place_bin_edges = np.arange(np.nanmin(pos.data[0]), np.nanmax(pos.data[0]) + dp, dp)
    place_bin_centers = place_bin_edges[:-1] + np.diff(place_bin_edges) / 2
    
    traj_dist = []
    traj_speed = []
    traj_step = []
    replay_type = []
    dist_rat_start = []
    dist_rat_end = []
    position = []

    for idx in range(bst_placecells.n_epochs):

        x = bst_placecells[idx].bin_centers

        y = mode_pth[bdries[idx]:bdries[idx+1]]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        velocity, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
        y = x*velocity+intercept
        position.append(y)

        # get spatial difference between bins
        dy = np.abs(np.diff(y))
        # get cumulative distance 
        traj_dist.append(np.nansum(dy))
        # calculate avg speed of trajectory (dist(cm) / time(sec))
        traj_speed.append(np.nansum(dy) / (np.nanmax(x) - np.nanmin(x)))
        # get mean step size 
        traj_step.append(np.nanmean(dy))
        
        # check if current event is outside beh epoch
        if all(x.max() < pos.abscissa_vals) | all(x.min() > pos.abscissa_vals):
            replay_type.append(np.nan)
            dist_rat_start.append(np.nan)
            dist_rat_end.append(np.nan)
            continue

        rat_event_pos = np.interp(x,pos.abscissa_vals,pos.data[0])
        rat_x_position = np.nanmean(rat_event_pos)

        # get dist of the start & end of trajectory to rat
        dist_rat_start.append(rat_x_position - y[0])
        dist_rat_end.append(rat_x_position - y[-1])

        # what side of the track is the rat on ? 
        min_max_env = [np.nanmin(pos.data[0]),np.nanmax(pos.data[0])]
        side = np.argmin(np.abs(min_max_env- rat_x_position))

        if (side == 1) & (velocity < 0):
            replay_type.append('forward')
        elif (side == 1) & (velocity > 0):
            replay_type.append('reverse')
        elif (side == 0) & (velocity < 0):
            replay_type.append('reverse')
        elif (side == 0) & (velocity > 0):
            replay_type.append('forward')
        else:
            replay_type.append(np.nan)

        if figs:
            fig = plt.figure(figsize=(4,3))
            ax = plt.gca()
            npl.plot(x,rat_event_pos,"^",color='brown',linewidth=10,ax=ax)
            ax.plot(x,y,'k',linewidth=2)
            ax.scatter(x[0],y[0],color='g')
            ax.scatter(x[-1],y[-1],color='r')
            ax.set_title(replay_type[idx])

    return traj_dist,traj_speed,traj_step,replay_type,dist_rat_start,dist_rat_end,position

def handle_behavior(basepath,epoch_df):
    beh_df = loading.load_animal_behavior(basepath)

    beh_df = beh_df[~np.isnan(beh_df.linearized)]

    pos = nel.AnalogSignalArray(data=np.array(beh_df.linearized),
                                timestamps=beh_df.time,
                                fs=1/statistics.mode(np.diff(beh_df.time)))
    return pos

def run_all(
    basepath, # basepath to session
    traj_shuff=1500, # number of shuffles to determine sig replay
    ds_run=0.5,
    ds_50ms=0.05,
    s_binsize=3, # spatial bins in tuning curve
    speed_thres=4, # running threshold to determine tuning curves
    min_rip_dur=0.08 # min ripple duration for replay
):
    """
    Main function that conducts the replay analysis
    """

    # maze_size_cm,pos,st_all = get_base_data(data_path,spike_path,session)
    cell_metrics,data,ripples,fs_dat = loading.load_basic_data(basepath)

    restrict_idx = ((cell_metrics.putativeCellType == "Pyramidal Cell") &
                        ((cell_metrics.brainRegion=="CA1") |
                        (cell_metrics.brainRegion=="rCA1") |
                        (cell_metrics.brainRegion=="lCA1")) &
                        (cell_metrics.bad_unit==False) & 
                        (cell_metrics.total >= 100))
    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]
    
    if cell_metrics.shape[0] == 0:
        return
    try:
        st_all = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx], fs=fs_dat)
    except:
        st_all = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx][0], fs=fs_dat)

    epoch_df = loading.load_epoch(basepath)
    # remove sleep and wheel running
    epoch_df = epoch_df[(epoch_df.environment != 'sleep') & (epoch_df.environment != 'wheel')]
    # remove sessions < 5 minutes
    epoch_df = epoch_df[(epoch_df.stopTime - epoch_df.startTime)/60 > 5]
    beh_epochs = nel.EpochArray([np.array([epoch_df.startTime,epoch_df.stopTime]).T])

    pos = handle_behavior(basepath,epoch_df)
    pos = pos[beh_epochs[0]]

    # compute and smooth speed
    speed1 = nel.utils.ddt_asa(pos, smooth=True, sigma=0.1, norm=True)
 
    # find epochs where the animal ran > 4cm/sec
    run_epochs = nel.utils.get_run_epochs(speed1, v1=speed_thres, v2=speed_thres)

    # restrict spike trains to those epochs during which the animal was running
    st_run = st_all[beh_epochs[0]][run_epochs] 
    
    # smooth and re-bin:
    # 300 ms spike smoothing
    bst_run = st_run.bin(ds=ds_50ms).smooth(sigma=0.3 , inplace=True).rebin(w=ds_run/ds_50ms)

    sigma = 3 # smoothing std dev in cm
    tc = nel.TuningCurve1D(bst=bst_run, 
                            extern=pos,
                            n_extern=int((np.nanmax(pos.data)-np.nanmin(pos.data))/s_binsize),
                            extmin=np.nanmin(pos.data),
                            extmax=np.nanmax(pos.data),
                            sigma=sigma,
                            min_duration=0)

    # locate pyr cells with >= 100 spikes, peak rate >= 1 Hz, peak/mean ratio >=1.5
    peak_firing_rates = tc.max(axis=1)
    mean_firing_rates = tc.mean(axis=1)
    ratio = peak_firing_rates/mean_firing_rates
  
    idx = (st_run.n_events >= 100) & (tc.ratemap.max(axis=1) >= 1) & (ratio>=1.5)
    unit_ids_to_keep = (np.where(idx)[0]+1).squeeze().tolist()

    sta_placecells = st_all._unit_subset(unit_ids_to_keep)
    tc = tc._unit_subset(unit_ids_to_keep)
    total_units = sta_placecells.n_active
    
    # access decoding accuracy on behavioral time scale 
    decoding_r2, median_error, decoding_r2_shuff, _ = decode_and_shuff(bst_run.loc[:,unit_ids_to_keep],
                                                                        tc,
                                                                        pos,
                                                                        n_shuffles=1000)
    # check decoding quality against chance distribution
    _, decoding_r2_pval = get_significant_events(decoding_r2, decoding_r2_shuff)
    
    # get ready to decode replay
    # restrict to events at least xxms
    ripples = ripples[ripples.duration >= min_rip_dur]
    
    # make epoch object
    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])

    # bin data into 20ms 
    bst_placecells = sta_placecells[ripple_epochs].bin(ds=0.02)

    # count units per event
    n_active = [bst.n_active for bst in bst_placecells]
    n_active = np.array(n_active) 
    # also count the proportion of bins in each event with 0 activity
    inactive_bin_prop = [sum(bst.n_active_per_bin == 0) / bst.lengths[0] for bst in bst_placecells]
    inactive_bin_prop = np.array(inactive_bin_prop) 
    # restrict bst to instances with >= 5 active units and < 50% inactive bins
    idx = (n_active >= 5) & (inactive_bin_prop < .5)
    bst_placecells = bst_placecells[np.where(idx)[0]]
    # restrict to instances with >= 5 active units
    n_active = n_active[idx]
    inactive_bin_prop = inactive_bin_prop[idx]

    # decode each ripple event
    posteriors, bdries, mode_pth, mean_pth = nel.decoding.decode1D(bst_placecells,
                                                                    tc,
                                                                    xmin=np.nanmin(pos.data),
                                                                    xmax=np.nanmax(pos.data))
            
    # score each event using trajectory_score_bst (sums the posterior probability in a range (w) from the LS line)
    scores, scores_time_swap, scores_col_cycle = replay.trajectory_score_bst(bst_placecells,
                                                                                tc,
                                                                                w=3,
                                                                                n_shuffles=traj_shuff,
                                                                                normalize=True)
    
    # find sig events using time and column shuffle distributions
    _,score_pval_time_swap = get_significant_events(scores, scores_time_swap)
    _,score_pval_col_cycle = get_significant_events(scores, scores_col_cycle)

    (
        traj_dist,
        traj_speed,
        traj_step,
        replay_type,
        dist_rat_start,
        dist_rat_end,
        position
    ) = get_features(bst_placecells,posteriors,bdries,mode_pth,pos)
    
    slope, intercept, r2values = replay.linregress_bst(bst_placecells, tc)

    # package data into results dictionary
    results = {}

    results['sta_placecells'] = sta_placecells
    results['bst_placecells'] = bst_placecells
    results['tc'] = tc
    results['posteriors'] = posteriors
    results['bdries'] = bdries
    results['mode_pth'] = mode_pth
    results['position'] = position
    
    temp_df = pd.DataFrame()
    # add event by event metrics to df
    temp_df['n_active'] = n_active
    temp_df['inactive_bin_prop'] = inactive_bin_prop
    temp_df['trajectory_score'] = scores
    temp_df['r_squared'] = r2values
    temp_df['slope'] = slope
    temp_df['intercept'] = intercept
    temp_df['score_pval_time_swap'] = score_pval_time_swap
    temp_df['score_pval_col_cycle'] = score_pval_col_cycle
    temp_df['traj_dist'] = traj_dist
    temp_df['traj_speed'] = traj_speed
    temp_df['traj_step'] = traj_step
    temp_df['replay_type'] = replay_type
    temp_df['dist_rat_start'] = dist_rat_start
    temp_df['dist_rat_end'] = dist_rat_end
    results['df'] = temp_df

    results['session'] = basepath
    results['decoding_r2'] = decoding_r2
    results['decoding_r2_pval'] = decoding_r2_pval
    results['decoding_median_error'] = median_error
    results['total_units'] = total_units

    return results

def main_loop(basepath,save_path):
    '''
    main_loop: file management 
    '''
    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.pkl')
    if os.path.exists(save_file):
        return
        
    # calc some features
    results = run_all(basepath)
    # save file
    with open(save_file, 'wb') as f:
        pickle.dump(results, f)

def main(df,save_path,parallel=True):
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()         
        processed_list = Parallel(n_jobs=num_cores)(delayed(main_loop)(basepath,save_path) for basepath in basepaths)
    else:    
        for basepath in basepaths:
            print(basepath)
            main_loop(basepath,save_path)   