import copy
import sys
sys.path.append(r'D:\github\ripple_heterogeneity\python')
import functions,loading

import pandas as pd
import numpy as np 
import glob
import os
import pickle
import seaborn as sns
import nelpy as nel
from sklearn.decomposition import PCA
import multiprocessing
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression

# @jit(nopython=True)
def SVCA_(X):
    """
    computes a cross-validated form of PCA: Shared Variance Component
    Analysis. Components are extracted from the covariance between neuron
    sets "ntrain" and "ntest" on training timepts "itrain". The variance of
    these covariance components is computed on testing timepts "itest". 
    This variance is the amount of reliable variance of that component 
    (because it's consistent across timepts).
    Inputs:
        X: (neurons x timepts)
    Outputs:
        scov (shared variance of each covariance component)
        varcov (total variance of each covariance component)

    From: https://github.com/MouseLand/stringer-et-al-2019/blob/master/utils.py    
    """
    
    NN,NT = X.shape

    # split cells into test and train
    norder = np.random.permutation(NN)
    nhalf = int(norder.size/2)
    ntrain = norder[:nhalf]
    ntest = norder[nhalf:]

    # split time into test and train
    torder = np.random.permutation(NT)
    thalf = int(torder.size/2)
    ttrain = torder[:thalf]
    ttest = torder[thalf:]

    cov = X[np.ix_(ntrain, ttrain)] @ X[np.ix_(ntest, ttrain)].T
    u = PCA(n_components=min(1024, nhalf-1), svd_solver='randomized').fit_transform(cov)
    u /= (u**2).sum(axis=0)**0.5
    v = cov.T @ u
    v /= (v**2).sum(axis=0)**0.5

    strain = u.T @ X[np.ix_(ntrain,ttest)]
    stest = v.T @ X[np.ix_(ntest,ttest)]

    # covariance k is uk.T * F * G.T * vk / npts
    scov = (strain * stest).mean(axis=1)
    varcov = (strain**2 + stest**2).mean(axis=1) / 2

    return scov,varcov,strain,stest

def SVCA(X,folds=10):
    """ calls SVCA_ """

    #  do multiple subsets and avg
    scov = []
    varcov = []
    strain = []
    stest = []
    svc_neur = []
    for i in range(folds):
        scov_,varcov_,strain_,stest_ = SVCA_(X)
        scov.append(scov_)
        varcov.append(varcov_)
        strain.append(strain_)
        stest.append(stest_)
        svc_neur.append(scov_ / varcov_)

    return svc_neur,scov,varcov,strain,stest

def main_analysis(unit_mat,beh_epochs,epoch_df,nrem_epochs,wake_epochs,restrict_to_state=True):

  # iter through each behavioral epoch
  scov = []
  varcov = []
  strain = []
  stest = []
  svc_neur = []
  for i,ep in enumerate(beh_epochs):
    if restrict_to_state:
        if epoch_df.environment.iloc[i] == 'sleep':
            temp_st = unit_mat[nrem_epochs][ep]
        else:
            temp_st = unit_mat[wake_epochs][ep]
    else:
        temp_st = unit_mat[ep]
            
    svc_neur_,scov_,varcov_,strain_,stest_ = SVCA(temp_st.data)
    scov.append(scov_)
    varcov.append(varcov_)
    strain.append(strain_)
    stest.append(stest_)
    svc_neur.append(svc_neur_)

  results = {}
  results['scov'] = scov
  results['varcov'] = varcov
  results['strain'] = strain
  results['stest'] = stest
  results['svc_neur'] = svc_neur

  return results

def load_needed_data(basepath):
    """ gets and formats basic data"""

    cell_metrics,data,ripples,fs_dat = loading.load_basic_data(basepath)

    restrict_idx = ((cell_metrics.putativeCellType == "Pyramidal Cell") &
                        ((cell_metrics.brainRegion=="CA1") |
                        (cell_metrics.brainRegion=="rCA1") |
                        (cell_metrics.brainRegion=="lCA1")) &
                        (cell_metrics.bad_unit==False))

    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]

    # get spike train array
    try:
        st = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx], fs=fs_dat)
    except: # if only single cell... should prob just skip session
        st = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx][0], fs=fs_dat)

    # behavioral epochs
    epoch_df = loading.load_epoch(basepath)
    behavioral_epochs = nel.EpochArray([np.array([epoch_df.startTime,
                                                    epoch_df.stopTime]).T])

    # get brain states                                                
    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict['NREMstate'])
    wake_epochs = nel.EpochArray(state_dict['WAKEstate'])

    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])

    unit_mat = functions.get_participation(st.data,ripples.start,ripples.stop,par_type='firing_rate')
    unit_mat = nel.AnalogSignalArray(data=unit_mat,timestamps=ripples.peaks,support=ripple_epochs)

    return cell_metrics,st,epoch_df,behavioral_epochs,nrem_epochs,wake_epochs,ripple_epochs,unit_mat

def pooled_incoherent_shuffle(X):
    """Incoherent shuffle on X, circ shifting rows."""
    data = copy.deepcopy(X)
    out = data.data
    for uu in range(X.n_signals):
        segment = out[uu,:]
        segment = np.roll(segment, np.random.randint(len(segment)))
        out[uu,:] = segment
    data.data = out
    return data

def estimate_slope(svc_neur):
    """ 
    get slope from svc_neur
    x var in regression is in log10
    Input: svc_neur nested list, each level is an epoch
    Output: slope, intercept, r2
    """
    slope, intercept, r2 = [],[],[]

    for i in range(len(svc_neur)):

        y = np.array(svc_neur[i]).mean(axis=0)
        x = np.log10(np.arange(len(y)))

        bad_idx = np.isnan(x) | np.isinf(x)
        y = y[~bad_idx]
        x = x[~bad_idx]
        x = x[:,np.newaxis]

        reg = LinearRegression().fit(x, y)

        slope.append(reg.coef_)
        intercept.append(reg.intercept_)
        r2.append(reg.score(x, y)**2)

    return slope, intercept, r2

def set_up_and_do_analysis(basepath,n_shuffles=100):

    (cell_metrics,
    st,
    epoch_df,
    behavioral_epochs,
    nrem_epochs,
    wake_epochs,
    ripple_epochs,
    unit_mat) = load_needed_data(basepath)

    if cell_metrics.shape[0] == 0:
        return     

    # run main analysis and get SVC
    results = main_analysis(unit_mat,behavioral_epochs,epoch_df,nrem_epochs,wake_epochs)

    slope, intercept, r2 = estimate_slope(results['svc_neur'])
    results['slope'] = slope
    results['intercept'] = intercept
    results['r2'] = r2

    results['cell_metrics'] = cell_metrics
    results['basepath'] = basepath

    # make null dist from shuffling across 
    svc_neur_shuff = []
    for _ in range(n_shuffles):
        scov_,varcov_,_,_ = SVCA_(pooled_incoherent_shuffle(unit_mat).data)
        svc_neur_shuff.append(scov_ / varcov_)

    results['svc_neur_shuff'] = svc_neur_shuff

    slope, intercept, r2 = estimate_slope([results['svc_neur_shuff']])
    results['slope_shuff'] = slope
    results['intercept_shuff'] = intercept
    results['r2_shuff'] = r2

    return results

def session_loop(basepath,save_path):

    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.pkl')
    if os.path.exists(save_file):
        return

    results = set_up_and_do_analysis(basepath)

    # save file
    with open(save_file, 'wb') as f:
        pickle.dump(results, f)

def main_run(df,save_path,parallel=True):
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()         
        processed_list = Parallel(n_jobs=num_cores)(delayed(session_loop)(basepath,save_path) for basepath in basepaths)
    else:    
        for basepath in basepaths:
            print(basepath)
            session_loop(basepath,save_path)   