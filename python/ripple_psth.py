import nelpy as nel
import functions,loading
import sys
sys.path.append(r'D:\github\neurocode\reactivation\assemblies')
import numpy as np
import pickle
import pandas as pd
import os
import multiprocessing
from joblib import Parallel, delayed


def main_analysis(st,ripples,behavioral_epochs,epoch_df,nrem_epochs,wake_epochs):

    # iter through each behavioral epoch
    ccg = []

    for i,ep in enumerate(behavioral_epochs):
        if epoch_df.environment.iloc[i] == 'sleep':
            temp_st = st[nrem_epochs][ep]
        else:
            temp_st = st[wake_epochs][ep]

        ccg.append(functions.compute_psth(temp_st.data,ripples.peaks.values))
    
    return ccg

def session_loop(basepath,save_path):

    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.pkl')
    if os.path.exists(save_file):
        return

    cell_metrics,data,ripples,fs_dat = loading.load_basic_data(basepath)

    restrict_idx = ((cell_metrics.putativeCellType == "Pyramidal Cell") &
                        ((cell_metrics.brainRegion=="CA1") |
                        (cell_metrics.brainRegion=="rCA1") |
                        (cell_metrics.brainRegion=="lCA1")) &
                        (cell_metrics.bad_unit==False))

    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]

    if cell_metrics.shape[0] == 0:
        return
        
    # get ripple epochs
    # ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])

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

    ccg = main_analysis(st,ripples,behavioral_epochs,epoch_df,nrem_epochs,wake_epochs)

    results = {}
    results['ccg'] = ccg
    results['UID'] = cell_metrics.UID
    results['basepath'] = basepath
    results['deepSuperficial'] = cell_metrics.deepSuperficial

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