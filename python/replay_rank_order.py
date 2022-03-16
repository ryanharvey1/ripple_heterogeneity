import sys
sys.path.append(r'D:\github\ripple_heterogeneity\python')
import functions,loading
import pandas as pd
import numpy as np 
import glob
import os
import pickle
import nelpy as nel
import warnings
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="fs was not specified, so we try to estimate it from the data...")

import multiprocessing
from joblib import Parallel, delayed


def get_replay_epochs(results,replay_type,direction,alpha=0.05):

    if not np.any((replay_type == 'forward') | (replay_type == 'reverse')):
        raise Exception('wrong replay_type')

    if not np.any((direction == 'outbound_epochs') | (direction == 'inbound_epochs')):
        raise Exception('wrong direction')

    idx = (results[direction]['df'].score_pval_time_swap <= alpha) & (results[direction]['df'].replay_type == replay_type)
    starts = results[direction]['df'][idx].start.values
    stops = results[direction]['df'][idx].stop.values

    sort_idx = np.argsort(starts)

    starts = starts[sort_idx]
    stops = stops[sort_idx]

    return nel.EpochArray(np.array([starts,stops]).T)
    
def get_pre_linear_post(basepath):
    epoch_df = loading.load_epoch(basepath)
    pattern_idx,_ = functions.find_epoch_pattern(epoch_df.environment,['sleep','linear','sleep'])
    epoch_df = epoch_df[pattern_idx]
    return epoch_df,nel.EpochArray([np.array([epoch_df.startTime,epoch_df.stopTime]).T])


def run_all(basepath):
    # for session in sessions:
    with open(basepath, 'rb') as f:
        results = pickle.load(f)

    basepath = results['outbound_epochs']['session']

    st,cell_metrics = loading.load_spikes(basepath,putativeCellType='Pyr',brainRegion='CA1')

    epochs = {
            'forward_outbound_replay': get_replay_epochs(results,'forward','outbound_epochs'),
            'forward_inbound_replay': get_replay_epochs(results,'forward','inbound_epochs'),
            'reverse_outbound_replay': get_replay_epochs(results,'reverse','outbound_epochs'),
            'reverse_inbound_replay': get_replay_epochs(results,'reverse','inbound_epochs')
            }

    epoch_df,behavior_epochs = get_pre_linear_post(basepath)

    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])

    median_rank_order = []
    label = []
    deepSuperficial = []
    uid = []
    environment = []
    avg_fr = []
    particip = []

    par_mat = functions.get_participation(st.data,ripple_epochs.starts,ripple_epochs.stops)
    par_mat = nel.AnalogSignalArray(data=par_mat,timestamps=ripple_epochs.centers)

    for beh_i,beh_epoch in enumerate(behavior_epochs):
        for key_ in  epochs.keys():
            temp_rank_order,rank_order = functions.get_rank_order(st[beh_epoch],epochs[key_])
            avg_fr.append(st[beh_epoch].n_spikes / beh_epoch.duration)
            particip.append(par_mat[beh_epoch].mean(axis=1))
            median_rank_order.append(temp_rank_order)
            label.append([key_]*len(temp_rank_order))
            deepSuperficial.append(cell_metrics.deepSuperficial)
            uid.append(cell_metrics.UID)
            environment.append([beh_i]*len(temp_rank_order))
    
    df_rank_order = pd.DataFrame()
    df_rank_order['median_rank_order'] = np.hstack(median_rank_order)
    df_rank_order['avg_fr'] = np.hstack(avg_fr)
    df_rank_order['particip'] = np.hstack(particip)
    df_rank_order['label'] = np.hstack(label)
    df_rank_order['deepSuperficial'] = np.hstack(deepSuperficial)
    df_rank_order['UID'] = np.hstack(uid)
    df_rank_order['environment'] = np.hstack(environment)
    df_rank_order['basepath'] = basepath

    df_rank_order['label_lay'] = df_rank_order['label'] + df_rank_order['deepSuperficial']

    return df_rank_order

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

def load_results(save_path):
    sessions = glob.glob(save_path +os.sep+ '*.pkl')
    df = pd.DataFrame()
    for session in sessions:
        with open(session, 'rb') as f:
            results = pickle.load(f)
        df = pd.concat([df,results],ignore_index=True)
    return df