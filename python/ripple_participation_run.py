import pandas as pd
import numpy as np 
import functions,loading
import nelpy as nel  # main nelpy imports
import os
import multiprocessing
from joblib import Parallel, delayed
import compress_repeated_epochs

    
def get_participation(st,ripple_epochs):
    # get participation prob.
    # make matrix n rows (units) by n cols (ripple epochs)
    unit_mat = np.zeros((st.n_units,ripple_epochs.n_intervals))
    for i,event in enumerate(st):
        unit_mat[:,i] = (event.n_events>0)*1
    return unit_mat

def get_epoched_values(st_unit,st_unit_rip,st_unit_no_rip,ep,ripple_epochs,state_epoch=False):
    
    avg_fr_not_rip = np.nan
    avg_fr_in_rip = np.nan
    participation_prob = np.nan
    n_spikes = np.nan

    if state_epoch == False:
        unit_mat = get_participation(st_unit_rip[ep],
                                        ripple_epochs[ep])
        participation_prob = np.sum(unit_mat,axis=1) / unit_mat.shape[1]
        if len(participation_prob) == 0:
            participation_prob = np.nan
        try:
            avg_fr_not_rip = st_unit_no_rip[ep].n_spikes / ep[~ripple_epochs].duration
        except:
            pass
        try:
            avg_fr_in_rip = st_unit_rip[ep].n_spikes / st_unit_rip[ep].support.duration
        except:
            pass
        try:
            n_spikes = st_unit[ep].n_spikes
        except:
            pass
    else:
        try:
            unit_mat = get_participation(st_unit_rip[ep][state_epoch],
                                            ripple_epochs[ep][state_epoch])
            participation_prob = np.sum(unit_mat,axis=1) / unit_mat.shape[1]
            if len(participation_prob) == 0:
                participation_prob = np.nan
        except:
            pass
        try:
            avg_fr_not_rip = (st_unit_no_rip[ep][state_epoch].n_spikes /
                                ep[~ripple_epochs][state_epoch].duration)
        except:
            pass
        try:    
            avg_fr_in_rip = (st_unit_rip[ep][state_epoch].n_spikes /
                                st_unit_rip[ep][state_epoch].support.duration)
        except:
            pass
        try:
            n_spikes = st_unit[ep][state_epoch].n_spikes
        except:
            pass

    return participation_prob,avg_fr_not_rip,avg_fr_in_rip,n_spikes

def main_analysis(
                    basepath,
                    cell_metrics,
                    st_unit,
                    ripple_epochs,
                    behavioral_epochs,
                    epoch_df,
                    nrem_epochs,
                    wake_epochs
                    ):

    # create empty dataframe to add metrics from each epoch
    df_save = pd.DataFrame()

    # create spike train object with ripple epochs
    st_unit_rip = st_unit[ripple_epochs]
    st_unit_no_rip = st_unit[~ripple_epochs]

    epochs = epoch_df.environment
    familiarity = epoch_df.behavioralParadigm

    for i,ep in enumerate(behavioral_epochs):
        temp_save = pd.DataFrame()

        (
        participation_prob,
        avg_fr_not_rip,
        avg_fr_in_rip,
        n_spikes
        ) = get_epoched_values(st_unit,st_unit_rip,st_unit_no_rip,ep,ripple_epochs)

        (
        participation_prob_nrem,
        avg_fr_not_rip_nrem,
        avg_fr_in_rip_nrem,
        n_spikes_nrem
        ) = get_epoched_values(st_unit,st_unit_rip,st_unit_no_rip,ep,ripple_epochs,nrem_epochs)

        (
        participation_prob_wake,
        avg_fr_not_rip_wake,
        avg_fr_in_rip_wake,
        n_spikes_wake
        ) = get_epoched_values(st_unit,st_unit_rip,st_unit_no_rip,ep,ripple_epochs,wake_epochs)

        # package results and previously saved metrics
        temp_save['UID'] = cell_metrics.UID
        temp_save['basepath'] = basepath
        temp_save['epoch'] = epochs[i]
        temp_save['epoch_n'] = i
        temp_save['familiarity'] = familiarity[i]
        temp_save['deepSuperficial'] = cell_metrics.deepSuperficial
        temp_save['brainRegion'] = cell_metrics.brainRegion
        temp_save['putativeCellType'] = cell_metrics.putativeCellType
        temp_save['participation_prob'] = participation_prob
        temp_save['participation_prob_nrem'] = participation_prob_nrem
        temp_save['participation_prob_wake'] = participation_prob_wake
        temp_save['avg_fr_not_rip'] = avg_fr_not_rip
        temp_save['avg_fr_in_rip'] = avg_fr_in_rip
        temp_save['n_spikes'] = n_spikes
        temp_save['avg_fr_not_rip_nrem'] = avg_fr_not_rip_nrem
        temp_save['avg_fr_in_rip_nrem'] = avg_fr_in_rip_nrem
        temp_save['n_spikes_nrem'] = n_spikes_nrem
        temp_save['avg_fr_not_rip_wake'] = avg_fr_not_rip_wake
        temp_save['avg_fr_in_rip_wake'] = avg_fr_in_rip_wake
        temp_save['n_spikes_wake'] = n_spikes_wake
        temp_save['session_dur'] = ep.duration
        temp_save['session_dur_nrem'] = ep[nrem_epochs].duration
        temp_save['session_dur_wake'] = ep[wake_epochs].duration
        temp_save['n_ripples'] = ripple_epochs[ep].n_intervals
        temp_save['n_ripples_nrem'] = ripple_epochs[ep][nrem_epochs].n_intervals
        temp_save['n_ripples_wake'] = ripple_epochs[ep][wake_epochs].n_intervals

        df_save = df_save.append(temp_save,ignore_index=True)

    return df_save

def session_loop(basepath,save_path):
    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.csv')
    if os.path.exists(save_file):
        return

    cell_metrics,data,ripples,fs_dat = loading.load_basic_data(basepath)

    restrict_idx = (
                    (cell_metrics.putativeCellType == "Pyramidal Cell") &
                    ((cell_metrics.brainRegion=="CA1") |
                    (cell_metrics.brainRegion=="rCA1") |
                    (cell_metrics.brainRegion=="lCA1")) &
                    (cell_metrics.bad_unit == False) 
                    )

    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]
    
    if cell_metrics.shape[0] == 0:
        return

    # get ripple epochs
    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])
    try:
        st_unit = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx], fs=fs_dat)
    except:
        st_unit = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx][0], fs=fs_dat)

    # behavioral epochs
    epoch_df = loading.load_epoch(basepath)
    
    if not ('behavioralParadigm' in epoch_df.columns):
        epoch_df['behavioralParadigm'] = 'unknown'

    # some epochs will have repeating back to back sleep sessions that are actually the same session
    epoch_df = compress_repeated_epochs.compress_repeated_epochs(epoch_df)

    behavioral_epochs = nel.EpochArray([np.array([epoch_df.startTime,
                                                    epoch_df.stopTime]).T])

    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict['NREMstate'])
    wake_epochs = nel.EpochArray(state_dict['WAKEstate'])

    df_save = main_analysis(
                            basepath,
                            cell_metrics,
                            st_unit,
                            ripple_epochs,
                            behavioral_epochs,
                            epoch_df,
                            nrem_epochs,
                            wake_epochs
                            )
    df_save.to_csv(save_file)

def participation_run(df,save_path,parallel=True):
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