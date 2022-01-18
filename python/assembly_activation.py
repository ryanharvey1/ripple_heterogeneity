import multiprocessing
from joblib import Parallel, delayed
import os
import pandas as pd
import numpy as np
import nelpy as nel
import pickle
import assembly_run
import loading
import sys
sys.path.append(r'D:\github\neurocode\reactivation\assemblies')
import assembly

def session_loop_activation(basepath,save_path,save_path_assembly):

    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.pkl')
    if os.path.exists(save_file):
        return

    assembly_file = os.path.join(save_path_assembly,basepath.replace(os.sep, "_").replace(":", "_")  + '.pkl')
    with open(assembly_file, 'rb') as f:
            results = pickle.load(f)

    cell_metrics,data,ripples,fs_dat = loading.load_basic_data(basepath)

    restrict_idx = ((cell_metrics.putativeCellType == "Pyramidal Cell") &
                        ((cell_metrics.brainRegion=="CA1") |
                        (cell_metrics.brainRegion=="rCA1") |
                        (cell_metrics.brainRegion=="lCA1")) &
                        (cell_metrics.bad_unit == False))

    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]

    st_unit = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx], fs=fs_dat)
    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])

    epoch_df = loading.load_epoch(basepath)
    beh_epochs = nel.EpochArray([np.array([epoch_df.startTime,epoch_df.stopTime]).T])

    st_unit_rip = st_unit[ripple_epochs]

    assembly_act = []
    for ep in beh_epochs:
        z_mat,ts = assembly_run.get_z_t(st_unit_rip[ep],ds=0.002)

        assembly_act.append( 
                            nel.AnalogSignalArray(
                                    data=assembly.computeAssemblyActivity(results['patterns_inside_ripples'], z_mat),
                                    timestamps=ts,
                                    fs=500
                                    )
                            )

    results['assembly_act_inside_ripples'] = assembly_act
    results['basepath'] = basepath
    results['cell_metrics'] = cell_metrics
    results['st_unit'] = st_unit
    results['epoch_df'] = epoch_df

    # save file
    with open(save_file, 'wb') as f:
        pickle.dump(results, f)

def assembly_run_activation(df,save_path,save_path_assembly,parallel=True):
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()         
        processed_list = Parallel(n_jobs=num_cores)(delayed(session_loop_activation)(basepath,save_path) for basepath in basepaths)
    else:    
        for basepath in basepaths:
            print(basepath)
            session_loop_activation(basepath,save_path,save_path_assembly)