import pandas as pd
import numpy as np 
import os
import functions,loading
import nelpy as nel 
import multiprocessing
from joblib import Parallel, delayed
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols

def load_basic_data(basepath):
    nChannels, fs, fs_dat, shank_to_channel = functions.loadXML(basepath)
    ripples = loading.load_ripples_events(basepath)
    cell_metrics,data = loading.load_cell_metrics(basepath)
    return cell_metrics,data,ripples,fs_dat

def get_participation(st,ripple_epochs):
    # get participation prob.
    # make matrix n rows (units) by n cols (ripple epochs)
    unit_mat = np.zeros((st.n_units,ripple_epochs.n_intervals))
    for i,event in enumerate(st):
        unit_mat[:,i] = (event.n_events>0)*1
    return unit_mat

def get_ripple_fr(st,ripple_epochs):
    # get participation prob.
    # make matrix n rows (units) by n cols (ripple epochs)
    unit_mat = np.zeros((st.n_units,ripple_epochs.n_intervals))
    for i,event in enumerate(st):
        unit_mat[:,i] = event.n_events/ripple_epochs[i].length
    return unit_mat

def bin_unit_mat(unit_mat,behavioral_epochs,bin_width = 120):
    
    bins = np.arange(behavioral_epochs.start,behavioral_epochs.stop,bin_width)
    bin_centers = bins[0:-1] + bin_width/2
    particip = []
    n_ripples = []
    for i in range(bins.shape[0]-1):
        idx = (unit_mat.abscissa_vals >= bins[i]) & (unit_mat.abscissa_vals <= bins[i+1])
        n_ripples.append(sum(idx))
        particip.append(np.sum(unit_mat.data[:,idx] == 1,axis=1)  / unit_mat.data[:,idx].shape[1])

    particip = np.vstack(particip)
    particip[np.isnan(particip)] = 0

    # make sure each bin has at least 1 ripples
    keep_idx = np.array(n_ripples)>0
    particip = particip[keep_idx,:]
    bin_centers = bin_centers[keep_idx]
    
    return nel.AnalogSignalArray(data=particip.T,timestamps=bin_centers)

def anova_table(aov):
    """
    The function was created specifically for the one-way ANOVA table results returned for Type II sum of squares
    """
    aov['mean_sq'] = aov[:]['sum_sq']/aov[:]['df']
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*aov['mean_sq'][-1]))/(sum(aov['sum_sq'])+aov['mean_sq'][-1])

    cols = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[cols]
    return aov

def find_drift(unit_mat_binned,epoch_df):
    """
    Use linear model to detect main effect of epoch
    A significant model will indicate that ripple participation probability is 
    different between epochs
    """
    pval = []
    eta_sq = []
    omega_sq = []
    first_ep_mean_particip = []
    for i in range(unit_mat_binned.data.shape[0]):

        df = pd.DataFrame()
        df['ts'] = unit_mat_binned.abscissa_vals
        df['constant'] = np.ones_like(unit_mat_binned.abscissa_vals)[0]
        df['y'] = unit_mat_binned.data[i]

        for ep in epoch_df.itertuples():
            idx = (df.ts >= ep.startTime) & (df.ts <= ep.stopTime)
            df.loc[idx,'ep'] = ep.name

        m01 = ols('y ~ ep', data=df).fit()
        anovaResults = anova_lm(m01, typ=2, robust="hc3")

        anovaResults = anova_table(anovaResults)

        first_ep_mean_particip.append(df[df.ep == df.ep.unique()[0]].y.mean())

        pval.append(anovaResults['PR(>F)']['ep'])
        eta_sq.append(anovaResults['eta_sq']['ep'])
        omega_sq.append(anovaResults['omega_sq']['ep'])

    df = pd.DataFrame()
    df['pval'] = pval
    df['eta_sq'] = eta_sq
    df['omega_sq'] = omega_sq
    df['first_ep_mean_particip'] = first_ep_mean_particip

    return df

def session_loop(basepath,save_path):

    save_file = os.path.join(save_path,basepath.replace(os.sep, "_").replace(":", "_")  + '.csv')
    if os.path.exists(save_file):
        return

    cell_metrics,data,ripples,fs_dat = load_basic_data(basepath)

    restrict_idx = (
                        (cell_metrics.putativeCellType == "Pyramidal Cell") &
                        ((cell_metrics.brainRegion=="CA1") |
                        (cell_metrics.brainRegion=="rCA1") |
                        (cell_metrics.brainRegion=="lCA1")) &
                        (cell_metrics.bad_unit == False) 
                        )

    # restrict cell metrics                      
    cell_metrics = cell_metrics[restrict_idx]

    # behavioral epochs
    epoch_df = loading.load_epoch(basepath)
    # some epochs will have repeating back to back sessions that are actually the same session
    epoch_df = functions.compress_repeated_epochs(epoch_df)

    # make sure there are enough epochs
    epoch_types = epoch_df.environment.unique()
    if ~(sum(epoch_types != 'sleep') > 1):
        return

    behavioral_epochs = nel.EpochArray([np.array([epoch_df.startTime,
                                                    epoch_df.stopTime]).T])

    # get ripple epochs
    ripple_epochs = nel.EpochArray([np.array([ripples.start,ripples.stop]).T])
    st_unit = nel.SpikeTrainArray(timestamps=np.array(data['spikes'],dtype=object)[restrict_idx], fs=fs_dat)

    # unit_mat = get_ripple_fr(st_unit[ripple_epochs],ripple_epochs)
    unit_mat = get_participation(st_unit[ripple_epochs],ripple_epochs)

    unit_mat = nel.AnalogSignalArray(data=unit_mat,timestamps=ripple_epochs.starts)
    
    unit_mat_binned = bin_unit_mat(unit_mat,behavioral_epochs)

    df_save = find_drift(unit_mat_binned,epoch_df) 

    df_save['basepath'] = basepath

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