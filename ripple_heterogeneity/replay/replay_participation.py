import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import (
    functions,
    loading,
    compress_repeated_epochs
)
import nelpy as nel

def get_replay_epochs(results, direction, alpha=0.05):
    """
    Get replay epochs from results
    Inputs:
        results: pandas dataframe of results
        direction: 'forward' or 'reverse'
        alpha: p value threshold
    Outputs:
        replay_epochs: list of replay epochs
    
    """
    # check replay direction input
    if not np.any((direction == "forward") | (direction == "reverse")):
        print("wrong direction")
        return

    starts = []
    stops = []
    # get sig replay epochs from outbound
    idx = (results["outbound_epochs"]["df"].score_pval_time_swap <= alpha) & (
        results["outbound_epochs"]["df"].replay_type == direction
    )
    starts.append(results["outbound_epochs"]["df"][idx].start)
    stops.append(results["outbound_epochs"]["df"][idx].stop)

    # get sig replay epochs from inbound
    idx = (results["inbound_epochs"]["df"].score_pval_time_swap <= alpha) & (
        results["inbound_epochs"]["df"].replay_type == direction
    )
    starts.append(results["inbound_epochs"]["df"][idx].start)
    stops.append(results["inbound_epochs"]["df"][idx].stop)

    # get sorted replay epochs and sort
    sort_idx = np.argsort(np.hstack(starts))
    starts = np.hstack(starts)[sort_idx]
    stops = np.hstack(stops)[sort_idx]

    return nel.EpochArray(np.array([starts, stops]).T)

def run(session):
    """
    Compile replay participation for a session
    Inputs:
        session: session directory
    Outputs:
        temp_df: pandas dataframe of replay and ripple participation
    
    """
    # load results
    with open(session, "rb") as f:
        results = pickle.load(f)

    # skip if no results
    if results is None:
        return pd.DataFrame()

    # retrive basepath
    try:
        basepath = results["outbound_epochs"]['session']
    except:
        basepath = results["inbound_epochs"]['session']

    # get forward and reverse replay epochs
    forward_replay = get_replay_epochs(results, "forward")
    reverse_replay = get_replay_epochs(results, "reverse")

    # get session epochs
    epoch_df = loading.load_epoch(basepath)
    # compress back to back sleep epochs
    epoch_df = compress_repeated_epochs.main(epoch_df)
    # add session epochs into nelpy epoch array
    behavior_epochs = nel.EpochArray(
        [np.array([epoch_df.startTime, epoch_df.stopTime]).T]
    )
    # get ripple epochs
    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])

    # state_dict = loading.load_SleepState_states(basepath)
    # nrem_epochs = nel.EpochArray(state_dict['NREMstate'])
    # wake_epochs = nel.EpochArray(state_dict['WAKEstate'])

    # locate active units that were used in anlysis
    # because out/inbound templates were used seperately, we need to include both
    # always load cell metrics from source to get most up to date data
    st, cell_metrics = loading.load_spikes(basepath)

    # get active units from cell metrics
    uid = pd.unique(
        np.hstack(
            [
                results["outbound_epochs"]["cell_metrics"].UID,
                results["inbound_epochs"]["cell_metrics"].UID,
            ]
        )
    )

    # remove uids with bad waveforms as we can not estimate deep/sup
    if "tags_bad_waveform" in cell_metrics.columns:
        a = set(cell_metrics[cell_metrics.tags_bad_waveform].UID.values)
        b = set(uid)
        c = b.difference(a)
        uid = np.sort(np.array(list(c)))
    
    _, x_ind, _ = np.intersect1d(cell_metrics.UID, uid, return_indices=True)
    unit_ids_to_keep = (x_ind + 1).squeeze().tolist()
    sta_placecells = st._unit_subset(unit_ids_to_keep)
    cell_metrics = cell_metrics.iloc[x_ind]

    # epoch array of all replay epochs
    all_replay = forward_replay | reverse_replay

    # epoch array of all ripple epochs that are not replay
    ripple_outside_replay = ripple_epochs[~all_replay]

    epoch = []
    epoch_i = []
    replay_par = []
    ripple_par = []
    UID = []
    n_ripples = []
    n_replays = []
    # iterate through all behavioral epochs and get ripple and replay participation
    for beh_ep_i, beh_ep in enumerate(behavior_epochs):

        replay_par.append(functions.get_participation(
            sta_placecells[beh_ep].data,
            all_replay[beh_ep].starts,
            all_replay[beh_ep].stops,
        ).mean(axis=1))

        ripple_par.append(functions.get_participation(
            sta_placecells[beh_ep].data,
            ripple_outside_replay[beh_ep].starts,
            ripple_outside_replay[beh_ep].stops,
        ).mean(axis=1))

        epoch.append([epoch_df.environment.values[beh_ep_i]]*sta_placecells.data.shape[0])
        epoch_i.append(np.tile(beh_ep_i,sta_placecells.data.shape[0]))

        UID.append(cell_metrics.UID.values)

        n_replays.append(np.tile(all_replay[beh_ep].n_intervals,sta_placecells.data.shape[0]))
        n_ripples.append(np.tile(ripple_outside_replay[beh_ep].n_intervals,sta_placecells.data.shape[0]))
        
    # stack all data
    temp_df = pd.DataFrame()
    temp_df["replay_par"] = np.hstack(replay_par)
    temp_df["ripple_par"] = np.hstack(ripple_par)
    temp_df["epoch"] = np.hstack(epoch)
    temp_df["epoch_i"] = np.hstack(epoch_i)
    temp_df["UID"] = np.hstack(UID)
    temp_df["n_replays"] = np.hstack(n_replays)
    temp_df["n_ripples"] = np.hstack(n_ripples)
    temp_df["basepath"] = basepath

    return temp_df

