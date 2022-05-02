import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, compress_repeated_epochs
import nelpy as nel


def run(basepath, replay_df=None, replay_save_path=None, alpha=0.05):
    """
    Compile replay participation for a session
    Inputs:
        session: session directory
    Outputs:
        temp_df: pandas dataframe of replay and ripple participation

    """
    # select current session from replay_df
    replay_df = replay_df[replay_df.basepath == basepath]

    starts = replay_df[
        (replay_df.score_pval_time_swap <= alpha) & (replay_df.replay_type == "forward")
    ].start
    stops = replay_df[
        (replay_df.score_pval_time_swap <= alpha) & (replay_df.replay_type == "forward")
    ].stop
    forward_replay = nel.EpochArray(np.array([starts, stops]).T)

    starts = replay_df[
        (replay_df.score_pval_time_swap <= alpha) & (replay_df.replay_type == "reverse")
    ].start
    stops = replay_df[
        (replay_df.score_pval_time_swap <= alpha) & (replay_df.replay_type == "reverse")
    ].stop
    reverse_replay = nel.EpochArray(np.array([starts, stops]).T)

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

    session = os.path.join(
        replay_save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )

    with open(session, "rb") as f:
        results = pickle.load(f)

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
    deepSuperficialDistance = []
    replay_fr = []
    ripple_fr = []
    avg_fr = []
    # iterate through all behavioral epochs and get ripple and replay participation
    for beh_ep_i, beh_ep in enumerate(behavior_epochs):

        current_st = sta_placecells[beh_ep]

        # get avg firing rate over epoch
        avg_fr.append(current_st.n_events / beh_ep.length)
        
        # get ripple firing rate
        # check if any spikes left in epoch
        if sta_placecells[beh_ep][all_replay].data is None:
            replay_fr.append(current_st.n_events * np.nan)
        else:
            replay_fr.append(
                sta_placecells[beh_ep][all_replay].n_events / beh_ep.length
            )

        # get ripple firing rate
        if sta_placecells[beh_ep][ripple_outside_replay].data is None:
            ripple_fr.append(current_st.n_events * np.nan)
        else:
            ripple_fr.append(
                sta_placecells[beh_ep][ripple_outside_replay].n_events / beh_ep.length
            )

        # get replay participation
        replay_par.append(
            functions.get_participation(
                current_st.data,
                all_replay[beh_ep].starts,
                all_replay[beh_ep].stops,
            ).mean(axis=1)
        )

        # get ripple participation
        ripple_par.append(
            functions.get_participation(
                current_st.data,
                ripple_outside_replay[beh_ep].starts,
                ripple_outside_replay[beh_ep].stops,
            ).mean(axis=1)
        )

        # get epoch info
        epoch.append(
            [epoch_df.environment.values[beh_ep_i]] * sta_placecells.data.shape[0]
        )
        epoch_i.append(np.tile(beh_ep_i, sta_placecells.data.shape[0]))
        # get UID
        UID.append(cell_metrics.UID.values)
        # get deep superficial distance
        deepSuperficialDistance.append(cell_metrics.deepSuperficialDistance.values)
        # get number of ripples and replays
        n_replays.append(
            np.tile(all_replay[beh_ep].n_intervals, sta_placecells.data.shape[0])
        )
        n_ripples.append(
            np.tile(
                ripple_outside_replay[beh_ep].n_intervals, sta_placecells.data.shape[0]
            )
        )

    # stack all data
    temp_df = pd.DataFrame()
    temp_df["avg_fr"] = np.hstack(avg_fr)
    temp_df["replay_fr"] = np.hstack(replay_fr)
    temp_df["ripple_fr"] = np.hstack(ripple_fr)
    temp_df["replay_par"] = np.hstack(replay_par)
    temp_df["ripple_par"] = np.hstack(ripple_par)
    temp_df["epoch"] = np.hstack(epoch)
    temp_df["epoch_i"] = np.hstack(epoch_i)
    temp_df["UID"] = np.hstack(UID)
    temp_df["deepSuperficialDistance"] = np.hstack(deepSuperficialDistance)
    temp_df["n_replays"] = np.hstack(n_replays)
    temp_df["n_ripples"] = np.hstack(n_ripples)
    temp_df["basepath"] = basepath

    return temp_df
