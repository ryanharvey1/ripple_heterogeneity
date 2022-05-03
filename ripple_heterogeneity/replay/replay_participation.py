import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, compress_repeated_epochs
import nelpy as nel


def run(basepath, replay_df=None, replay_save_path=None, alpha=0.05, partic_pad=0.05):
    """
    Compile replay participation for a session
    Inputs:
        session: session directory
    Outputs:
        temp_df: pandas dataframe of replay and ripple participation

    """
    # select current session from replay_df
    replay_df = replay_df[replay_df.basepath == basepath]

    sig_replay_idx = replay_df.score_pval_time_swap <= alpha

    starts = replay_df[sig_replay_idx & (replay_df.replay_type == "forward")].start
    stops = replay_df[sig_replay_idx & (replay_df.replay_type == "forward")].stop
    forward_replay = nel.EpochArray(np.array([starts, stops]).T)

    starts = replay_df[sig_replay_idx & (replay_df.replay_type == "reverse")].start
    stops = replay_df[sig_replay_idx & (replay_df.replay_type == "reverse")].stop
    reverse_replay = nel.EpochArray(np.array([starts, stops]).T)

    non_sig_replay_idx = replay_df.score_pval_time_swap > alpha
    starts = replay_df[non_sig_replay_idx].start
    stops = replay_df[non_sig_replay_idx].stop
    canidate_non_replay = nel.EpochArray(np.array([starts, stops]).T)

    # get session epochs
    epoch_df = loading.load_epoch(basepath)
    # compress back to back sleep epochs
    epoch_df = compress_repeated_epochs.main(epoch_df)
    # get index of linear track relative to other (pre is negative, post is positive)
    idx = epoch_df.environment == "linear"
    if sum(idx) > 1:
        duration = epoch_df.stopTime - epoch_df.startTime
        idx = idx & (duration == max(duration[idx]))
    beh_idx_rel_to_linear = np.arange(epoch_df.shape[0]) - np.where(idx)[0][0]

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

    ripple_outside_replay = ripple_outside_replay.expand(partic_pad)
    canidate_non_replay = canidate_non_replay.expand(partic_pad)
    if not all_replay.isempty:
        all_replay = all_replay.expand(partic_pad)
    if not forward_replay.isempty:
        forward_replay = forward_replay.expand(partic_pad)
    if not reverse_replay.isempty:
        reverse_replay = reverse_replay.expand(partic_pad)

    epoch = []
    epoch_i = []
    epoch_rel_i = []
    replay_par = []
    ripple_par = []
    UID = []
    n_ripples = []
    n_replays = []
    deepSuperficialDistance = []
    replay_fr = []
    ripple_fr = []
    avg_fr = []
    non_replay_par = []
    forward_replay_par = []
    reverse_replay_par = []
    non_replay_fr = []

    # iterate through all behavioral epochs and get ripple and replay participation
    for beh_ep_i, beh_ep in enumerate(behavior_epochs):

        current_st = sta_placecells[beh_ep]

        # get avg firing rate over epoch
        avg_fr.append(current_st.n_events / beh_ep.length)

        # get ripple firing rate
        # check if any spikes left in epoch
        if current_st[all_replay].isempty:
            replay_fr.append(current_st.n_events * np.nan)
        else:
            replay_fr.append(
                current_st[all_replay].n_events / all_replay[beh_ep].length
            )

        if current_st[all_replay].isempty:
            non_replay_fr.append(current_st.n_events * np.nan)
        else:
            non_replay_fr.append(
                current_st[canidate_non_replay].n_events
                / canidate_non_replay[beh_ep].length
            )

        # get ripple firing rate
        if current_st[ripple_outside_replay].isempty:
            ripple_fr.append(current_st.n_events * np.nan)
        else:
            ripple_fr.append(
                current_st[ripple_outside_replay].n_events
                / ripple_outside_replay[beh_ep].length
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

        # get canidate participation
        non_replay_par.append(
            functions.get_participation(
                current_st.data,
                canidate_non_replay[beh_ep].starts,
                canidate_non_replay[beh_ep].stops,
            ).mean(axis=1)
        )

        # get forward replay participation
        if not forward_replay.isempty:
            forward_replay_par.append(
                functions.get_participation(
                    current_st.data,
                    forward_replay[beh_ep].starts,
                    forward_replay[beh_ep].stops,
                ).mean(axis=1)
            )
        else:
            forward_replay_par.append(current_st.n_events * np.nan)

        # get reverse replay participation
        if not forward_replay.isempty:
            reverse_replay_par.append(
                functions.get_participation(
                    current_st.data,
                    reverse_replay[beh_ep].starts,
                    reverse_replay[beh_ep].stops,
                ).mean(axis=1)
            )
        else:
            reverse_replay_par.append(current_st.n_events * np.nan)

        # get epoch info
        epoch.append(
            [epoch_df.environment.values[beh_ep_i]] * sta_placecells.data.shape[0]
        )
        epoch_i.append(np.tile(beh_ep_i, sta_placecells.data.shape[0]))
        epoch_rel_i.append(
            np.tile(beh_idx_rel_to_linear[beh_ep_i], sta_placecells.data.shape[0])
        )

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
    temp_df["non_replay_fr"] = np.hstack(non_replay_fr)
    temp_df["replay_par"] = np.hstack(replay_par)
    temp_df["ripple_par"] = np.hstack(ripple_par)
    temp_df["non_replay_par"] = np.hstack(non_replay_par)
    temp_df["forward_replay_par"] = np.hstack(forward_replay_par)
    temp_df["reverse_replay_par"] = np.hstack(reverse_replay_par)
    temp_df["epoch"] = np.hstack(epoch)
    temp_df["epoch_i"] = np.hstack(epoch_i)
    temp_df["epoch_rel_i"] = np.hstack(epoch_rel_i)
    temp_df["UID"] = np.hstack(UID)
    temp_df["deepSuperficialDistance"] = np.hstack(deepSuperficialDistance)
    temp_df["n_replays"] = np.hstack(n_replays)
    temp_df["n_ripples"] = np.hstack(n_ripples)
    temp_df["basepath"] = basepath

    return temp_df
