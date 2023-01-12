import copy
import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    well_traversal_classification,
    compress_repeated_epochs,
)
from scipy import stats
import pandas as pd
import numpy as np
import nelpy as nel
import os
from scipy import stats
import warnings
from nelpy.analysis import replay

warnings.filterwarnings("ignore")


def get_w_maze_trajectories(
    position_df, max_distance_from_well=20, min_distance_traveled=50
):

    well_locations = np.array(
        [
            [
                position_df.query("states == 0").projected_x.mean(),
                position_df.query("states == 0").projected_y.max(),
            ],
            [
                position_df.query("states == 2").projected_x.mean(),
                position_df.query("states == 2").projected_y.max(),
            ],
            [
                position_df.query("states == 1").projected_x.mean(),
                position_df.query("states == 1").projected_y.max(),
            ],
        ]
    )

    temp_df = position_df[~np.isnan(position_df.x)]
    segments_df, _ = well_traversal_classification.segment_path(
        temp_df["timestamps"].values,
        temp_df[["projected_x", "projected_y"]].values,
        well_locations,
        max_distance_from_well=max_distance_from_well,
    )

    segments_df = well_traversal_classification.score_inbound_outbound(
        segments_df, min_distance_traveled=min_distance_traveled
    )
    conditions = [
        "from_well == 'Center' & to_well == 'Left'",
        "from_well == 'Left' & to_well == 'Center'",
        "from_well == 'Center' & to_well == 'Right'",
        "from_well == 'Right' & to_well == 'Center'",
    ]
    condition_labels = [
        "center_left",
        "left_center",
        "center_right",
        "right_center",
    ]
    trajectories = {}
    for con, con_label in zip(conditions, condition_labels):
        trajectories[con_label] = nel.EpochArray(
            np.array(
                [segments_df.query(con).start_time, segments_df.query(con).end_time]
            ).T
        )

    return trajectories


def flip_pos_within_epoch(pos, dir_epoch):
    """
    flip_pos_within_epoch: flips x coordinate within epoch
        Made to reverse x coordinate within nelpy array for replay analysis
    Input:
        pos: nelpy analog array with single dim
        dir_epoch: epoch to flip
    Output:
        pos: original pos, but fliped by epoch
    """

    def flip_x(x):
        return (x * -1) - np.nanmin(x * -1)

    # make pos df
    pos_df = pd.DataFrame()
    pos_df["ts"] = pos.abscissa_vals
    pos_df["x"] = pos.data.T
    pos_df["dir"] = False

    # make index within df of epoch
    for ep in dir_epoch:
        pos_df.loc[pos_df["ts"].between(ep.starts[0], ep.stops[0]), "dir"] = True

    # flip x within epoch
    pos_df.loc[pos_df.dir == True, "x"] = flip_x(pos_df[pos_df.dir == True].x)

    # add position back to input pos
    pos._data = np.expand_dims(pos_df.x.values, axis=0)

    return pos


def get_tuning_curves(
    pos, st_all, dir_epoch, speed_thres, ds_50ms, s_binsize, tuning_curve_sigma
):
    # compute and smooth speed
    speed1 = nel.utils.ddt_asa(pos, smooth=True, sigma=0.1, norm=True)

    # find epochs where the animal ran > 4cm/sec
    run_epochs = nel.utils.get_run_epochs(
        speed1, v1=speed_thres, v2=speed_thres
    ).merge()

    # restrict spike trains to those epochs during which the animal was running
    st_run = st_all[dir_epoch][run_epochs]

    # smooth and re-bin:
    bst_run = st_run.bin(ds=ds_50ms)

    x_max = np.ceil(np.nanmax(pos[dir_epoch].data))
    x_min = np.floor(np.nanmin(pos[dir_epoch].data))

    n_bins = int((x_max - x_min) / s_binsize)

    tc = nel.TuningCurve1D(
        bst=bst_run,
        extern=pos[dir_epoch][run_epochs],
        n_extern=n_bins,
        extmin=x_min,
        extmax=x_max,
        sigma=tuning_curve_sigma,
        min_duration=0,
    )
    return tc, st_run, bst_run


def restrict_to_place_cells(
    tc,
    st_run,
    bst_run,
    st_all,
    cell_metrics,
    place_cell_min_spks,
    place_cell_min_rate,
    place_cell_peak_mean_ratio,
):
    # locate pyr cells with >= 100 spikes, peak rate >= 1 Hz, peak/mean ratio >=1.5
    peak_firing_rates = tc.max(axis=1)
    mean_firing_rates = tc.mean(axis=1)
    ratio = peak_firing_rates / mean_firing_rates

    idx = (
        (st_run.n_events >= place_cell_min_spks)
        & (tc.ratemap.max(axis=1) >= place_cell_min_rate)
        & (ratio >= place_cell_peak_mean_ratio)
    )
    unit_ids_to_keep = (np.where(idx)[0] + 1).squeeze().tolist()

    sta_placecells = st_all._unit_subset(unit_ids_to_keep)
    tc = tc._unit_subset(unit_ids_to_keep)
    total_units = sta_placecells.n_active
    bst_run = bst_run.loc[:, unit_ids_to_keep]

    # restrict cell_metrics to place cells
    cell_metrics_ = cell_metrics[idx]

    return sta_placecells, tc, bst_run, cell_metrics_, total_units


def decode_and_score(bst, tc, pos):
    # access decoding accuracy on behavioral time scale
    posteriors, lengths, mode_pth, mean_pth = nel.decoding.decode1D(
        bst, tc, xmin=np.nanmin(pos.data), xmax=np.nanmax(pos.data)
    )
    actual_pos = pos(bst.bin_centers)
    bad_idx = np.isnan(actual_pos) | np.isnan(mode_pth)
    actual_pos = actual_pos[~bad_idx]
    mode_pth = mode_pth[~bad_idx]
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(actual_pos, mode_pth)
    median_error = np.nanmedian(np.abs(actual_pos - mode_pth))

    return rvalue, median_error


def pooled_incoherent_shuffle_bst(bst):
    out = copy.deepcopy(bst)
    data = out._data

    for uu in range(bst.n_units):
        segment = np.atleast_1d(np.squeeze(data[uu, :]))
        segment = np.roll(segment, np.random.randint(len(segment)))
        data[uu, :] = segment
    return out


def decode_and_shuff(bst, tc, pos, n_shuffles=500):
    """ """
    rvalue, median_error = decode_and_score(bst, tc, pos)

    scores = np.zeros(bst.n_epochs)
    if n_shuffles > 0:
        rvalue_time_swap = np.zeros((n_shuffles, 1))
        median_error_time_swap = np.zeros((n_shuffles, 1))

    for shflidx in range(n_shuffles):
        bst_shuff = pooled_incoherent_shuffle_bst(bst)
        rvalue_time_swap[shflidx], median_error_time_swap[shflidx] = decode_and_score(
            bst_shuff, tc, pos
        )

    return rvalue, median_error, rvalue_time_swap, median_error_time_swap


def get_features(bst_placecells, posteriors, bdries, mode_pth, pos, dp=3):
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

    # find direction of movement from position
    x_slope = []
    for p in pos:
        good_idx = ~np.isnan(p.data[0])
        cur_x = p.data[0][good_idx]
        if len(cur_x) == 0:
            x_slope.append(np.nan)
        else:
            b1, _, _, _, _ = stats.linregress(np.arange(len(cur_x)), cur_x)
            x_slope.append(b1)
    # if the majority (>.5) of laps have x coords that increase (positive slopes)
    if np.mean(np.array(x_slope) > 0) > 0.5:
        outbound = True
    else:
        outbound = False

    for idx in range(bst_placecells.n_epochs):

        x = bst_placecells[idx].bin_centers

        y = mode_pth[bdries[idx] : bdries[idx + 1]]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        velocity, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
        y = x * velocity + intercept
        position.append(y)

        # get spatial difference between bins
        dy = np.abs(np.diff(y))
        # get cumulative distance
        traj_dist.append(np.nansum(dy))
        # calculate avg speed of trajectory (dist(cm) / time(sec))
        traj_speed.append(np.nansum(dy) / (np.nanmax(x) - np.nanmin(x)))
        # get mean step size
        traj_step.append(np.nanmean(dy))

        if velocity > 0:
            replay_type.append("forward")
        elif velocity < 0:
            replay_type.append("reverse")
        else:
            replay_type.append("unknown")

    return traj_dist, traj_speed, traj_step, replay_type, position


def handle_replay_canidates(basepath, beh_epochs, min_rip_dur):
    # find canidate replay events
    ripples = loading.load_ripples_events(basepath)
    manipulation_df = loading.load_manipulation(
        basepath, struct_name="optoStim", return_epoch_array=False
    )
    closed_loop_idx = manipulation_df.ev_label == "closed_loop"
    manipulation_df.loc[closed_loop_idx, "start"] = (
        manipulation_df[closed_loop_idx].start.values - 0.1
    )
    manipulation_df.loc[closed_loop_idx, "duration"] = (
        manipulation_df[closed_loop_idx].duration.values + 0.1
    )

    ripple_epochs = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)
    manip_epochs = nel.EpochArray(
        np.array([manipulation_df.start, manipulation_df.stop]).T
    )

    # remove ripples that are in the manipulation epochs
    _, idx_overlap_manip = functions.overlap_intersect(
        manip_epochs,
        ripple_epochs,
    )
    ripples = ripples.drop(idx_overlap_manip).reset_index(drop=True)

    # keep ripples that are inside the behavior epochs
    ripple_epochs = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)
    ripples = ripples[
        np.in1d(ripples.start.values, ripple_epochs[beh_epochs].starts)
    ].reset_index(drop=True)

    # combine the ripples and manipulation_df epochs
    ripples = pd.concat([ripples, manipulation_df], ignore_index=True)
    # sort by start time
    ripples.sort_values(by=["start"], inplace=True)
    ripples = ripples.reset_index(drop=True)

    ripple_epochs = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)

    # remove ripples that are too short
    ripples = ripples[ripples.duration > min_rip_dur].reset_index(drop=True)
    ripple_epochs = nel.EpochArray(np.array([ripples.start, ripples.stop]).T)

    # Merge intervals that are close or overlapping
    ripple_epochs = ripple_epochs.merge()

    _, x_ind, _ = np.intersect1d(
        ripples.start.values, ripple_epochs.starts, return_indices=True
    )
    ripples = ripples.loc[x_ind].reset_index(drop=True)

    return ripples, ripple_epochs


def run(
    basepath,
    max_distance_from_well=20,  # in cm, max distance from well to consider a well traversal
    min_distance_traveled=50,  # in cm, min distance traveled to consider a well traversal
    traj_shuff=1500,  # number of shuffles to determine sig replay
    behav_shuff=250,  # number of shuffles to determine sig decoding
    ds_50ms=0.05,  # bin width to bin st for tuning curve
    s_binsize=3,  # spatial bins in tuning curve
    speed_thres=4,  # running threshold to determine tuning curves
    min_rip_dur=0.08,  # min ripple duration for replay
    place_cell_min_rate=1,  # min peak rate of tuning curve
    place_cell_min_spks=100,  # at least 100 spikes while running above speed_thres
    place_cell_peak_mean_ratio=1.5,  # peak firing rate / mean firing rate
    replay_binsize=0.02,  # bin size to decode replay
    tuning_curve_sigma=3,  # 3 cm sd of smoothing on tuning curve
    expand_canidate_by_mua=False,  # whether to expand candidate units by mua (note: will only take rips with mua)
    restrict_manipulation=False,  # whether to restrict manipulation epochs
    shuffle_parallel=False,  # whether to shuffle in parallel
    ds_beh_decode=0.2,  # bin width to bin st for decoding behavior
    add_opto_stim_induced_ripples=True,  # whether to add opto-stim induced ripples
    extend_ripples_dur=0,  # extend ripples by this amount (in sec)
    putativeCellType="Pyr",  # cell type to use for putative cells
    brainRegion="CA1",  # brain region to use
    restrict_ripples_epoch="wmaze",
):

    # load session epoch data
    epoch_df = loading.load_epoch(basepath)
    epoch_df = epoch_df[epoch_df.environment == restrict_ripples_epoch]

    if epoch_df.shape[0] == 0:
        return None

    # get session bounds to provide support
    session_bounds = nel.EpochArray(
        [epoch_df.startTime.iloc[0], epoch_df.stopTime.iloc[-1]]
    )
    # compress repeated sleep sessions
    epoch_df = compress_repeated_epochs.main(epoch_df)
    # put into nel format
    beh_epochs = nel.EpochArray([np.array([epoch_df.startTime, epoch_df.stopTime]).T])

    position_df = loading.load_animal_behavior(basepath)

    if position_df is None:
        return None

    # remove nan values
    bad_idx = np.isnan(position_df.linearized)
    position_df = position_df[~bad_idx]

    # add position to nelpy array
    pos = nel.AnalogSignalArray(
        data=np.array(position_df.linearized),
        timestamps=position_df.time,
        fs=position_df.sr.iloc[0],
    )
    # restrict to wmaze
    pos = pos[beh_epochs[epoch_df.environment == "wmaze"]]

    # locate each trajectory start and end
    trajectories = get_w_maze_trajectories(
        position_df,
        max_distance_from_well=max_distance_from_well,
        min_distance_traveled=min_distance_traveled,
    )

    # flip the x coordinate so it is always increasing
    for con in trajectories.keys():
        x_slope = []
        for pos_seg in pos[trajectories[con]]:
            # use regression to find slope of x coordinate
            b1, _, _, _, _ = stats.linregress(
                np.arange(len(pos_seg.data[0])), pos_seg.data[0]
            )
            x_slope.append(b1)
        # if the majority (>.5) of laps have x coords that decrease
        if np.mean(np.array(x_slope) < 0) > 0.5:
            pos = flip_pos_within_epoch(pos, trajectories[con])

    st_all, cell_metrics = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegion
    )

    if st_all.n_active < 5:
        return

    ripples, ripple_epochs = handle_replay_canidates(basepath, beh_epochs, min_rip_dur)

    results = {}
    for dir_epoch_label in trajectories.keys():
        results[dir_epoch_label] = {}

    for dir_i, dir_epoch_label in enumerate(trajectories.keys()):

        dir_epoch = trajectories[dir_epoch_label]

        if dir_epoch.isempty:
            continue

        # construct tuning curves
        tc, st_run, bst_run = get_tuning_curves(
            pos, st_all, dir_epoch, speed_thres, ds_50ms, s_binsize, tuning_curve_sigma
        )

        # locate pyr cells with >= 100 spikes, peak rate >= 1 Hz, peak/mean ratio >=1.5
        (
            sta_placecells,
            tc,
            bst_run,
            cell_metrics_,
            total_units,
        ) = restrict_to_place_cells(
            tc,
            st_run,
            bst_run,
            st_all,
            cell_metrics,
            place_cell_min_spks,
            place_cell_min_rate,
            place_cell_peak_mean_ratio,
        )
        if tc.isempty:
            continue

        # access decoding accuracy on behavioral time scale
        bst_run_beh = sta_placecells[dir_epoch].bin(ds=ds_beh_decode)
        decoding_r2, median_error, decoding_r2_shuff, _ = decode_and_shuff(
            bst_run_beh, tc, pos[dir_epoch], n_shuffles=behav_shuff
        )

        # check decoding quality against chance distribution
        _, decoding_r2_pval, _ = functions.get_significant_events(
            [decoding_r2], decoding_r2_shuff
        )

        # get ready to decode replay
        # bin data for replay (20ms default)
        bst_placecells = sta_placecells.bin(ds=replay_binsize)[ripple_epochs]

        # count units per event
        n_active = [bst.n_active for bst in bst_placecells]
        n_active = np.array(n_active)
        # also count the proportion of bins in each event with 0 activity
        inactive_bin_prop = [
            sum(bst.n_active_per_bin == 0) / bst.lengths[0] for bst in bst_placecells
        ]
        inactive_bin_prop = np.array(inactive_bin_prop)
        # restrict bst to instances with >= 5 active units and < 50% inactive bins
        idx = (n_active >= 5) & (inactive_bin_prop < 0.5)
        bst_placecells = bst_placecells[np.where(idx)[0]]
        # restrict to instances with >= 5 active units
        n_active = n_active[idx]
        inactive_bin_prop = inactive_bin_prop[idx]
        current_ripples = ripples[idx]

        # bug patch, adding unit_ids and n_intervals for nel.decoding.decode1D
        bst_placecells.unit_ids = tc.unit_ids
        bst_placecells.n_epochs = bst_placecells.n_intervals

        # decode each ripple event
        posteriors, bdries, mode_pth, mean_pth = nel.decoding.decode1D(
            bst_placecells,
            tc,
            xmin=np.nanmin(pos[dir_epoch].data),
            xmax=np.nanmax(pos[dir_epoch].data),
        )

        # score each event using trajectory_score_bst (sums the posterior probability in a range (w) from the LS line)
        scores, scores_time_swap, scores_col_cycle = replay.trajectory_score_bst(
            bst_placecells, tc, w=3, n_shuffles=traj_shuff, normalize=True
        )

        _, score_pval_time_swap, _ = functions.get_significant_events(
            scores, scores_time_swap
        )
        _, score_pval_col_cycle, _ = functions.get_significant_events(
            scores, scores_col_cycle
        )

        (traj_dist, traj_speed, traj_step, replay_type, position) = get_features(
            bst_placecells, posteriors, bdries, mode_pth, pos[dir_epoch], dp=s_binsize
        )

        slope, intercept, r2values = replay.linregress_bst(bst_placecells, tc)

        # package data into results dictionary
        results[dir_epoch_label]["cell_metrics"] = cell_metrics_

        results[dir_epoch_label]["sta_placecells"] = sta_placecells
        results[dir_epoch_label]["bst_placecells"] = bst_placecells
        results[dir_epoch_label]["bst_run"] = bst_run
        results[dir_epoch_label]["bst_run_beh"] = bst_run_beh
        results[dir_epoch_label]["pos"] = pos[dir_epoch]
        results[dir_epoch_label]["tc"] = tc
        results[dir_epoch_label]["posteriors"] = posteriors
        results[dir_epoch_label]["bdries"] = bdries
        results[dir_epoch_label]["mode_pth"] = mode_pth
        results[dir_epoch_label]["position"] = position

        temp_df = current_ripples.copy()
        # add event by event metrics to df
        temp_df["n_active"] = n_active
        temp_df["inactive_bin_prop"] = inactive_bin_prop
        temp_df["trajectory_score"] = scores
        temp_df["r_squared"] = r2values
        temp_df["slope"] = slope
        temp_df["intercept"] = intercept
        temp_df["score_pval_time_swap"] = score_pval_time_swap
        temp_df["score_pval_col_cycle"] = score_pval_col_cycle
        temp_df["traj_dist"] = traj_dist
        temp_df["traj_speed"] = traj_speed
        temp_df["traj_step"] = traj_step
        temp_df["replay_type"] = replay_type

        results[dir_epoch_label]["df"] = temp_df
        results[dir_epoch_label]["session"] = basepath
        results[dir_epoch_label]["decoding_r2"] = decoding_r2
        results[dir_epoch_label]["decoding_r2_pval"] = decoding_r2_pval
        results[dir_epoch_label]["decoding_median_error"] = median_error
        results[dir_epoch_label]["total_units"] = total_units

    return results


def load_results(save_path):
    """
    Load results from a previous run.
    """

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    df = pd.DataFrame()

    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        # locate basepath
        for key_ in results.keys():
            try:
                basepath = results[key_]["session"]
                break
            except:
                continue

        for key_ in results.keys():
            try:
                results[key_]["sta_placecells"]
            except:
                continue

            # calc and add ripple participation
            st = results[key_]["sta_placecells"]
            bst = results[key_]["bst_placecells"]
            if len(bst.support.starts) == 0:
                results[key_]["df"]["pop_partic"] = 0
            else:
                particip_mat = functions.get_participation(
                    st.data, bst.support.starts, bst.support.stops
                )
                results[key_]["df"]["pop_partic"] = particip_mat.mean(axis=0)

            # add behavioral decoding quality
            results[key_]["df"]["decoding_r2"] = float(results[key_]["decoding_r2"])
            results[key_]["df"]["decoding_r2_pval"] = float(
                results[key_]["decoding_r2_pval"]
            )
            results[key_]["df"]["decoding_median_error"] = float(
                results[key_]["decoding_median_error"]
            )
            results[key_]["df"]["total_units"] = float(results[key_]["total_units"])
            results[key_]["df"]["direction"] = key_

            results[key_]["df"]["basepath"] = basepath
            df = pd.concat([df, results[key_]["df"]], ignore_index=True)

    return df
