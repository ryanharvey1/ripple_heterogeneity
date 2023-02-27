from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.readout import assembly_multi_region
import pandas as pd
import numpy as np
import os
import nelpy as nel
import glob
import pickle
import pynapple as nap
import logging

logging.getLogger().setLevel(logging.ERROR)


def locate_task_epoch(epoch_df, env, position_df):
    """
    locate_task_epoch: here, I'm trying to find the longest behavioral
        epoch that has linearized coordinates. The issue is that some older
        data has epochs labeled as tmaze, but no tracking data.
    """
    epoch_df = epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    # epoch_df.sort_values("duration", ascending=False, inplace=True)
    # return int(epoch_df[epoch_df.environment.str.contains(env)].index[0])

    task_idx_longest_epoch = (
        epoch_df.query("environment.str.contains(@env)")
        .sort_values(by="duration", ascending=False)
        .index
    )
    valid_epoch_names = position_df[~position_df.linearized.isnull()].epochs.unique()
    task_idx_valid_position = np.where(np.isin(epoch_df.name, valid_epoch_names))[0]

    # epoch intersection between longest and valid position
    possible_epochs = list(set(task_idx_longest_epoch) & set(task_idx_valid_position))

    # if longest is in possible, return
    if task_idx_longest_epoch[0] in possible_epochs:
        return int(task_idx_longest_epoch[0])
    else:
        # else return the first possible epoch
        return int(possible_epochs[0])


def dissociate_laps_by_states(states, dir_epoch, states_of_interest=[1, 2]):
    # unique_states = np.unique(states.data[~np.isnan(states.data)])
    lap_id = []
    for ep in dir_epoch:
        state_count = []
        for us in states_of_interest:
            state_count.append(np.nansum(states[ep].data == us))
        lap_id.append(states_of_interest[np.argmax(state_count)])
    return np.array(lap_id).astype(int)


def get_pos(basepath, epochs, epoch_df, task_idx):
    position_df = loading.load_animal_behavior(basepath)
    position_df_no_nan = position_df.query("not x.isnull() & not y.isnull()")

    if position_df_no_nan.shape[0] == 0:
        return None, None, None, None

    if "linearized" not in position_df_no_nan.columns:
        return None, None, None, None

    if "states" not in position_df_no_nan.columns:
        return None, None, None, None

    pos = nel.PositionArray(
        data=position_df_no_nan["linearized"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    # make min pos 2
    # pos._data = (pos.data - np.nanmin(pos.data)) + 2

    pos = pos[epochs[task_idx]]

    states = nel.AnalogSignalArray(
        data=position_df_no_nan["states"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    states = states[epochs[task_idx]]

    # get outbound and inbound epochs
    outbound_epochs, inbound_epochs = functions.get_linear_track_lap_epochs(
        pos.abscissa_vals, pos.data[0], newLapThreshold=20
    )
    outbound_epochs = functions.find_good_lap_epochs(pos, outbound_epochs)
    inbound_epochs = functions.find_good_lap_epochs(pos, inbound_epochs)

    if not inbound_epochs.isempty:
        raise TypeError("inbound_epochs should be empty for tmaze")

    # locate laps with the majority in state 1 or 2
    lap_id = dissociate_laps_by_states(
        states, outbound_epochs, states_of_interest=[1, 2]
    )

    right_epochs = nel.EpochArray(data=outbound_epochs.data[lap_id == 1, :])
    left_epochs = nel.EpochArray(data=outbound_epochs.data[lap_id == 2, :])

    position_df_no_nan = position_df_no_nan[
        position_df_no_nan["time"].between(
            epoch_df.iloc[task_idx].startTime, epoch_df.iloc[task_idx].stopTime
        )
    ]
    return pos, right_epochs, left_epochs, states, position_df_no_nan


def find_closest_position_index(position_df_no_nan, xy):
    idx = np.argmin(np.abs(position_df_no_nan.x - xy[0] + position_df_no_nan.y - xy[1]))
    return idx


def get_cross_region_assemblies(
    basepath,
    regions=None,  # brain regions to load
    deepSuperficial=None,
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,
    epoch=None,
):
    m1 = assembly_reactivation.AssemblyReact(
        basepath,
        brainRegion=regions[0] + "|" + regions[1],
        putativeCellType=putativeCellType,
        weight_dt=weight_dt,
        z_mat_dt=z_mat_dt,
    )
    st, cell_metrics = loading.load_spikes(
        basepath,
        brainRegion=regions[0] + "|" + regions[1],
        putativeCellType=putativeCellType,
    )
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)
    m1.st = st.iloc[
        :,
        cell_metrics.deepSuperficial.str.contains(deepSuperficial + "|unknown").values,
    ]
    m1.cell_metrics = cell_metrics[
        cell_metrics.deepSuperficial.str.contains(deepSuperficial + "|unknown").values
    ]
    # check for any cortex cells
    if (
        not m1.cell_metrics.brainRegion.str.contains(regions[0]).any()
        & m1.cell_metrics.brainRegion.str.contains(regions[1]).any()
    ):
        return None, None

    m1.get_weights(epoch=epoch)
    assembly_act = m1.get_assembly_act(epoch=epoch)
    return m1, assembly_act


def locate_t_maze_key_locations(current_pos):
    start_pos = (
        current_pos.dropna(subset=["linearized", "x", "y"])
        .query("states==0")
        .sort_values(by="linearized")[["x", "y"]][:10]
        .mean()
        .values
    )
    decision_pos = (
        current_pos.dropna(subset=["linearized", "x", "y"])
        .query("states==0")
        .sort_values(by="linearized", ascending=False)[["x", "y"]][:10]
        .mean()
        .values
    )
    reward_right_pos = start_pos.copy()
    reward_left_pos = start_pos.copy()
    return start_pos, decision_pos, reward_left_pos, reward_right_pos


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    cross_regions=(("CA1", "PFC"), ("CA1", "EC1|EC2|EC3|EC4|EC5|MEC")),
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,  # dt in seconds for binning st to get activation strength
    verbose=False,  # print out progress
    env="Mwheel|Tmaze|tmaze",  # enviroment you want to look at (current should only be tmaze)
    s_binsize=3,  # spatial bin size
    smooth_sigma=3,  # smoothing sigma in cm
    smooth_window=10,  # smoothing window in cm
    speed_thres=4,  # speed threshold for ratemap in cm/sec
):

    # locate and load linearization file to get key maze locations
    linearization_file = os.path.join(basepath, "linearization_nodes_edges.pkl")
    # if this file doesn't exist, skip
    if not os.path.exists(linearization_file):
        return None
    with open(linearization_file, "rb") as f:
        nodes_and_edges = pickle.load(f)

    epoch_df = loading.load_epoch(basepath)
    epochs = nel.EpochArray(np.array([epoch_df.startTime, epoch_df.stopTime]).T)

    # check if any env
    if not epoch_df.environment.str.contains(env).any():
        return None

    # find longest xx session
    position_df = loading.load_animal_behavior(basepath)

    task_idx = locate_task_epoch(epoch_df, env, position_df)

    # load theta epochs for later restriction when detecting assemblies
    state_dict = loading.load_SleepState_states(basepath)
    theta_epochs = nel.EpochArray(state_dict["THETA"])

    # locate key points (TODO: fix the hard coded values)
    # FujisawaS/AYA10 and Mwheel data will have these reward zones
    if (
        "FujisawaS" in basepath
        or "Mwheel" in epoch_df.iloc[task_idx]
        or "AYA10" in basepath
    ):
        start_pos = nodes_and_edges["node_positions"][0]
        decision_pos = nodes_and_edges["node_positions"][1]
        reward_left_pos = nodes_and_edges["node_positions"][3]
        reward_right_pos = nodes_and_edges["node_positions"][5]
    else:
        # locate key points from states
        current_pos = position_df.query(
            "time >= @epoch_df.iloc[@task_idx].startTime & time <= @epoch_df.iloc[@task_idx].stopTime"
        )
        (
            start_pos,
            decision_pos,
            reward_left_pos,
            reward_right_pos,
        ) = locate_t_maze_key_locations(current_pos)

    # get cross region assemblies
    assem_labels = []
    assembly_act = []
    abscissa_vals = []
    m1 = {}
    for cross_region in cross_regions:
        for deepSuperficial in ["Deep", "Superficial"]:
            m1_, assembly_act_ = get_cross_region_assemblies(
                basepath,
                regions=cross_region,
                deepSuperficial=deepSuperficial,
                putativeCellType=putativeCellType,
                weight_dt=weight_dt,
                z_mat_dt=z_mat_dt,
                epoch=epochs[task_idx][theta_epochs],
            )
            if assembly_act_ is None:
                continue

            # make sure assembly members represent cross region label
            keep_assembly = []
            _, _, keep_assembly_, is_member = find_sig_assembly.main(m1_.patterns)
            for assembly_i in range(m1_.n_assemblies()):
                member_idx = is_member[assembly_i, :]
                cortex_check = (
                    m1_.cell_metrics[member_idx]
                    .brainRegion.str.contains(cross_region[1])
                    .any()
                )
                ca1_check = (
                    m1_.cell_metrics[member_idx]
                    .brainRegion.str.contains(cross_region[0])
                    .any()
                )
                ca1_layer_check = (
                    m1_.cell_metrics[member_idx].deepSuperficial == deepSuperficial
                ).any()
                keep_assembly.append(cortex_check & ca1_check & ca1_layer_check)

            assembly_act.append(assembly_act_.data[keep_assembly, :])
            abscissa_vals.append(assembly_act_.abscissa_vals)
            assem_labels.append(
                [deepSuperficial + "_" + cross_region[1]] * sum(keep_assembly)
            )
            m1[deepSuperficial + "_" + cross_region[1]] = m1_

    if len(assembly_act) == 0:
        return None

    assembly_act = nel.AnalogSignalArray(
        data=np.vstack(assembly_act),
        timestamps=abscissa_vals[0],
    )
    
    if assembly_act.isempty:
        return None

    assem_labels = np.hstack(assem_labels)

    pos, right_epochs, left_epochs, states, position_df_no_nan = get_pos(
        basepath, epochs, epoch_df, task_idx
    )
    if pos is None:
        return None

    # locate key locations in linear coords
    # restrict to particular state when locating closed position
    # TODO: fix the hard coded values

    idx = find_closest_position_index(position_df_no_nan.query("states==0"), start_pos)
    x_start = np.interp(
        position_df_no_nan.query("states==0").iloc[idx].time,
        pos.abscissa_vals,
        pos.data[0],
    )

    idx = find_closest_position_index(
        position_df_no_nan.query("states==0"), decision_pos
    )
    x_decision = np.interp(
        position_df_no_nan.query("states==0").iloc[idx].time,
        pos.abscissa_vals,
        pos.data[0],
    )

    idx = find_closest_position_index(
        position_df_no_nan.query("states==1"), reward_left_pos
    )
    x_reward_left = np.interp(
        position_df_no_nan.query("states==1").iloc[idx].time,
        pos.abscissa_vals,
        pos.data[0],
    )

    idx = find_closest_position_index(
        position_df_no_nan.query("states==2"), reward_right_pos
    )
    x_reward_right = np.interp(
        position_df_no_nan.query("states==2").iloc[idx].time,
        pos.abscissa_vals,
        pos.data[0],
    )

    ext_xmin, ext_xmax = (
        np.floor(pos.data[0].min() / 10) * 10,
        np.ceil(pos.data[0].max() / 10) * 10,
    )
    n_extern = int((ext_xmax - ext_xmin) / s_binsize)

    label_df = pd.DataFrame()
    tc = pd.DataFrame()

    direction_label = ["right_epochs", "left_epochs"]

    for dir_epoch_i, dir_epoch in enumerate([right_epochs, left_epochs]):

        # compute and smooth speed
        speed1 = nel.utils.ddt_asa(pos[dir_epoch], smooth=True, sigma=0.1, norm=True)

        # find epochs where the animal ran > 4cm/sec
        run_epochs = nel.utils.get_run_epochs(speed1, v1=speed_thres, v2=speed_thres)

        # add assembly activations to tsdframe object
        tsdframe = nap.TsdFrame(
            t=assembly_act[dir_epoch][run_epochs].abscissa_vals,
            d=assembly_act[dir_epoch][run_epochs].data.T,
        )

        # add position to tsd object
        feature = nap.Tsd(
            t=pos[dir_epoch][run_epochs].abscissa_vals,
            d=pos[dir_epoch][run_epochs].data[0],
        )

        # compute tuning curves
        tc_ = nap.compute_1d_tuning_curves_continous(
            tsdframe, feature, nb_bins=n_extern, minmax=[ext_xmin, ext_xmax]
        )

        # smooth
        tc_ = tc_.rolling(
            window=smooth_window, win_type="gaussian", center=True, min_periods=1
        ).mean(std=smooth_sigma)

        # store tuning curves
        tc = pd.concat([tc, tc_], axis=1, ignore_index=True)

        # store associated metadata
        label_df_ = pd.DataFrame()
        label_df_["assembly_n"] = np.arange(assembly_act.n_signals).astype(int)
        label_df_["direction_label"] = direction_label[dir_epoch_i]
        label_df_["cross_region_label"] = assem_labels
        # label_df_ = label_df_.merge(counts_df.query("is_member_sig"))

        label_df = pd.concat([label_df, label_df_], ignore_index=True)

    label_df["x_start"] = x_start
    label_df["x_decision"] = x_decision
    label_df["x_reward_left"] = x_reward_left
    label_df["x_reward_right"] = x_reward_right
    label_df["basepath"] = basepath

    results = {
        "tc": tc,
        "label_df": label_df,
        "assembly_act_task": assembly_act,
        "react": m1,
        "pos": pos,
        "right_epochs": right_epochs,
        "left_epochs": left_epochs,
        "states": states,
        "task_epoch": epochs[task_idx],
    }

    return results


def load_results(save_path, verbose=False):
    """
    load_results: load results from a pickle file
    """
    print("Loading results...")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    tc = []
    label_df = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        tc.append(results["tc"])
        label_df = pd.concat([label_df, results["label_df"]], ignore_index=True)

    return tc, label_df
