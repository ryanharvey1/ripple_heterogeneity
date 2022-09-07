from ripple_heterogeneity.utils import (
    functions,
    loading,
)
from ripple_heterogeneity.assembly import assembly_reactivation
from ripple_heterogeneity.readout import assembly_multi_region
import pandas as pd
import numpy as np
import os
import nelpy as nel
import glob
import pickle
import pynapple as nap


def locate_task_epoch(m1, env):
    epoch_df = m1.epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    epoch_df.sort_values("duration", ascending=False, inplace=True)
    return int(epoch_df[epoch_df.environment.str.contains(env)].index[0])


def dissociate_laps_by_states(states, dir_epoch, states_of_interest=[1, 2]):
    # unique_states = np.unique(states.data[~np.isnan(states.data)])
    lap_id = []
    for ep in dir_epoch:
        state_count = []
        for us in states_of_interest:
            state_count.append(np.nansum(states[ep].data == us))
        lap_id.append(states_of_interest[np.argmax(state_count)])
    return np.array(lap_id).astype(int)


def get_pos(basepath, m1, task_idx):
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
    pos._data = (pos.data - np.nanmin(pos.data)) + 2

    pos = pos[m1.epochs[task_idx]]

    states = nel.AnalogSignalArray(
        data=position_df_no_nan["states"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    states = states[m1.epochs[task_idx]]

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

    return pos, right_epochs, left_epochs, states, position_df_no_nan


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.1,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,  # dt in seconds for binning st to get activation strength
    verbose=False,  # print out progress
    env="Mwheel|Tmaze|tmaze",  # enviroment you want to look at (current should only be linear)
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

    # locate key points (TODO: fix the hard coded values)
    start_pos = nodes_and_edges["node_positions"][0]
    decision_pos = nodes_and_edges["node_positions"][1]
    reward_left_pos = nodes_and_edges["node_positions"][3]
    reward_right_pos = nodes_and_edges["node_positions"][5]

    m1 = assembly_reactivation.AssemblyReact(
        basepath,
        brainRegion=regions,
        putativeCellType=putativeCellType,
        weight_dt=weight_dt,
        z_mat_dt=z_mat_dt,
    )
    m1.load_data()

    # check if no cells were found
    if m1.cell_metrics.shape[0] == 0:
        return None
    # check for any ca1 cells
    if not m1.cell_metrics.brainRegion.str.contains("CA1").any():
        return None

    # check for any cortex cells
    if not m1.cell_metrics.brainRegion.str.contains(
        "PFC|EC1|EC2|EC3|EC4|EC5|MEC"
    ).any():
        return None

    # check if any env
    if not m1.epoch_df.environment.str.contains(env).any():
        return None

    # find longest xx session
    task_idx = locate_task_epoch(m1, env)

    m1.get_weights(m1.epochs[task_idx])

    # check if any assemblies
    if len(m1.patterns) == 0:
        return None

    _, assembly_df, keep_assembly = assembly_multi_region.compile_results_df(
        {"react": m1}
    )

    # check if any sig members
    if not assembly_df.is_member_sig.any():
        return None

    counts_df = (
        assembly_df.groupby(["assembly_n", "is_member_sig"])
        .apply(
            lambda x: pd.Series(
                {
                    "n_deep": (x.deepSuperficial == "Deep").sum(),
                    "n_sup": (x.deepSuperficial == "Superficial").sum(),
                    "n_mec": (
                        x.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC")
                    ).sum(),
                    "n_pfc": (x.brainRegion.str.contains("PFC")).sum(),
                }
            )
        )
        .reset_index()
    )

    assembly_act_task = m1.get_assembly_act(epoch=m1.epochs[task_idx])
    assembly_act_task._data = assembly_act_task.data[keep_assembly]

    pos, right_epochs, left_epochs, states, position_df_no_nan = get_pos(
        basepath, m1, task_idx
    )
    if pos is None:
        return

    # TODO: locate key locations in linear coords


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
            t=assembly_act_task[dir_epoch][run_epochs].abscissa_vals,
            d=assembly_act_task[dir_epoch][run_epochs].data.T,
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
        label_df_["assembly_n"] = np.arange(assembly_act_task.n_signals).astype(int)
        label_df_["direction_label"] = direction_label[dir_epoch_i]
        label_df_ = label_df_.merge(counts_df.query("is_member_sig"))

        label_df = pd.concat([label_df, label_df_], ignore_index=True)

    results = {
        "tc": tc,
        "label_df": label_df,
        "assembly_act_task": assembly_act_task,
        "react": m1,
        "pos": pos,
        "right_epochs": right_epochs,
        "left_epochs": left_epochs,
        "states": states,
        "task_epoch": m1.epochs[task_idx],
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
        results["label_df"]["basepath"] = results["react"].basepath
        label_df = pd.concat([label_df, results["label_df"]], ignore_index=True)

    return tc, label_df
