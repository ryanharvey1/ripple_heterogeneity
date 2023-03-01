from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
from ripple_heterogeneity.assembly import assembly_reactivation,find_sig_assembly
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


def locate_task_epoch(epoch_df, env):
    epoch_df = epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    epoch_df.sort_values("duration", ascending=False, inplace=True)
    return int(epoch_df[epoch_df.environment.str.contains(env)].index[0])


def get_pos(basepath, epochs, epoch_df, task_idx):

    position_df = loading.load_animal_behavior(basepath)

    if "linearized" not in position_df.columns:
        return None, None, None, None
    
    position_df_no_nan = position_df.query(
        "not x.isnull() & not y.isnull() & not linearized.isnull()"
    )

    if position_df_no_nan.shape[0] == 0:
        return None, None, None, None

    pos = nel.PositionArray(
        data=position_df_no_nan["linearized"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    # make min pos 2
    pos._data = (pos.data - np.nanmin(pos.data)) + 2

    pos = pos[epochs[task_idx]]

    # get outbound and inbound epochs
    outbound_epochs, inbound_epochs = functions.get_linear_track_lap_epochs(
        pos.abscissa_vals, pos.data[0], newLapThreshold=20
    )

    position_df_no_nan = position_df_no_nan[
        position_df_no_nan["time"].between(
            epoch_df.iloc[task_idx].startTime, epoch_df.iloc[task_idx].stopTime
        )
    ]
    return pos, outbound_epochs, inbound_epochs, position_df_no_nan


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
    if len(m1.patterns) == 0:
        return None, None

    assembly_act = m1.get_assembly_act(epoch=epoch)
    return m1, assembly_act


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    cross_regions=(("CA1", "PFC"), ("CA1", "EC1|EC2|EC3|EC4|EC5|MEC")),
    putativeCellType="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt=0.05,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,  # dt in seconds for binning st to get activation strength
    verbose=False,  # print out progress
    env="linear",  # enviroment you want to look at (current should only be linear)
    s_binsize=3,  # spatial bin size
    smooth_sigma=3,  # smoothing sigma in cm
    smooth_window=10,  # smoothing window in cm
    speed_thres=4,  # speed threshold for ratemap in cm/sec
):

    epoch_df = loading.load_epoch(basepath)
    epochs = nel.EpochArray(np.array([epoch_df.startTime, epoch_df.stopTime]).T)

    # check if any env
    if not epoch_df.environment.str.contains(env).any():
        return None

    # find longest xx session
    # task_idx = locate_task_epoch(m1, env)
    # position_df = loading.load_animal_behavior(basepath)
    # task_idx = int(np.where(epoch_df.name == position_df[~position_df.linearized.isnull()].epochs.unique()[0])[0][0])

    # find longest xx session
    task_idx = locate_task_epoch(epoch_df, env)

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
                epoch=epochs[task_idx],
            )
            if assembly_act_ is None:
                continue
            if m1_.n_assemblies() == 0:
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
                [deepSuperficial + "_" + cross_region[1]] * assembly_act_.n_signals
            )
            m1[deepSuperficial + "_" + cross_region[1]] = m1_

    if len(assembly_act) == 0:
        return None

    assembly_act = nel.AnalogSignalArray(
        data=np.vstack(assembly_act),
        timestamps=abscissa_vals[0],
    )
    assem_labels = np.hstack(assem_labels)

    pos, outbound_epochs, inbound_epochs, position_df_no_nan = get_pos(
        basepath, epochs, epoch_df, task_idx
    )
    if pos is None:
        return None

    ext_xmin, ext_xmax = (
        np.floor(pos.data[0].min() / 10) * 10,
        np.ceil(pos.data[0].max() / 10) * 10,
    )
    n_extern = int((ext_xmax - ext_xmin) / s_binsize)

    label_df = pd.DataFrame()
    tc = pd.DataFrame()

    direction_label = ["outbound_epochs", "inbound_epochs"]

    for dir_epoch_i, dir_epoch in enumerate([outbound_epochs, inbound_epochs]):

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

        label_df = pd.concat([label_df, label_df_], ignore_index=True)

    label_df["basepath"] = basepath

    results = {
        "tc": tc,
        "label_df": label_df,
        "assembly_act_task": assembly_act,
        "react": m1,
        "pos": pos,
        "outbound_epochs": outbound_epochs,
        "inbound_epochs": inbound_epochs,
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
