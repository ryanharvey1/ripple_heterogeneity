from ripple_heterogeneity.utils import (
    functions,
    loading,
    batch_analysis,
    add_new_deep_sup,
    compress_repeated_epochs,
)
from ripple_heterogeneity.assembly import assembly_reactivation
from ripple_heterogeneity.readout import assembly_multi_region
import pandas as pd
import numpy as np
import os
import nelpy as nel
import glob
import pickle
import itertools
from ripple_heterogeneity.place_cells import maps
from scipy import stats
import pynapple as nap


def locate_task_epoch(m1, env):
    epoch_df = m1.epoch_df.reset_index(drop=True).copy()
    epoch_df["duration"] = epoch_df.stopTime - epoch_df.startTime
    epoch_df.sort_values("duration", ascending=False, inplace=True)
    return int(epoch_df[epoch_df.environment.str.contains(env)].index[0])


def get_pos(basepath, m1, task_idx):
    position_df = loading.load_animal_behavior(basepath)
    position_df_no_nan = position_df.query("not x.isnull() & not y.isnull()")

    pos = nel.PositionArray(
        data=position_df_no_nan["linearized"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    # make min pos 2
    pos._data = (pos.data - np.nanmin(pos.data)) + 2

    pos = pos[m1.epochs[task_idx]]

    # get outbound and inbound epochs
    outbound_epochs, inbound_epochs = functions.get_linear_track_lap_epochs(
        pos.abscissa_vals, pos.data[0], newLapThreshold=20
    )
    return pos, outbound_epochs, inbound_epochs


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",
    putativeCellType="Pyr",
    weight_dt=0.1,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03,  # dt in seconds for binning st to get activation strength
    verbose=False,  # print out progress
    env="linear",
    s_binsize=3,
    sigma=3,
    speed_thres=4,
):

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

    # find longest xx session
    if not m1.epoch_df.environment.str.contains(env).any():
        return None
    task_idx = locate_task_epoch(m1, env)

    m1.get_weights(m1.epochs[task_idx])
    _, assembly_df = assembly_multi_region.compile_results_df({"react": m1})

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

    pos, outbound_epochs, inbound_epochs = get_pos(basepath, m1, task_idx)

    ext_xmin, ext_xmax = (
        np.floor(pos.data[0].min() / 10) * 10,
        np.ceil(pos.data[0].max() / 10) * 10,
    )
    n_extern = int((ext_xmax - ext_xmin) / s_binsize)

    label_df = pd.DataFrame()

    direction_label = ["outbound_epochs", "inbound_epochs"]

    for dir_epoch_i, dir_epoch in enumerate([outbound_epochs, inbound_epochs]):

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
        tc_ = tc_.rolling(window=5, win_type="gaussian", center=True).mean(std=sigma)

        # store tuning curves
        tc = pd.concat([tc, tc_], axis=1, ignore_index=True)

        # store associated metadata
        label_df_ = pd.DataFrame()
        label_df_["assembly_n"] = np.arange(assembly_act_task.n_signals).astype(int)
        label_df_["direction_label"] = direction_label[dir_epoch_i]
        label_df_.merge(counts_df.query("is_member_sig"))

        label_df = pd.concat([label_df, label_df_], ignore_index=True)

    results = {
        "tc": tc,
        "label_df": label_df,
        "assembly_act_task": assembly_act_task,
        "react": m1,
    }

    return results

def load_results():
    pass
