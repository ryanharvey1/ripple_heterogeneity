from ripple_heterogeneity.utils import (
    loading,
    add_new_deep_sup,
)
import numpy as np
import nelpy as nel
import logging
logging.getLogger().setLevel(logging.ERROR)

def get_pos(basepath):
    position_df = loading.load_animal_behavior(basepath)

    if position_df is None:
        return None
    if position_df.shape[0] == 0:
        return None

    if "timestamps" not in position_df.columns and "time" in position_df.columns:
        position_df["timestamps"] = position_df.time.values

    # select position data
    if "x" in position_df.columns and "y" in position_df.columns:
        pos = nel.PositionArray(
            data=position_df[["x", "y"]].values.T,
            timestamps=position_df.timestamps.values,
        )
    elif "x" in position_df.columns:
        pos = nel.PositionArray(
            data=position_df[["x"]].values.T,
            timestamps=position_df.timestamps.values,
        )
    elif "linearized" in position_df.columns:
        pos = nel.PositionArray(
            data=position_df[["linearized"]].values.T,
            timestamps=position_df.timestamps.values,
        )
    else:
        return None
    return pos

def run(basepath, speed_thres=4):

    st, cell_metrics = loading.load_spikes(
        basepath,
        brainRegion="CA1",
        putativeCellType="Pyr",
    )
    if st.isempty:
        return None

    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)

    pos = get_pos(basepath)
    
    if pos is None:
        return None
    if np.isnan(pos.data).all():
        return None

    # calculate speed
    speed = nel.utils.ddt_asa(pos, smooth=True, sigma=0.1, norm=True)
    # get inactive epochs
    inactive_epochs = nel.utils.get_inactive_epochs(
        speed, v1=speed_thres, v2=speed_thres
    )
    if inactive_epochs.isempty:
        return None
    # get ripple epochs
    ripple_epochs = loading.load_ripples_events(basepath, return_epoch_array=True)
    if ripple_epochs.isempty:
        return None

    # get session epochs for epoch support
    epoch_df = loading.load_epoch(basepath)
    beh_epochs = nel.EpochArray(
        np.array([epoch_df.startTime.values, epoch_df.stopTime.values]).T
    )
    if beh_epochs.isempty:
        return None
    support = nel.EpochArray(
        np.array([epoch_df.startTime.values[0], epoch_df.stopTime.values[-1]]).T
    )
    inactive_epochs.domain = support
    ripple_epochs.domain = support

    # get behavior epochs that are not sleep
    beh_epochs = beh_epochs[~epoch_df.environment.str.contains("sleep")]
    if beh_epochs.isempty:
        return None

    # create immobile and mobile epochs that our outside of ripples
    immobile_non_ripple = inactive_epochs & ~ripple_epochs & beh_epochs
    mobile_non_ripple = ~inactive_epochs & ~ripple_epochs & beh_epochs

    if immobile_non_ripple.isempty | mobile_non_ripple.isempty:
        return None

    # get spike trains for immobile and mobile epochs
    st_immobile = st[immobile_non_ripple]
    st_mobile = st[mobile_non_ripple]

    # calculate firing rates
    immobile_fr = st_immobile.n_events / immobile_non_ripple.length
    mobile_fr = st_mobile.n_events / mobile_non_ripple.length

    # create results dataframe with results and metadata
    results = cell_metrics[
        ["UID", "brainRegion", "putativeCellType", "deepSuperficial"]
    ].copy()
    results["immobile_fr"] = immobile_fr
    results["mobile_fr"] = mobile_fr
    results["fr_diff_ratio_a"] = (immobile_fr - mobile_fr) / (immobile_fr + mobile_fr)
    results["fr_diff_ratio_b"] = (immobile_fr - mobile_fr) / mobile_fr
    results["fr_diff_ratio_c"] = immobile_fr / mobile_fr
    results["fr_diff"] = immobile_fr - mobile_fr
    results["basepath"] = basepath

    return results
