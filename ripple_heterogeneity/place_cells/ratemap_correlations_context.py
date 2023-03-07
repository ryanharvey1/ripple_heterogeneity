"""This script finds the correlation between ratemaps for each cell for each task type."""

from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
)
import nelpy as nel
from ripple_heterogeneity.place_cells import maps
import pandas as pd
import numpy as np
from itertools import combinations
from typing import Union

def find_repeated_strings(lst: list) -> dict:
    """
    Find all repeated strings in a list.
    Parameters
    ----------
    lst : list or list-like object
        List of strings to search for repeated strings.
    Returns
    -------
    repeated_strings : dict
        Dictionary with repeated strings as keys and list of indices where
        the string is found as values.
    """
    if not isinstance(lst, list):
        try:
            lst = list(lst)
        except:
            raise TypeError("lst must be a list or a list-like object")
    repeated_strings = {}
    for i in range(len(lst)):
        if lst[i] in repeated_strings:
            repeated_strings[lst[i]].append(i)
        elif lst.count(lst[i]) > 1:
            repeated_strings[lst[i]] = [i]
    return repeated_strings


def ratemap_cor(A: np.ndarray, B: np.ndarray) -> float:
    # for given ratemaps A & B, find their correlation
    x1 = A.flatten()
    x2 = B.flatten()
    bad_idx = np.isnan(x1) | np.isnan(x2)
    return np.corrcoef(x1[~bad_idx], x2[~bad_idx])[0, 1]


def run(basepath: str) -> Union[pd.DataFrame,None]:

    position_df = loading.load_animal_behavior(basepath)
    position_df_no_nan = position_df.query("not x.isnull() & not y.isnull()")

    if position_df_no_nan.shape[0] == 0:
        return pd.DataFrame()

    st, cm = loading.load_spikes(basepath, putativeCellType="Pyr", brainRegion="CA1")

    if st.isempty:
        return pd.DataFrame()

    cm = add_new_deep_sup.deep_sup_from_distance(cm)

    pos = nel.PositionArray(
        data=position_df_no_nan[["x", "y"]].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    beh_df = loading.load_epoch(basepath)
    beh_df = beh_df.query("environment != 'sleep'")

    beh_epochs = nel.EpochArray(np.array([beh_df.startTime, beh_df.stopTime]).T)

    # find all repeated strings in beh_df.environment
    tasks = find_repeated_strings(list(beh_df.environment.values))

    results_df = pd.DataFrame()
    for task in tasks:

        results_df_temp = pd.DataFrame()
        # linearize position if necessary
        if "linear" in task:
            for i in tasks[task]:
                idx = functions.in_intervals(pos.abscissa_vals, beh_epochs[i].data)
                x, y = pos.data[:, idx]
                pos._data[:, idx] = functions.linearize_position(x, y)

        # pull out coordinates for for this task type to get spatial extent
        current_pos = pos[beh_epochs[tasks[task]]]
        if current_pos.isempty:
            continue
        xminmax = [np.nanmin(current_pos.data[0, :]), np.nanmax(current_pos.data[0, :])]
        yminmax = [np.nanmin(current_pos.data[1, :]), np.nanmax(current_pos.data[1, :])]

        # find ratemaps for each cell
        # for reduced complexity, using 2d spatial maps and not splitting by direction (linear track)
        ratemaps = []
        spatial_pval = []
        for i in tasks[task]:

            if pos[beh_epochs[i]].isempty | st[beh_epochs[i]].isempty:
                continue

            spatial_maps = maps.SpatialMap(
                pos[beh_epochs[i]],
                st[beh_epochs[i]],
                dim=2,
                s_binsize=3,
                x_minmax=xminmax,
                y_minmax=yminmax,
                tuning_curve_sigma=3,
                n_shuff=500,
                parallel_shuff=True,
            )
            # store ratemaps
            ratemaps.append(spatial_maps.ratemap)
            # store p-values for spatial information shuffle
            spatial_pval.append(spatial_maps.shuffle_spatial_information())

        # find correlation between all pairs of ratemaps for each cell
        # find all pairs of ratemaps
        pairs = np.array(list(combinations(np.arange(0, len(ratemaps)), 2)))
        # set up array to store correlations
        within_corr = np.zeros((len(st.data), len(pairs)))
        # loop through pairs
        for pair_i, s in enumerate(pairs):
            # loop through cells
            for cell_i in range(len(st.data)):
                # find correlation between ratemaps
                within_corr[cell_i, pair_i] = ratemap_cor(
                    ratemaps[s[0]][cell_i, :, :], ratemaps[s[1]][cell_i, :, :]
                )
        # take mean across pairs
        within_corr = np.nanmean(within_corr, axis=1)
        results_df_temp["corr"] = within_corr
        results_df_temp["n_pass_spatial_shuff"] = np.sum(
            np.vstack(spatial_pval) < 0.05, axis=0
        )
        results_df_temp["UID"] = cm.UID.values
        results_df_temp["deepSuperficial"] = cm.deepSuperficial.values
        results_df_temp["task"] = task
        results_df = pd.concat([results_df, results_df_temp], ignore_index=True)

    return results_df