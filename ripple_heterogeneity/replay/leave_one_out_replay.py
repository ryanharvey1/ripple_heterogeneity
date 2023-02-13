from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs,
)
import pandas as pd
import numpy as np
import os
import pickle
from nelpy.analysis import replay
import copy
import random
import logging

logging.getLogger().setLevel(logging.ERROR)


def leave_one_out_score(results, direction, idx_replay):
    """
    Leave one out score for a given direction.
    Inputs:
        results: dict of results from replay_run.run_replay
        direction: string, 'inbound_epochs' or 'outbound_epochs'
    Outputs:
        score: float, score (cell x replay) for the given direction
    """
    scores = []
    units = np.array(results[direction]["bst_placecells"].series_ids).T
    for i in units:

        bst = copy.deepcopy(results[direction]["bst_placecells"])
        tc = copy.deepcopy(results[direction]["tc"])

        tc = tc._unit_subset(units[units != i])
        bst = bst._unit_subset(units[units != i])

        scores.append(
            replay.trajectory_score_bst(
                bst[idx_replay], tc, w=3, n_shuffles=0, normalize=True
            )
        )
    return np.array(scores)


def get_leave_one_out_score(results, direction):
    """
    Get leave one out score for a given direction.
    Inputs:
        results: dict of results from replay_run.run_replay
        direction: string, 'inbound_epochs' or 'outbound_epochs'
    Outputs:
        left_out_score: float, score (n cell) for the given direction
        left_out_score_cell_count_norm: float, normed score (n cell) for the given direction
        len(idx_replay): int, number of replays
    """
    # locate significant replay events
    idx_replay = np.where(results[direction]["df"].score_pval_time_swap < 0.05)[0]
    if len(idx_replay) == 0:
        return np.nan, np.nan, np.nan
    # find observed score for each replay event
    obs_scores = replay.trajectory_score_bst(
        results[direction]["bst_placecells"][idx_replay],
        results[direction]["tc"],
        w=3,
        n_shuffles=0,
        normalize=True,
    )

    # find score for each replay event and leave one cell out for each replay
    scores = leave_one_out_score(results, direction, idx_replay)
    left_out_score = (scores - obs_scores)
    # get number of cells in each replay
    n_cells_per_replay = [
        bst.n_active for bst in results[direction]["bst_placecells"][idx_replay]
    ]
    # normed score by number of cells in each replay
    left_out_score_cell_count_norm = (
        scores * n_cells_per_replay - obs_scores * n_cells_per_replay
    )

    return left_out_score, left_out_score_cell_count_norm, len(idx_replay)

def run(basepath:str,replay_save_path: str = None)-> pd.core.frame.DataFrame:

    # locate saved replay file and load it
    save_file = os.path.join(
        replay_save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )

    with open(save_file, "rb") as f:
        results = pickle.load(f)

    if results is None:
        return None
    
    epoch_df = loading.load_epoch(basepath)
    epoch_df = compress_repeated_epochs.main(epoch_df, epoch_name="sleep")

    results_df = pd.DataFrame()

    for direction in ["inbound_epochs", "outbound_epochs"]:
        if results[direction] == {}:
            continue
        (
            left_out_score,
            left_out_score_cell_count_norm,
            n_replay
        ) = get_leave_one_out_score(results, direction)

        if left_out_score is None or np.isnan(left_out_score).all():
            continue

        temp_df = pd.DataFrame()
        temp_df["left_out_score"] = left_out_score.T.flatten()
        temp_df["left_out_score_cell_count_norm"] = left_out_score_cell_count_norm.T.flatten()

        temp_df["replay_n"] = (
            (np.ones((left_out_score.shape[0], left_out_score.shape[1])) * np.arange(left_out_score.shape[1]))
            .T.astype(int)
            .ravel()
        )

        # add metadata from cell metrics
        cell_metrics = results[direction]["cell_metrics"]
        temp_df["UID"] = np.tile(cell_metrics.UID.values, left_out_score.shape[1])
        temp_df["deepSuperficialDistance"] = np.tile(cell_metrics.deepSuperficialDistance.values, left_out_score.shape[1])
        temp_df["deepSuperficial"] = np.tile(cell_metrics.deepSuperficial.values, left_out_score.shape[1])
        temp_df["brainRegion"] = np.tile(cell_metrics.brainRegion.values, left_out_score.shape[1])

        # add metadata from replay_df
        current_replay_df = results[direction]["df"].query(
            "score_pval_time_swap < 0.05"
        )
        keys = current_replay_df.keys()
        for replay_i, df in enumerate(current_replay_df.iterrows()):
            for key in keys:
                temp_df.loc[temp_df.replay_n == replay_i, key] = df[1][key]


        # add epoch metadata
        for t in range(epoch_df.shape[0]):
            idx = (temp_df.start >= epoch_df.startTime.iloc[t]) & (
                temp_df.stop <= epoch_df.stopTime.iloc[t]
            )
            temp_df.loc[idx, "epochs"] = epoch_df.name.iloc[t]
            temp_df.loc[idx, "environment"] = epoch_df.environment.iloc[t]
            temp_df.loc[idx, "epoch_n"] = t

        temp_df["direction"] = direction
        temp_df["basepath"] = basepath
        results_df = pd.concat([results_df, temp_df], ignore_index=True)

    if results_df.shape[0] == 0:
        return None

    results_df = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(results_df)

    return results_df

