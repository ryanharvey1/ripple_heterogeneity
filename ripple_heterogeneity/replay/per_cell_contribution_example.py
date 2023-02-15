from ripple_heterogeneity.utils import (
    functions,
    add_new_deep_sup,
)
import pandas as pd
import numpy as np
import os
import pickle
import nelpy as nel
import copy
import os
from nelpy.analysis import replay


def get_scores(bst, tc, n_shuffles=400):

    # find observed score for each replay event
    scores, scores_time_swap, _ = replay.trajectory_score_bst(
        bst,
        tc,
        w=3,
        n_shuffles=n_shuffles,
        normalize=True,
    )
    _, _, std = functions.get_significant_events(scores, scores_time_swap)

    std = np.nanmean(std)
    scores = np.nanmean(scores)

    return std, scores


def run(
    basepath: str, replay_save_path: str, n_shuffles: int = 400
) -> pd.core.frame.DataFrame:

    # locate saved replay file and load it
    save_file = os.path.join(
        replay_save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )

    if not os.path.exists(save_file):
        return None

    # load replay results
    with open(save_file, "rb") as f:
        results = pickle.load(f)

    if results is None:
        return None
    
    # iterate over the two template types
    results_df = pd.DataFrame()
    for direction in ["inbound_epochs", "outbound_epochs"]:

        if results[direction] == {}:
            continue

        # locate significant replay events
        idx_replay = np.where(results[direction]["df"].score_pval_time_swap < 0.05)[0]

        # get new copy of bst and tc with all neurons
        bst = copy.deepcopy(results[direction]["bst_placecells"])
        tc = copy.deepcopy(results[direction]["tc"])

        # find observed score for each replay event
        std_obs, obs_scores = get_scores(bst[idx_replay], tc, n_shuffles=n_shuffles)

        # pull out cell metrics for use
        cell_metrics = results[direction]["cell_metrics"]
        cell_metrics = add_new_deep_sup.deep_sup_from_distance(cell_metrics)

        # iterate over deep and superficial
        # in this loop, cells will be dropped one by one and replay quality will be assessed each time
        for deepsup in ["Deep", "Superficial"]:
            
            # get new copy of bst and tc with all neurons
            bst = copy.deepcopy(results[direction]["bst_placecells"])
            tc = copy.deepcopy(results[direction]["tc"])

            # get current series ids for particular deep/sup
            units_ = np.array(bst.series_ids)[
                cell_metrics.deepSuperficial.values == deepsup
            ]
            scores = []
            n_dropped = []
            std_dropped = []
            for i, units_i in enumerate(units_):
                # get current unit series ids (will decrease each iteration)
                units = np.array(bst.series_ids).T

                # note, middle layer will always remain as we are only removing deep/sup
                tc = tc._unit_subset(units[units != units_i])
                bst = bst._unit_subset(units[units != units_i])

                # save n dropped, plus 1 as we drop 1 in first iter
                n_dropped.append(i + 1)

                # if we have removed all from bst, continue
                if bst.isempty:
                    scores.append(np.nan)
                    std_dropped.append(np.nan)
                    continue

                # find score for each replay event
                std_dropped_, scores_ = get_scores(
                    bst[idx_replay], tc, n_shuffles=n_shuffles
                )

                # get mean scores over replays
                std_dropped.append(np.nanmean(std_dropped_))
                scores.append(np.nanmean(scores_))

            temp_df = pd.DataFrame()
            temp_df["scores_dropped"] = scores
            temp_df["std_dropped"] = std_dropped
            temp_df["scores_obs"] = obs_scores
            temp_df["std_obs"] = std_obs
            temp_df["n_dropped"] = n_dropped
            temp_df["UID"] = cell_metrics.query("deepSuperficial==@deepsup").UID.values
            temp_df["deepSuperficial"] = cell_metrics.query("deepSuperficial==@deepsup").deepSuperficial.values
            temp_df["putativeCellType"] = cell_metrics.query("deepSuperficial==@deepsup").putativeCellType.values
            temp_df["total_n_units"] = cell_metrics.shape[0]
            temp_df["n_deep"] = cell_metrics.query("deepSuperficial=='Deep'").shape[0]
            temp_df["n_sup"] = cell_metrics.query("deepSuperficial=='Superficial'").shape[0]
            temp_df["n_middle"] = cell_metrics.query("deepSuperficial=='middle'").shape[0]
            temp_df["direction"] = direction

            results_df = pd.concat([results_df, temp_df], ignore_index=True)

    results_df["basepath"] = basepath

    return results_df
