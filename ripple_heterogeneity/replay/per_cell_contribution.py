"""analysis to do per cell contribution analysis from Grosmark & Buzsaki 2016"""
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs
)
import pandas as pd
import numpy as np
import glob
import os
import pickle
import nelpy as nel
from nelpy.analysis import replay
import copy
import random


def weighted_correlation(posterior, time=None, place_bin_centers=None):
    def _m(x, w):
        """Weighted Mean"""
        return np.sum(x * w) / np.sum(w)

    def _cov(x, y, w):
        """Weighted Covariance"""
        return np.sum(w * (x - _m(x, w)) * (y - _m(y, w))) / np.sum(w)

    def _corr(x, y, w):
        """Weighted Correlation"""
        return _cov(x, y, w) / np.sqrt(_cov(x, x, w) * _cov(y, y, w))

    if time is None:
        time = np.arange(posterior.shape[1])
    if place_bin_centers is None:
        place_bin_centers = np.arange(posterior.shape[0])

    place_bin_centers = place_bin_centers.squeeze()
    posterior[np.isnan(posterior)] = 0.0

    return _corr(time[:, np.newaxis], place_bin_centers[np.newaxis, :], posterior.T)


def circular_place_field_shuffle(tuningcurve):
    """
    circularly shift each place field individually by a random amount
    """
    # make copy of tuning curves so manipulation doesn't effect original
    tuningcurve_new = copy.deepcopy(tuningcurve)
    # make random ints to circularly shift tuning curves
    shuff_idx = random.sample(
        range(tuningcurve_new.ratemap.shape[1]), tuningcurve_new.ratemap.shape[0]
    )
    # randomly shift each tuning curve
    x = [
        np.roll(tc, shuff_val)
        for tc, shuff_val in zip(tuningcurve_new.ratemap, shuff_idx)
    ]
    # add data back to copy
    tuningcurve_new._ratemap = np.array(x)

    return tuningcurve_new


def weighted_correlation_score_bst(bst, tuningcurve, n_shuffles=1000):

    posterior, bdries, mode_pth, mean_pth = replay.decode(bst=bst, ratemap=tuningcurve)
    scores = np.zeros(bst.n_epochs)
    if n_shuffles > 0:
        scores_circ_place_shuff = np.zeros((n_shuffles, bst.n_epochs))

    for idx in range(bst.n_epochs):
        posterior_array = posterior[:, bdries[idx] : bdries[idx + 1]]
        scores[idx] = weighted_correlation(posterior_array)

        for shflidx in range(n_shuffles):

            posterior_shuff, _, _, _ = replay.decode(
                bst=bst[idx], ratemap=circular_place_field_shuffle(tuningcurve)
            )
            scores_circ_place_shuff[shflidx, idx] = weighted_correlation(
                posterior_shuff
            )
    if n_shuffles > 0:
        return scores, scores_circ_place_shuff
    else:
        return scores


def shuffle_and_score_single_cell(tc_new, bst, cell_id, n_shuffles=1000):

    scores_circ_place_shuff = np.zeros((n_shuffles, bst.n_epochs))

    for shflidx in range(n_shuffles):
        shuff_idx = random.sample(range(tc_new.ratemap.shape[1]), 1)
        tc_new._ratemap[cell_id, :] = np.roll(tc_new.ratemap[cell_id, :], shuff_idx)

        scores_circ_place_shuff[shflidx, :] = weighted_correlation_score_bst(
            bst, tc_new, n_shuffles=0
        )

    return np.array(scores_circ_place_shuff)


def get_pcc_score(
    results, direction, n_shuffles_single_cell=1000, n_shuffles_corr=1000
):
    """
    get_pcc_score: will return a (n cells x n events) matrix of contribution of that cell to the event
    """

    # pull out data
    # only look at pre-defined replay events
    idx_replay = np.where(results[direction]["df"].score_pval_time_swap < 0.05)[0]
    bst = copy.deepcopy(results[direction]["bst_placecells"][idx_replay])
    tc = copy.deepcopy(results[direction]["tc"])

    # get the number of units that participated in each event
    n_active = [bst_.n_active for bst_ in bst]
    n_active = np.array(n_active)

    # get observed scores
    scores, scores_circ_place_shuff = weighted_correlation_score_bst(
        bst, tc, n_shuffles=n_shuffles_corr
    )
    _, pval, rZ_obs = functions.get_significant_events(scores, scores_circ_place_shuff)

    pcc = []
    # iter over each cell
    for cell_id in range(len(bst.series_labels)):
        tc_new = copy.deepcopy(tc)

        # shuffle single cell and return shuffled weighted correlations
        scores_circ_place_shuff = shuffle_and_score_single_cell(
            tc_new,
            bst,
            cell_id,
            n_shuffles=n_shuffles_single_cell,
        )
        _, pval, rZ_shuff = functions.get_significant_events(
            scores, scores_circ_place_shuff
        )

        # calculate per cell contribution
        pcc.append((rZ_obs - rZ_shuff) * n_active)

    pcc = np.array(pcc)

    # find events that the cells did not participation in and make pcc nan
    events_active = [np.any(bst_.data > 0, axis=1) for bst_ in bst]
    events_active = np.array(events_active).T
    pcc[events_active == False] = np.nan

    return pcc


def run(
    basepath: str,
    replay_save_path: str = None,
    n_shuffles_single_cell: int = 1000,
    n_shuffles_corr: int = 1000,
) -> pd.core.frame.DataFrame:

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

    # saved replay files have two major keys associated with the direction of travel
    # iterate over them
    results_df = pd.DataFrame()
    for direction in ["inbound_epochs", "outbound_epochs"]:

        pcc = get_pcc_score(
            results,
            direction,
            n_shuffles_single_cell=n_shuffles_single_cell,
            n_shuffles_corr=n_shuffles_corr,
        )
        # add to dataframe
        temp_df = pd.DataFrame()
        temp_df["pcc"] = pcc.T.flatten()
        temp_df["replay_n"] = (
            (np.ones((pcc.shape[0], pcc.shape[1])) * np.arange(pcc.shape[1]))
            .T.astype(int)
            .ravel()
        )
        # add metadata from cell metrics
        cell_metrics = results[direction]["cell_metrics"]
        temp_df["UID"] = np.tile(cell_metrics.UID.values, pcc.shape[1])
        temp_df["deepSuperficialDistance"] = np.tile(
            cell_metrics.deepSuperficialDistance.values, pcc.shape[1]
        )
        temp_df["deepSuperficial"] = np.tile(
            cell_metrics.deepSuperficial.values, pcc.shape[1]
        )
        temp_df["brainRegion"] = np.tile(cell_metrics.brainRegion.values, pcc.shape[1])

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
            idx = (temp_df.start >= epoch_df.startTime.iloc[t]) & (temp_df.stop <= epoch_df.stopTime.iloc[t])
            temp_df.loc[idx,'epochs'] = epoch_df.name.iloc[t] 
            temp_df.loc[idx,'environment'] = epoch_df.environment.iloc[t] 
            temp_df.loc[idx,'epoch_n'] = t

        temp_df["direction"] = direction
        temp_df["basepath"] = basepath

        results_df = pd.concat([results_df, temp_df], ignore_index=True)

    results_df = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(results_df)

    return results_df
