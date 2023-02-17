import numpy as np
from nelpy.analysis import replay
from nelpy.decoding import decode1D as decode
from nelpy.decoding import get_mode_pth_from_array
import multiprocessing
from joblib import Parallel, delayed

# from replay_trajectory_classification.standard_decoder import detect_line_with_radon


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

    return _corr(time[:, np.newaxis],place_bin_centers[np.newaxis, :], posterior.T)


def shuffle_and_score(posterior_array, w, normalize, tc, ds, dp):

    posterior_ts = replay.time_swap_array(posterior_array)
    posterior_cs = replay.column_cycle_array(posterior_array)

    scores_time_swap = replay.trajectory_score_array(
        posterior=posterior_ts, w=w, normalize=normalize
    )
    scores_col_cycle = replay.trajectory_score_array(
        posterior=posterior_cs, w=w, normalize=normalize
    )

    weighted_corr_time_swap = weighted_correlation(posterior_ts)
    weighted_corr_col_cycle = weighted_correlation(posterior_cs)

    return (
        scores_time_swap,
        scores_col_cycle,
        weighted_corr_time_swap,
        weighted_corr_col_cycle,
    )


def trajectory_score_bst(
    bst,
    tuningcurve,
    w=None,
    n_shuffles=1000,
    weights=None,
    normalize=False,
    parallel=True,
):

    if w is None:
        w = 0
    if not float(w).is_integer:
        raise ValueError("w has to be an integer!")

    if float(n_shuffles).is_integer:
        n_shuffles = int(n_shuffles)
    else:
        raise ValueError("n_shuffles must be an integer!")

    posterior, bdries, _, _ = decode(bst=bst, ratemap=tuningcurve)

    scores = np.zeros(bst.n_epochs)
    weighted_corr = np.zeros(bst.n_epochs)

    if n_shuffles > 0:
        scores_time_swap = np.zeros((n_shuffles, bst.n_epochs))
        scores_col_cycle = np.zeros((n_shuffles, bst.n_epochs))
        weighted_corr_time_swap = np.zeros((n_shuffles, bst.n_epochs))
        weighted_corr_col_cycle = np.zeros((n_shuffles, bst.n_epochs))

    if parallel:
        num_cores = multiprocessing.cpu_count()

    ds, dp = bst.ds, np.diff(tuningcurve.bins)[0]

    for idx in range(bst.n_epochs):
        posterior_array = posterior[:, bdries[idx] : bdries[idx + 1]]
        scores[idx] = replay.trajectory_score_array(
            posterior=posterior_array, w=w, normalize=normalize
        )
        weighted_corr[idx] = weighted_correlation(posterior_array)

        if parallel:
            (
                scores_time_swap[:, idx],
                scores_col_cycle[:, idx],
                weighted_corr_time_swap[:, idx],
                weighted_corr_col_cycle[:, idx],
            ) = zip(
                *Parallel(n_jobs=num_cores)(
                    delayed(shuffle_and_score)(
                        posterior_array, w, normalize, tuningcurve, ds, dp
                    )
                    for _ in range(n_shuffles)
                )
            )
        else:
            (
                scores_time_swap[:, idx],
                scores_col_cycle[:, idx],
                weighted_corr_time_swap[:, idx],
                weighted_corr_col_cycle[:, idx],
            ) = zip(
                *[
                    shuffle_and_score(
                        posterior_array, w, normalize, tuningcurve, ds, dp
                    )
                    for _ in range(n_shuffles)
                ]
            )

    if n_shuffles > 0:
        return (
            scores,
            weighted_corr,
            scores_time_swap,
            scores_col_cycle,
            weighted_corr_time_swap,
            weighted_corr_col_cycle,
        )
    return scores, weighted_corr
