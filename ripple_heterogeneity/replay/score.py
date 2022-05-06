import numpy as np
from nelpy.analysis import replay
from nelpy.decoding import decode1D as decode
from nelpy.decoding import get_mode_pth_from_array
import multiprocessing
from joblib import Parallel, delayed


def shuffle_and_score(posterior_array, w, normalize, tc):

    posterior_ts = replay.time_swap_array(posterior_array)
    posterior_cs = replay.column_cycle_array(posterior_array)

    scores_time_swap = replay.trajectory_score_array(
        posterior=posterior_ts, w=w, normalize=normalize
    )
    scores_col_cycle = replay.trajectory_score_array(
        posterior=posterior_cs, w=w, normalize=normalize
    )
    diff_mode_pth_cs = np.diff(get_mode_pth_from_array(posterior_cs, tuningcurve=tc))

    return (
        scores_time_swap,
        scores_col_cycle,
        np.abs(np.nanmean(diff_mode_pth_cs)),
    )


def trajectory_score_bst(
    bst, tuningcurve, w=None, n_shuffles=1000, weights=None, normalize=False
):
    """Compute the trajectory scores from Davidson et al. for each event
    in the BinnedSpikeTrainArray.

    This function returns the trajectory scores by decoding all the
    events in the BinnedSpikeTrainArray, and then calling an external
    function to determine the slope and intercept for each event, and
    then finally computing the scores for those events.

    If n_shuffles > 0, then in addition to the trajectory scores,
    shuffled scores will be returned for both column cycle shuffling, as
    well as posterior time bin shuffling (time swap).

    NOTE1: this function does NOT attempt to find the line that
    maximizes the trajectory score. Instead, it delegates the
    determination of the line to an external function (which currently
    is called from trajectory_score_array), and at the time of writing
    this documentation, is simply the best line fit to the modes of the
    decoded posterior distribution.

    NOTE2: the score is then the sum of the probabilities in a band of
    w bins around the line, ignoring bins that are NaNs. Even when w=0
    (only the sum of peak probabilities) this is different from the r^2
    coefficient of determination, in that here more concentrated
    posterior probabilities contribute more heavily than weaker ones.

    NOTE3: the returned scores are NOT normalized, but if desired, they
    can be normalized by dividing by the number of non-NaN bins in each
    event.

    Reference(s)
    ------------
    Davidson TJ, Kloosterman F, Wilson MA (2009)
        Hippocampal replay of extended experience. Neuron 63:497-507

    Parameters
    ----------
    bst : BinnedSpikeTrainArray
        BinnedSpikeTrainArray containing all the candidate events to
        score.
    tuningcurve : TuningCurve1D
        Tuning curve to decode events in bst.
    w : int, optional (default is 0)
        Half band width for calculating the trajectory score. If w=0,
        then only the probabilities falling directly under the line are
        used. If w=1, then a total band of 2*w+1 = 3 will be used.
    n_shuffles : int, optional (default is 250)
        Number of times to perform both time_swap and column_cycle
        shuffles.
    weights : not yet used, but idea is to assign weights to the bands
        surrounding the line
    normalize : bool, optional (default is False)
        If True, the scores will be normalized by the number of non-NaN
        bins in each event.

    Returns
    -------
    scores, [scores_time_swap, scores_col_cycle]
        scores is of size (bst.n_epochs, )
        scores_time_swap and scores_col_cycle are each of size
            (n_shuffles, bst.n_epochs)
    """

    if w is None:
        w = 0
    if not float(w).is_integer:
        raise ValueError("w has to be an integer!")

    if float(n_shuffles).is_integer:
        n_shuffles = int(n_shuffles)
    else:
        raise ValueError("n_shuffles must be an integer!")

    posterior, bdries, _, _ = decode(bst=bst, ratemap=tuningcurve)

    # idea: cycle each column so that the top w rows are the band
    # surrounding the regression line

    scores = np.zeros(bst.n_epochs)
    avg_jump = np.zeros(bst.n_epochs)

    if n_shuffles > 0:
        scores_time_swap = np.zeros((n_shuffles, bst.n_epochs))
        scores_col_cycle = np.zeros((n_shuffles, bst.n_epochs))
        jump_time_swap = np.zeros((n_shuffles, bst.n_epochs))
        jump_col_cycle = np.zeros((n_shuffles, bst.n_epochs))

    num_cores = multiprocessing.cpu_count()

    for idx in range(bst.n_epochs):
        posterior_array = posterior[:, bdries[idx] : bdries[idx + 1]]
        scores[idx] = replay.trajectory_score_array(
            posterior=posterior_array, w=w, normalize=normalize
        )
        avg_jump[idx] = np.abs(
            np.nanmean(
                np.diff(
                    get_mode_pth_from_array(posterior_array, tuningcurve=tuningcurve)
                )
            )
        )

        (
            scores_time_swap[:, idx],
            scores_col_cycle[:, idx],
            jump_col_cycle[:, idx],
        ) = zip(
            *Parallel(n_jobs=num_cores)(
                delayed(shuffle_and_score)(posterior_array, w, normalize, tuningcurve)
                for i in range(n_shuffles)
            )
        )

    if n_shuffles > 0:
        return (
            scores,
            avg_jump,
            scores_time_swap,
            scores_col_cycle,
            jump_col_cycle,
        )
    return scores, avg_jump
