import numpy as np
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel


def shuffle_labels(labels):
    """
    Shuffle labels
    Input:
        labels: labels
    Output:
        shuffled_labels: shuffled labels
    """
    return np.random.permutation(labels).reshape(labels.shape)


def get_shuffled_labels(cell_metrics, replay_par_mat, ep, n_shuffles=1000):
    """
    Shuffle labels and get the number of deep and superficial cells
    Input:
        cell_metrics: cell metrics dataframe
        replay_par_mat: participation matrix
        ep: epoch number, current replay event
        n_shuffles: number of shuffles
    Output:
        n_deep: shuffled number of deep cells
        n_sup: shuffled number of superficial cells
    """
    n_deep = []
    n_sup = []
    for _ in range(n_shuffles):
        labels = shuffle_labels(cell_metrics.deepSuperficial)
        n_deep.append(sum(labels[replay_par_mat[:, ep] == 1] == "Deep"))
        n_sup.append(sum(labels[replay_par_mat[:, ep] == 1] == "Superficial"))
    return n_deep, n_sup


def get_significant_events(cell_metrics, replay_par_mat, n_shuffles=1000, q_perc=90):
    """
    Get the number of significant events
    Input:
        cell_metrics: cell metrics dataframe
        replay_par_mat: participation matrix
        n_shuffles: number of shuffles
        q_perc: percentile
    Output:
        sig_idx_deep: index of significant events in deep cells
        sig_idx_sup: index of significant events in superficial cells
    """
    # initialize variables
    n_deep = np.zeros([n_shuffles, replay_par_mat.shape[1]])
    n_sup = np.zeros([n_shuffles, replay_par_mat.shape[1]])

    # get the shuffled number of deep and superficial cells
    for i in range(replay_par_mat.shape[1]):
        n_deep[:, i], n_sup[:, i] = get_shuffled_labels(
            cell_metrics, replay_par_mat, i, n_shuffles=n_shuffles
        )

    # get the number of observed deep and superficial cells
    n_deep_obs = []
    n_sup_obs = []
    for i in range(replay_par_mat.shape[1]):
        cur_idx = replay_par_mat[:, i] == 1
        n_deep_obs.append(sum(cell_metrics[cur_idx].deepSuperficial == "Deep"))
        n_sup_obs.append(sum(cell_metrics[cur_idx].deepSuperficial == "Superficial"))

    # get the index of significant events
    sig_idx_deep, _ = functions.get_significant_events(
        np.hstack(n_deep_obs), n_deep, q=q_perc
    )
    sig_idx_sup, _ = functions.get_significant_events(
        np.hstack(n_sup_obs), n_sup, q=q_perc
    )

    return sig_idx_deep, sig_idx_sup


def run(
    basepath=None,
    replay_df=None,
    n_shuffles=1000,
    q_perc=90,
    expand_replay_epoch=0.05,
    brainRegion="CA1",
    putativeCellType="Pyr",
):
    """
    Run the analysis to locate replay events with significant bias
        of the number of deep and superficial cells
    Input:
        basepath: base path to the data
        replay_df: dataframe with the data from replay_run.load_results
        n_shuffles: number of shuffles for the label shuffling
        q_perc: percentile for the significance test
    Output:
        temp_df: dataframe with the results of the analysis

    """
    if basepath is None:
        raise ValueError("basepath is required")
    if replay_df is None:
        raise ValueError("df is required")

    df_idx = (replay_df.score_pval_time_swap < 0.05) & (replay_df.basepath == basepath)
    # get the replay epochs for this basepath
    start = replay_df[df_idx].start
    stop = replay_df[df_idx].stop
    replay_epochs = nel.EpochArray(np.array([start, stop]).T)
    # expand the epochs by xxms
    replay_epochs = replay_epochs.expand(expand_replay_epoch)
    # load the spikes
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion=brainRegion, putativeCellType=putativeCellType
    )
    # add deep superficial labels
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)
    # get the participation matrix
    replay_par_mat = functions.get_participation(
        st.data, replay_epochs.starts, replay_epochs.stops, par_type="binary"
    )
    # get the significant events
    sig_idx_deep, sig_idx_sup = get_significant_events(
        cell_metrics, replay_par_mat, n_shuffles=n_shuffles, q_perc=q_perc
    )
    # add the significant event identifiers to the dataframe
    temp_df = replay_df[df_idx]
    temp_df.reset_index(inplace=True)
    temp_df["sig_unit_bias"] = "unknown"
    temp_df.loc[sig_idx_sup, "sig_bias"] = "sup"
    temp_df.loc[sig_idx_deep, "sig_bias"] = "deep"

    return temp_df
