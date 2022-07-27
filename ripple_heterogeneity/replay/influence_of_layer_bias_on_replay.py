import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel
import glob


def sample_weight_calculator(data, sample_weighting_method, beta=None):
    """
    Calculate sample weights for a given dataframe.
    Inputs:
        data: pandas dataframe of categorical data
        sample_weighting_method: string, one of "effective", "ins", "balance", "equal"
        beta: float, only used if sample_weighting_method is
            "effective" default: ((len(data) - 1) / len(data))
    Returns:
        sample_weights: numpy array of sample weights

    example:
        data = pd.DataFrame()
        data["classes"] = ["Deep", "Superficial", "Deep", "Deep", "Superficial", "Deep"]

        print(sample_weight_calculator(data, "effective"))
        print(sample_weight_calculator(data, "ins"))
        print(sample_weight_calculator(data, "isns"))
        print(sample_weight_calculator(data, "balance"))
        print(sample_weight_calculator(data, "equal"))

        >> [0.742 1.258]
        >> [0.667 1.333]
        >> [0.828 1.172]
        >> [0.75, 1.5]
        >> [1 1]

    https://www.kaggle.com/code/pnarerdoan/xgboost-imbalance-data/notebook

    """

    def get_weights_inverse_num_of_samples(data, power):
        samples_per_class = list(data.value_counts(sort=False))
        number_of_class = len(samples_per_class)
        weights_for_samples = 1.0 / np.array(np.power(list(samples_per_class), power))
        weights_for_samples = (
            weights_for_samples / np.sum(weights_for_samples) * number_of_class
        )
        return weights_for_samples

    def get_weights_effective_num_of_samples(data, beta):
        if beta is None:
            beta = (len(data) - 1) / len(data)

        samples_per_class = list(data.value_counts(sort=False))
        number_of_class = len(samples_per_class)

        effective_num = 1.0 - np.power(beta, list(samples_per_class))
        weights_for_samples = (1.0 - beta) / np.array(effective_num)
        weights_for_samples = (
            weights_for_samples / np.sum(weights_for_samples) * number_of_class
        )

        return weights_for_samples

    def get_balance_sample_weights(data):

        total_number_instace = len(data)
        samples_per_class = list(data.value_counts(sort=False))
        number_of_class = len(samples_per_class)

        weights_for_samples = []

        for i in range(len(samples_per_class)):
            weights_for_samples.append(
                total_number_instace / (number_of_class * list(samples_per_class)[i])
            )

        return weights_for_samples

    def get_equal_weights(data):
        samples_per_class = list(data.value_counts(sort=False))
        number_of_class = len(samples_per_class)
        weights_for_samples = np.ones((number_of_class,), dtype=int)

        return weights_for_samples

    if sample_weighting_method == "effective":
        weights_for_samples = get_weights_effective_num_of_samples(data, beta)
    elif sample_weighting_method == "ins":
        weights_for_samples = get_weights_inverse_num_of_samples(data, power=1)
    elif sample_weighting_method == "isns":
        weights_for_samples = get_weights_inverse_num_of_samples(data, power=0.5)
    elif sample_weighting_method == "balance":
        weights_for_samples = get_balance_sample_weights(data)
    elif sample_weighting_method == "equal":
        weights_for_samples = get_equal_weights(data)
    else:
        raise ValueError("sample_weighting_method not recognized")

    return weights_for_samples


def get_weighted_avg_pyr_dist(cell_metrics, replay_par_mat):
    """
    Calculate weighted average pyr distance
    Input:
        cell_metrics: cell metrics dataframe
        replay_par_mat: binary participation matrix
    Output:
        df: weighted average pyr distance
    """

    avg_pyr_dist = []
    weight_methods = []
    event_id = []
    weight_methods_options = ["equal", "balance", "ins", "isns", "effective"]

    for i in range(replay_par_mat.shape[1]):

        # pull out index of active cells
        cur_idx = replay_par_mat[:, i] == 1

        # get active cells
        temp_df = cell_metrics[cur_idx]

        # keep only deep and superficial cells
        temp_df = temp_df[
            (temp_df.deepSuperficial == "Deep")
            | (temp_df.deepSuperficial == "Superficial")
        ]

        # initiate weights
        weights = np.tile(1.0, len(temp_df))

        # make df of class labels
        data = pd.DataFrame()
        data["classes"] = temp_df.deepSuperficial.values

        # if no data, skip, but add nans
        if len(data) == 0:
            for method in weight_methods_options:
                avg_pyr_dist.append(np.nan)
                weight_methods.append(method)
                event_id.append(i)
            continue

        # iterate through weight methods
        for method in weight_methods_options:

            # add weights to respective locations
            for weight, classes in zip(
                sample_weight_calculator(data, method),
                data.value_counts(sort=False).reset_index().classes.values,
            ):
                weights[np.where(temp_df.deepSuperficial.values == classes)[0]] = weight

            if len(weights) == 0:
                weight_methods.append(np.nan)
                event_id.append(i)
                continue

            # calculate average distance
            avg_pyr_dist.append(
                np.average(
                    temp_df.deepSuperficialDistance.values,
                    weights=weights,
                )
            )
            # store weight method
            weight_methods.append(method)
            event_id.append(i)

    df = pd.DataFrame()
    df["avg_pyr_dist"] = np.hstack(avg_pyr_dist)
    df["weight_methods"] = np.hstack(weight_methods)
    df["event_id"] = np.hstack(event_id)

    # make df long to wide format
    df = pd.pivot(
        df, index="event_id", columns="weight_methods", values="avg_pyr_dist"
    ).reset_index()

    return df


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
        n_middle: shuffled number of middle cells
    """
    n_deep = []
    n_sup = []
    n_middle = []
    for _ in range(n_shuffles):
        labels = shuffle_labels(cell_metrics.deepSuperficial)
        n_deep.append(sum(labels[replay_par_mat[:, ep] == 1] == "Deep"))
        n_sup.append(sum(labels[replay_par_mat[:, ep] == 1] == "Superficial"))
        n_middle.append(sum(labels[replay_par_mat[:, ep] == 1] == "middle"))
    return n_deep, n_sup, n_middle


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
    n_middle = np.zeros([n_shuffles, replay_par_mat.shape[1]])

    # get the shuffled number of deep and superficial cells
    for i in range(replay_par_mat.shape[1]):
        n_deep[:, i], n_sup[:, i], n_middle[:, i] = get_shuffled_labels(
            cell_metrics, replay_par_mat, i, n_shuffles=n_shuffles
        )

    # get the number of observed deep and superficial cells
    n_deep_obs = []
    n_sup_obs = []
    n_middle_obs = []
    for i in range(replay_par_mat.shape[1]):
        cur_idx = replay_par_mat[:, i] == 1
        deep_idx = cell_metrics[cur_idx].deepSuperficial.values == "Deep"
        sup_idx = cell_metrics[cur_idx].deepSuperficial.values == "Superficial"
        middle_idx = cell_metrics[cur_idx].deepSuperficial.values == "middle"

        n_deep_obs.append(sum(deep_idx))
        n_sup_obs.append(sum(sup_idx))
        n_middle_obs.append(sum(middle_idx))

    # get the index of significant events
    sig_idx_deep, pval_deep = functions.get_significant_events(
        np.hstack(n_deep_obs), n_deep, q=q_perc
    )
    sig_idx_sup, pval_sup = functions.get_significant_events(
        np.hstack(n_sup_obs), n_sup, q=q_perc
    )
    sig_idx_middle, pval_middle = functions.get_significant_events(
        np.hstack(n_middle_obs), n_middle, q=q_perc
    )

    return (
        sig_idx_deep,
        sig_idx_sup,
        sig_idx_middle,
        n_deep_obs,
        n_sup_obs,
        n_middle_obs,
        pval_deep,
        pval_sup,
        pval_middle,
    )


def run(
    basepath=None,
    replay_df=None, # provide if you want to just look at replay events
    use_replay_df=False, # if true, use replay_df instead of all ripples
    n_shuffles=500, # number of shuffles
    q_perc=95, # percentile of significant events threshold
    expand_replay_epoch=0.05, # how much to expand the replay/ripple epoch
    brainRegion="CA1", # brain region
    putativeCellType="Pyr", # cell type
    only_significant_replay=False, # if true, only look at significant replay events
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

    if use_replay_df:
        if only_significant_replay:
            df_idx = (replay_df.score_pval_time_swap < 0.05) & (
                replay_df.basepath == basepath
            )
            replay_df = replay_df[df_idx]
        else:
            df_idx = replay_df.basepath == basepath
            replay_df = replay_df[df_idx]
    else:
        replay_df = loading.load_ripples_events(basepath)


    # get the replay epochs for this basepath
    replay_epochs = nel.EpochArray(np.array([replay_df.start, replay_df.stop]).T)
    # expand the epochs by xxms
    replay_epochs = replay_epochs.expand(expand_replay_epoch)

    # load the spikes
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion=brainRegion, putativeCellType=putativeCellType
    )
    
    if st.isempty:
        return None

    # add deep superficial labels
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)

    # get the participation matrix
    replay_par_mat = functions.get_participation(
        st.data, replay_epochs.starts, replay_epochs.stops, par_type="binary"
    )

    # get the significant events
    (
        sig_idx_deep,
        sig_idx_sup,
        sig_idx_middle,
        n_deep_obs,
        n_sup_obs,
        n_middle_obs,
        pval_deep,
        pval_sup,
        pval_middle,
    ) = get_significant_events(
        cell_metrics, replay_par_mat, n_shuffles=n_shuffles, q_perc=q_perc
    )

    weighted_df = get_weighted_avg_pyr_dist(cell_metrics, replay_par_mat)

    # add the significant event identifiers to the dataframe
    temp_df = replay_df
    temp_df.reset_index(drop=True, inplace=True)
    temp_df["sig_unit_bias"] = "unknown"
    temp_df.loc[sig_idx_sup, "sig_unit_bias"] = "sup"
    temp_df.loc[sig_idx_deep, "sig_unit_bias"] = "deep"
    temp_df.loc[sig_idx_middle, "sig_unit_bias"] = "middle"
    temp_df["n_deep_obs"] = n_deep_obs
    temp_df["n_sup_obs"] = n_sup_obs
    temp_df["n_middle_obs"] = n_middle_obs
    temp_df["pval_deep"] = pval_deep
    temp_df["pval_sup"] = pval_sup
    temp_df["pval_middle"] = pval_middle

    temp_df = pd.concat([temp_df, weighted_df], axis=1)

    return temp_df


def load_results(save_path):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    df = pd.DataFrame()
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        df = pd.concat([df, results], ignore_index=True)
    return df
