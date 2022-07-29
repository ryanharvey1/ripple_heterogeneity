import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
    compress_repeated_epochs,
)
import pandas as pd
import numpy as np
import nelpy as nel
import os
import itertools
from pyinform import conditional_entropy, mutual_info
import multiprocessing
from joblib import Parallel, delayed


def pairwise_info(X, return_index=True):
    """
    Compute the pairwise mutual information between all pairs of variables in X.
    inputs:
        X: a numpy array of shape (n,d) where d is the number of samples and n is the number of variables
    outputs:
        a numpy array of shape (n) where the (i,j) entry is mutual_info between the ith and jth variables
        pairs: a numpy array of shape (n,2) where the (i,j) entry is the pair of variables

    """
    x = np.arange(0, X.shape[0])
    pairs = np.array(list(itertools.combinations(x, 2)))
    mi = []
    for pair in pairs:
        mi.append(mutual_info(X[pair[0], :], X[pair[1], :]))
    if return_index:
        return np.hstack(mi), pairs
    else:
        return np.hstack(mi)


def pairwise_conditional_entropy(X):
    """
    Compute the pairwise mutual information between all pairs of variables in X.
    inputs:
        X: a numpy array of shape (n,d) where d is the number of samples and n is the number of variables
    outputs:
        a numpy array of shape (n) where the (i,j) entry is conditional_entropy between the ith and jth variables
        pairs: a numpy array of shape (n,2) where the (i,j) entry is the pair of variables

    """
    x = np.arange(0, X.shape[0])
    pairs = np.array(list(itertools.product(x, repeat=2)))
    ce = []
    for pair in pairs:
        ce.append(conditional_entropy(X[pair[0], :], X[pair[1], :]))
    return np.hstack(ce), pairs


def find_pairs(cm, pairs, pair_1, pair_2, layer):
    """
    Find the index for var pairs using in1d that match pair_1, pair_2, and layer
    inputs:
        cm: cell_metrics dataframe
        pairs: pairs index of variables
        pair_1: first pair of variables (e.g. "PFC", "EC1")
        pair_2: second pair of variables (e.g. "PFC", "EC1")
        layer: layer of variables (e.g. "Deep", "Superficial")
    outputs:
        index of pairs that match pair_1, pair_2, and layer
    """
    pairs_idx = (
        np.in1d(
            pairs[:, 0],
            np.where(
                cm.brainRegion.str.contains(pair_1)
                & cm.deepSuperficial.str.contains(layer)
            )[0],
        )
        & np.in1d(pairs[:, 1], np.where(cm.brainRegion.str.contains(pair_2))[0])
    ) | (
        np.in1d(pairs[:, 0], np.where(cm.brainRegion.str.contains(pair_2))[0])
        & np.in1d(
            pairs[:, 1],
            np.where(
                cm.brainRegion.str.contains(pair_1)
                & cm.deepSuperficial.str.contains(layer)
            )[0],
        )
    )
    return pairs_idx


def make_df(
    ce,
    pairs,
    cm,
    y_label=None,
    reference_region="CA1",
    target_region=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],
    optional_pval=None,
):
    """
    make a dataframe from the input ce and map to cell_metrics(cm)
    inputs:
        ce: metric between pairs of variables
        pairs: pairs index of variables
        cm: cell_metrics dataframe
        y_label: label for ce
        reference_region: reference region for ce
        target_region: target region for ce
    outputs:
        df: dataframe
    """
    deep_pfc = find_pairs(cm, pairs, reference_region, target_region[0], "Deep")
    deep_mec = find_pairs(cm, pairs, reference_region, target_region[1], "Deep")
    sup_pfc = find_pairs(cm, pairs, reference_region, target_region[0], "Superficial")
    sup_mec = find_pairs(cm, pairs, reference_region, target_region[1], "Superficial")

    df = pd.DataFrame()

    df[y_label] = np.hstack([ce[deep_pfc], ce[deep_mec], ce[sup_pfc], ce[sup_mec]])

    if optional_pval is not None:
        df["pval"] = np.hstack(
            [
                optional_pval[deep_pfc],
                optional_pval[deep_mec],
                optional_pval[sup_pfc],
                optional_pval[sup_mec],
            ]
        )

    df["label"] = np.hstack(
        [
            ["deep_pfc"] * len(ce[deep_pfc]),
            ["deep_mec"] * len(ce[deep_mec]),
            ["sup_pfc"] * len(ce[sup_pfc]),
            ["sup_mec"] * len(ce[sup_mec]),
        ]
    )
    df["reference_id"] = np.hstack(
        [
            pairs[deep_pfc, 0],
            pairs[deep_mec, 0],
            pairs[sup_pfc, 0],
            pairs[sup_mec, 0],
        ]
    )

    df["target_id"] = np.hstack(
        [
            pairs[deep_pfc, 1],
            pairs[deep_mec, 1],
            pairs[sup_pfc, 1],
            pairs[sup_mec, 1],
        ]
    )

    df["reference"] = np.hstack(
        [
            cm.brainRegion.iloc[pairs[deep_pfc, 0]],
            cm.brainRegion.iloc[pairs[deep_mec, 0]],
            cm.brainRegion.iloc[pairs[sup_pfc, 0]],
            cm.brainRegion.iloc[pairs[sup_mec, 0]],
        ]
    )

    df["target"] = np.hstack(
        [
            cm.brainRegion.iloc[pairs[deep_pfc, 1]],
            cm.brainRegion.iloc[pairs[deep_mec, 1]],
            cm.brainRegion.iloc[pairs[sup_pfc, 1]],
            cm.brainRegion.iloc[pairs[sup_mec, 1]],
        ]
    )
    df["reference_uid"] = np.hstack(
        [
            cm.UID.iloc[pairs[deep_pfc, 0]],
            cm.UID.iloc[pairs[deep_mec, 0]],
            cm.UID.iloc[pairs[sup_pfc, 0]],
            cm.UID.iloc[pairs[sup_mec, 0]],
        ]
    )

    df["target_uid"] = np.hstack(
        [
            cm.UID.iloc[pairs[deep_pfc, 1]],
            cm.UID.iloc[pairs[deep_mec, 1]],
            cm.UID.iloc[pairs[sup_pfc, 1]],
            cm.UID.iloc[pairs[sup_mec, 1]],
        ]
    )
    return df


def shuffle_mi(ripple_mat, n_shuffles=500, parallel=True):
    def shuffle_patterns(patterns):
        return np.random.permutation(patterns.flatten()).reshape(patterns.shape)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        shuffle_val = Parallel(n_jobs=num_cores)(
            delayed(pairwise_info)(shuffle_patterns(ripple_mat), return_index=False)
            for _ in range(n_shuffles)
        )
    else:
        shuffle_val = [
            pairwise_info(shuffle_patterns(ripple_mat), return_index=False)
            for _ in range(n_shuffles)
        ]
    return shuffle_val

def get_epochs(basepath):
    # load session epoch data
    epoch_df = loading.load_epoch(basepath)

    # compress repeated sleep sessions
    epoch_df = compress_repeated_epochs.main(epoch_df)
    # search for pre task post epochs
    idx, _ = functions.find_pre_task_post(epoch_df.environment)
    if idx is None:
        return None, None

    epoch_df = epoch_df[idx]
    epochs = nel.EpochArray(
        [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
        label=epoch_df.environment.values,
    )
    return epochs, epoch_df

def run(
    basepath,
    putativeCellType="Pyr",  # type of cell to use for the analysis
    reference_region="CA1",  # reference region
    target_region=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # downstream of reference_region
    brainRegions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to include
    rip_exp_start=0.05,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.2,  # ripple expansion stop, in seconds, how much to expand ripples
    parallel_shuffle=True,  # whether to run shuffle in parallel
    n_shuffles=500,  # number of shuffles to run
):

    epochs, epoch_df = get_epochs(basepath)
    if epochs is None:
        return None

    # load in spike data
    st, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegions
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

    if st.isempty:
        return None
    if st.n_active < 2:
        return None

    # load ripples and apply ripple expansion
    ripples_df = loading.load_ripples_events(basepath)
    ripples = (
        nel.EpochArray(np.array([ripples_df.start, ripples_df.stop]).T)
        .expand(rip_exp_start, direction="start")
        .expand(rip_exp_stop, direction="stop")
    )

    results = pd.DataFrame()
    # iterate over pre/task/post epochs
    for ep, env_label, ep_label in zip(
        epochs, epoch_df.environment.values, ["pre", "task", "post"]
    ):
        # bin spikes into ripples and get firing rate
        ripple_mat = functions.get_participation(
            st[ep].data, ripples[ep].starts, ripples[ep].stops, par_type="firing_rate"
        )

        # get mutual information
        mi, pairs_mi = pairwise_info(ripple_mat)
        # get conditional entropy
        ce, pairs_ce = pairwise_conditional_entropy(ripple_mat)

        # get shuffle values
        mi_shuff = shuffle_mi(ripple_mat, n_shuffles=n_shuffles, parallel=parallel_shuffle)
        _, pvalues = functions.get_significant_events(mi, np.array(mi_shuff))

        # add multual information and conditional entropy to dataframe
        mutual_info_df = make_df(
            mi,
            pairs_mi,
            cm,
            "mutual_info",
            reference_region,
            target_region,
            optional_pval=pvalues,
        )

        conditional_entropy_df = make_df(
            ce, pairs_ce, cm, "conditional_entropy", reference_region, target_region
        )
        # because conditional entropy is directional, remove the pairs where the reference is the target
        conditional_entropy_df.loc[
            ~conditional_entropy_df.reference.str.contains(reference_region),
            "conditional_entropy",
        ] = np.nan
        conditional_entropy_df.dropna(inplace=True)

        results_temp = pd.concat(
            [mutual_info_df, conditional_entropy_df], ignore_index=True
        )

        results_temp["environment"] = env_label
        results_temp["epoch"] = ep_label
        results_temp["basepath"] = basepath

        results = pd.concat([results, results_temp], ignore_index=True)

    return results


def load_results(save_path, verbose=False):
    """
    load_results: load results from a pickle file
    """

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    results = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results_ = pickle.load(f)
        if results_ is None:
            continue

        results = pd.concat([results, results_], ignore_index=True)

    return results
