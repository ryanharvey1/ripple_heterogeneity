import glob
import os
import pickle
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score
from sklearn.tree import ExtraTreeClassifier
from ripple_heterogeneity.utils import (
    loading,
    add_new_deep_sup,
)
import nelpy as nel
import pandas as pd
import numpy as np


def get_data(basepath, putativeCellType="Pyr", brainRegion="CA1", min_epoch_dur=5):
    # get behavioral epochs
    epoch_df = loading.load_epoch(basepath)
    # remove sleep epochs
    epoch_df = epoch_df[epoch_df.environment != "sleep"]
    # remove epochs that were too short
    epoch_df = epoch_df[((epoch_df.stopTime - epoch_df.startTime) / 60) > min_epoch_dur]
    epochs = nel.EpochArray(
        [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
        label=epoch_df.environment.values,
    )
    # load spikes
    spikes, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegion
    )

    return epochs, spikes, cm


def run(
    basepath,
    bin_size=30,
    smooth_sigma=15,
    putativeCellType="Pyr",
    brainRegion="CA1",
    min_epoch_dur=5,
):

    epochs, spikes, cm = get_data(
        basepath,
        putativeCellType=putativeCellType,
        brainRegion=brainRegion,
        min_epoch_dur=min_epoch_dur,
    )

    if epochs.n_intervals < 2:
        return None

    # bin spikes and smooth
    bst = spikes[epochs].bin(ds=bin_size).smooth(sigma=smooth_sigma)
    # construct y for classification (vector of numeric labels for each epoch)
    y = np.hstack([np.zeros(bst_.n_bins) + i_epoch for i_epoch, bst_ in enumerate(bst)])

    # pull out the data for the training and test set
    X = bst.data.T
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.33, random_state=42
    )

    # construct the pipeline for the classifier
    clf = make_pipeline(StandardScaler(), ExtraTreeClassifier())

    scores = []
    bal_acc = []

    # iterate over each cell to see how well it can predict the environment
    for i in range(X_train.shape[1]):
        # train the classifier
        clf.fit(X_train[:, i].reshape(-1, 1), y_train)
        # calculate accuracy of the classifier
        scores.append(clf.score(X_test[:, i].reshape(-1, 1), y_test))
        bal_acc.append(
            balanced_accuracy_score(y_test, clf.predict(X_test[:, i].reshape(-1, 1)))
        )

    # store the results in a dataframe
    results = pd.DataFrame(
        {
            "UID": cm.UID,
            "accuracy": scores,
            "bal_accuracy": bal_acc,
        }
    )
    results["chance_accuracy"] = 1 / epochs.n_intervals
    results["n_intervals"] = epochs.n_intervals
    results["deepSuperficial"] = cm.deepSuperficial
    results["deepSuperficialDistance"] = cm.deepSuperficialDistance
    results["brainRegion"] = cm.brainRegion
    results["putativeCellType"] = cm.putativeCellType
    results["basepath"] = basepath
    results = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(results)

    return results

def load_results(save_path):
    """
    load_results: load results from a pickle file
    """

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    results_df = pd.DataFrame()
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        results_df = pd.concat([results_df, results], ignore_index=True)