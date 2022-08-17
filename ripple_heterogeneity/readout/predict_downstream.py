import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel
from ripple_heterogeneity.utils import compress_repeated_epochs
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import (
    mean_squared_error,
    mean_squared_log_error,
    median_absolute_error,
)
from sklearn.decomposition import PCA
from sklearn import preprocessing


def predict(X, y):
    """
    predict: predict the downstream activity of the target regions
    """
    # remove nan and inf
    bad_idx = np.hstack(np.isinf(y) | np.isnan(y))
    y = y[~bad_idx]
    X = X[~bad_idx, :]

    scaler = preprocessing.StandardScaler()
    X = scaler.fit_transform(X)

    # split into train and test
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.4, random_state=0
    )
    # train a linear regression model
    return RidgeCV().fit(X_train, y_train), X_train, X_test, y_train, y_test


def shuffle_data(X, y, n_shuff=1000):
    mse = []
    medse = []
    for i in range(n_shuff):
        idx = np.random.permutation(len(X))
        reg, X_train, X_test, y_train, y_test = predict(X[idx], y)
        # get model performance
        mse.append(mean_squared_error(y_test, reg.predict(X_test)))
        medse.append(median_absolute_error(y_test, reg.predict(X_test)))
    return mse, medse


def run(
    basepath,  # path to data folder
    reference_region=["CA1"],  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    ripple_expand=0.1,  # in seconds, how much to expand ripples
    ev_thres=0.8,  # explained variance threshold for PCA
    min_ripples=10,  # minimum number of ripples per epoch
    n_shuff=1000,  # number of shuffles to do
    par_type="binary"
):

    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion=[*target_regions, *reference_region]
    )
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)

    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T]).expand(
        ripple_expand
    )

    ep_df = loading.load_epoch(basepath)
    ep_df = compress_repeated_epochs.main(ep_df, epoch_name="sleep")
    # locate pre task post structure
    idx, _ = functions.find_pre_task_post(ep_df.environment)
    if idx is None:
        return None

    ep_df = ep_df[idx]
    ep_epochs = nel.EpochArray([np.array([ep_df.startTime, ep_df.stopTime]).T])

    ca1_deep_idx = (
        cell_metrics.brainRegion.str.contains("CA1").values
        & (cell_metrics.deepSuperficial == "Deep")
        & (cell_metrics.putativeCellType.str.contains("Pyr"))
    )
    ca1_sup_idx = (
        cell_metrics.brainRegion.str.contains("CA1").values
        & (cell_metrics.deepSuperficial == "Superficial")
        & (cell_metrics.putativeCellType.str.contains("Pyr"))
    )
    if (sum(ca1_deep_idx) < min_cells) | (sum(ca1_sup_idx) < min_cells):
        return None

    epoch = []
    epoch_i = []
    targ_reg = []
    n_x_components = []
    mse = []
    mean_score = []
    std_score = []
    n_deep = []
    n_sup = []
    n_target_cells = []
    test_score = []
    train_score = []
    medse = []
    median_medse_shuff = []
    median_mse_shuff = []
    for ep_i, ep in enumerate(ep_epochs):
        # get participation for every cell
        st_par = functions.get_participation(
            st[ep].data,
            ripple_epochs[ep].starts,
            ripple_epochs[ep].stops,
            par_type=par_type,
        )
        if st_par.shape[1] < min_ripples:
            continue
        # calculate ratio of n deep and n superficial cells
        deep_sup_ratio = (st_par[ca1_deep_idx, :] > 0).sum(axis=0) / (
            st_par[ca1_sup_idx, :] > 0
        ).sum(axis=0)

        for region in target_regions:
            if sum(cell_metrics.brainRegion.str.contains(region).values) < min_cells:
                continue
            # log transform to get better predictions
            y = np.log(deep_sup_ratio + 1).reshape(-1, 1)
            # get target participation data
            X = st_par[cell_metrics.brainRegion.str.contains(region).values, :]
            # # get pca dims that explain XX of the variance
            X = PCA(n_components=ev_thres, svd_solver="full").fit_transform(X.T)

            reg, X_train, X_test, y_train, y_test = predict(X, y)

            # # get the predicted values
            pred = reg.predict(X_test)

            # get model performance
            mse.append(mean_squared_error(y_test, pred))
            medse.append(median_absolute_error(y_test, pred))
            scores = cross_val_score(reg, X, y, cv=5)
            mean_score.append(scores.mean())
            std_score.append(scores.std())
            test_score.append(reg.score(X_test, y_test))
            train_score.append(reg.score(X_train, y_train))
            # get metadata
            n_x_components.append(X.shape[1])
            epoch.append(ep_df.environment.iloc[ep_i])
            epoch_i.append(ep_i)
            targ_reg.append(region)
            n_deep.append(sum(ca1_deep_idx))
            n_sup.append(sum(ca1_sup_idx))
            n_target_cells.append(
                sum(cell_metrics.brainRegion.str.contains(region).values)
            )
            # get vars for prediction gain
            mse_shuff, medse_shuff = shuffle_data(X, y, n_shuff=n_shuff)
            median_medse_shuff.append(np.median(medse_shuff))
            median_mse_shuff.append(np.median(mse_shuff))

    if len(epoch) == 0:
        return pd.DataFrame()

    # create a dataframe
    df = pd.DataFrame()
    df["epoch"] = np.hstack(epoch)
    df["epoch_i"] = np.hstack(epoch_i)
    df["targ_reg"] = np.hstack(targ_reg)
    df["n_x_components"] = np.hstack(n_x_components)
    df["mse"] = np.hstack(mse)
    df["medse"] = np.hstack(medse)
    df["median_mse_shuff"] = np.hstack(median_mse_shuff)
    df["median_medse_shuff"] = np.hstack(median_medse_shuff)
    # df["rmsle"] = np.hstack(rmsle)
    df["mean_score"] = np.hstack(mean_score)
    df["std_score"] = np.hstack(std_score)
    df["test_score"] = np.hstack(test_score)
    df["train_score"] = np.hstack(train_score)
    df["n_deep"] = np.hstack(n_deep)
    df["n_sup"] = np.hstack(n_sup)
    df["n_target_cells"] = np.hstack(n_target_cells)
    df["basepath"] = basepath

    return df


def load_results(save_path, verbose=False):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    df = pd.DataFrame()
    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        df = pd.concat([df, results], ignore_index=True)

    return df
