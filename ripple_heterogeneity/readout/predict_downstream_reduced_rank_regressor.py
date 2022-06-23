import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import (
    functions,
    compress_repeated_epochs,
    loading,
    add_new_deep_sup,
    reduced_rank_regressor,
    kernel_reduced_rank_ridge_regression,
)
import nelpy as nel
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from scipy import around
from scipy import size
from scipy.linalg import norm
from sklearn.cross_decomposition import CCA, PLSCanonical, PLSRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error


def sqerr(matrix1, matrix2):
    """Squared error (frobenius norm of diff) between two matrices."""
    return around(pow(norm(matrix1 - matrix2, "fro"), 2) / size(matrix2, 0), 5)


def shuffle_data(X, y, rank, reg, n_shuff=1000):
    testing_error = []
    r2_test = []
    for i in range(n_shuff):
        idx = np.random.permutation(X.shape[0])
        X_train, X_test, y_train, y_test = train_test_split(
            X[idx, :], y, test_size=0.4, random_state=42
        )
        regressor = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()
        regressor.rank = int(rank)
        regressor.reg = reg
        regressor.fit(X_train, y_train)

        # get model performance
        testing_error.append(sqerr(regressor.predict(X_test), y_test))
        r2_test.append(regressor.score(X_test, y_test))

    return testing_error, r2_test


def get_data(basepath, target_regions, reference_region, ripple_expand):
    st, cm = loading.load_spikes(
        basepath, brainRegion=[*target_regions, *reference_region]
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T]).expand(
        ripple_expand
    )

    ep_df = loading.load_epoch(basepath)
    ep_df = compress_repeated_epochs.main(ep_df, epoch_name="sleep")
    session_epoch = nel.EpochArray(
        [np.array([ep_df.startTime.iloc[0], ep_df.stopTime.iloc[-1]]).T]
    )

    # locate pre task post structure
    idx, _ = functions.find_pre_task_post(ep_df.environment)
    if idx is None:
        return None, None, None, None, None

    ep_df = ep_df[idx]
    ep_epochs = nel.EpochArray([np.array([ep_df.startTime, ep_df.stopTime]).T])
    return st, cm, ripple_epochs, ep_epochs, ep_df, session_epoch


def run_grid_search(X_train, y_train, n_grid=10, cv=5):
    """
    grid_search: grid search for the reduced rank regressor
    """
    rank_grid = np.linspace(
        1, min(X_train.shape[1], y_train.shape[1]), num=n_grid
    ).astype(int)

    reg_grid = np.power(10, np.linspace(-20, 20, num=n_grid + 1))

    parameters_grid_search = {"reg": reg_grid, "rank": rank_grid}

    rrr = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()

    grid_search = GridSearchCV(
        rrr,
        parameters_grid_search,
        cv=cv,
        scoring="neg_mean_squared_error",
        n_jobs=-1,
    )
    return grid_search.fit(X_train, y_train)


def run(
    basepath,  # path to data folder
    reference_region=["CA1"],  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    ripple_expand=0.2,  # in seconds, how much to expand ripples
    min_ripples=10,  # minimum number of ripples per epoch
    n_shuff=1000,  # number of shuffles to do
    rank=10,  # rank of the reduced rank regressor (not used)
    reg=1e-6,  # regularization parameter (not used)
    source_cell_type="Pyr",  # source cell type
    target_cell_type=None,  # cell type to use for target cells
    n_grid=10,  # number of grid search parameters to use
    cv=5,  # number of cross validation folds
    use_entire_session=False,  # use entire session or just pre task post
):

    st, cm, ripple_epochs, ep_epochs, ep_df, session_epoch = get_data(
        basepath, target_regions, reference_region, ripple_expand
    )
    if st is None:
        return None

    epoch = []
    epoch_i = []
    targ_reg = []
    n_x_components = []
    n_target_cells = []
    training_error = []
    testing_error = []
    median_error_shuff = []
    mean_error_shuff = []
    n_ca1 = []
    ca1_sub_layer = []
    r2_rrr_train = []
    r2_rrr_test = []
    mean_r2_shuff = []
    median_r2_shuff = []
    r2_shuffles = []
    r2_cca = []
    r2_plsr = []
    r2_plsc = []
    rrr_rank = []
    rrr_reg = []
    mse_cca = []
    mse_plsc = []
    mse_plsr = []
    scaler = preprocessing.StandardScaler()

    if use_entire_session:
        ep_epochs = session_epoch

    # iterate over all epochs
    for ep_i, ep in enumerate(ep_epochs):
        # continue if there are too few ripples
        if len(ripple_epochs[ep].starts) < min_ripples:
            continue

        # get participation for every cell
        # this will continue if there is insufficient data
        try:
            st_par = functions.get_participation(
                st[ep].data,
                ripple_epochs[ep].starts,
                ripple_epochs[ep].stops,
                par_type="firing_rate",
            )
        except:
            continue

        # rescale using standard scaler
        X = scaler.fit_transform(st_par)

        # iterate over ca1 sublayers regions
        for ca1_sub in ["Deep", "Superficial"]:
            # iterate over target regions
            for region in target_regions:
                if sum(cm.brainRegion.str.contains(region).values) < min_cells:
                    continue

                ca1_idx = (
                    cm.brainRegion.str.contains("CA1").values
                    & (cm.deepSuperficial == ca1_sub)
                    & (cm.putativeCellType.str.contains(source_cell_type))
                )
                if sum(ca1_idx) < min_cells:
                    continue

                if target_cell_type is not None:
                    target_idx = (
                        cm.brainRegion.str.contains(region).values
                        & cm.putativeCellType.str.contains(target_cell_type).values
                    )
                else:
                    target_idx = cm.brainRegion.str.contains(region).values

                X_train, X_test, y_train, y_test = train_test_split(
                    X[ca1_idx, :].T,
                    X[target_idx, :].T,
                    test_size=0.4,
                    random_state=42,
                )
                grid_search = run_grid_search(X_train, y_train, n_grid=n_grid, cv=cv)

                regressor = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()
                regressor.rank = int(grid_search.best_params_["rank"])
                regressor.reg = grid_search.best_params_["reg"]
                regressor.fit(X_train, y_train)

                mdl = CCA().fit(X_train, y_train)
                r2_cca.append(mdl.score(X_test, y_test))
                mse_cca.append(mean_squared_error(y_test, mdl.predict(X_test)))

                mdl = PLSCanonical().fit(X_train, y_train)
                r2_plsc.append(mdl.score(X_test, y_test))
                mse_plsc.append(mean_squared_error(y_test, mdl.predict(X_test)))

                mdl = PLSRegression().fit(X_train, y_train)
                r2_plsr.append(mdl.score(X_test, y_test))
                mse_plsr.append(mean_squared_error(y_test, mdl.predict(X_test)))

                # get model performance
                training_error.append(
                    mean_squared_error(y_train, regressor.predict(X_train))
                )
                testing_error.append(
                    mean_squared_error(y_test, regressor.predict(X_test))
                )
                r2_rrr_train.append(regressor.score(X_train, y_train))
                r2_rrr_test.append(regressor.score(X_test, y_test))

                # get metadata
                rrr_rank.append(regressor.rank)
                rrr_reg.append(regressor.reg)
                n_x_components.append(X.shape[1])
                epoch.append(ep_df.environment.iloc[ep_i])
                epoch_i.append(ep_i)
                targ_reg.append(region)
                ca1_sub_layer.append(ca1_sub)
                n_ca1.append(sum(ca1_idx))
                n_target_cells.append(sum(cm.brainRegion.str.contains(region).values))

                # get vars for shuffles
                if n_shuff > 0:
                    error_shuff, r2_shuff = shuffle_data(
                        X[ca1_idx, :].T,
                        X[target_idx, :].T,
                        regressor.rank,
                        regressor.reg,
                        n_shuff=n_shuff,
                    )
                    median_error_shuff.append(np.median(error_shuff))
                    mean_error_shuff.append(np.mean(error_shuff))
                    mean_r2_shuff.append(np.mean(r2_shuff))
                    median_r2_shuff.append(np.median(r2_shuff))
                    r2_shuffles.append(r2_shuff)

    if len(epoch) == 0:
        return pd.DataFrame()

    # create a dataframe
    df = pd.DataFrame()
    df["epoch"] = np.hstack(epoch)
    df["epoch_i"] = np.hstack(epoch_i)
    df["targ_reg"] = np.hstack(targ_reg)
    df["ca1_sub_layer"] = np.hstack(ca1_sub_layer)
    df["n_x_components"] = np.hstack(n_x_components)
    df["training_error"] = np.hstack(training_error)
    df["testing_error"] = np.hstack(testing_error)
    df["mse_cca"] = np.hstack(mse_cca)
    df["mse_plsc"] = np.hstack(mse_plsc)
    df["mse_plsr"] = np.hstack(mse_plsr)
    df["r2_rrr_train"] = np.hstack(r2_rrr_train)
    df["r2_rrr_test"] = np.hstack(r2_rrr_test)
    df["rrr_rank"] = np.hstack(rrr_rank)
    df["rrr_reg"] = np.hstack(rrr_reg)
    df["r2_cca"] = np.hstack(r2_cca)
    df["r2_plsc"] = np.hstack(r2_plsc)
    df["r2_plsr"] = np.hstack(r2_plsr)
    if n_shuff > 0:
        df["mean_error_shuff"] = np.hstack(mean_error_shuff)
        df["median_error_shuff"] = np.hstack(median_error_shuff)
        df["mean_r2_shuff"] = np.hstack(mean_r2_shuff)
        df["median_r2_shuff"] = np.hstack(median_r2_shuff)
        _, pval = functions.get_significant_events(
            np.hstack(r2_rrr_train), np.vstack(r2_shuffles)
        )
        df["pvalues"] = pval
    df["n_ca1"] = np.hstack(n_ca1)
    df["n_target_cells"] = np.hstack(n_target_cells)
    df["basepath"] = basepath
    df["use_entire_session"] = use_entire_session


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
