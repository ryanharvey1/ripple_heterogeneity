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
from sklearn.model_selection import cross_validate
from sklearn.model_selection import TimeSeriesSplit


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
        r2_test.append(regressor.score(X_test, y_test))
        testing_error.append(mean_squared_error(y_test, regressor.predict(X_test)))
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
        return None, None, None, None, None, None, None, None

    ep_df = ep_df[idx]
    ep_epochs = nel.EpochArray([np.array([ep_df.startTime, ep_df.stopTime]).T])

    state_dict = loading.load_SleepState_states(basepath)
    nrem_epochs = nel.EpochArray(state_dict["NREMstate"])
    wake_epochs = nel.EpochArray(state_dict["WAKEstate"])

    return (
        st,
        cm,
        ripple_epochs,
        ep_epochs,
        ep_df,
        session_epoch,
        nrem_epochs,
        wake_epochs,
    )


def run_grid_search(X_train, y_train, n_grid=10, cv=5, max_rank=64):
    """
    grid_search: grid search for the reduced rank regressor
    """

    rank_grid = np.arange(1, min(X_train.shape[1], y_train.shape[1], max_rank)).astype(
        int
    )

    parameters_grid_search = {"rank": rank_grid}

    rrr = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()

    # folds = TimeSeriesSplit(n_splits=cv)

    grid_search = GridSearchCV(
        rrr,
        parameters_grid_search,
        cv=cv,
        scoring="neg_mean_squared_error",
        n_jobs=-2,
    )
    return grid_search.fit(X_train, y_train)


# def evaluate(model, X, y, cv, verbose=False):
#     cv_results = cross_validate(
#         model,
#         X,
#         y,
#         cv=cv,
#         scoring=["neg_mean_absolute_error", "neg_root_mean_squared_error", "r2"],
#     )
#     mae = -cv_results["test_neg_mean_absolute_error"]
#     rmse = -cv_results["test_neg_root_mean_squared_error"]
#     r2 = cv_results["test_r2"]
#     if verbose:

#         print(
#             f"Mean Absolute Error:     {mae.mean():.3f} +/- {mae.std():.3f}\n"
#             f"Root Mean Squared Error: {rmse.mean():.3f} +/- {rmse.std():.3f}\n"
#             f"R2:                      {r2.mean():.3f} +/- {r2.std():.3f}"
#         )
#     return mae, rmse, r2


def run(
    basepath,  # path to data folder
    reference_region=["CA1"],  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    ripple_expand=0.2,  # in seconds, how much to expand ripples
    min_ripples=10,  # minimum number of ripples per epoch
    n_shuff=1000,  # number of shuffles to do
    rank=10,  # rank of the reduced rank regressor (not used)
    reg=1,  # regularization parameter
    source_cell_type="Pyr",  # source cell type
    target_cell_type=None,  # cell type to use for target cells
    n_grid=20,  # number of grid search parameters to use
    cv=5,  # number of cross validation folds
    max_rank=64,  # maximum rank to use in the reduced rank regressor
    use_entire_session=False,  # use entire session or just pre task post
):

    (
        st,
        cm,
        ripple_epochs,
        ep_epochs,
        ep_df,
        session_epoch,
        nrem_epochs,
        wake_epochs,
    ) = get_data(basepath, target_regions, reference_region, ripple_expand)
    if st is None:
        return None

    # initialize output vars
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
    states = []
    train_error_units = []
    test_error_units = []
    rrr_rank_units = []
    rrr_reg_units = []
    n_x_components_units = []
    epoch_units = []
    epoch_i_units = []
    states_units = []
    targ_reg_units = []
    ca1_sub_layer_units = []
    n_ca1_units = []
    n_target_cells_units = []
    mse_cca_units = []
    mse_plsc_units = []
    mse_plsr_units = []
    target_uid = []

    scaler = preprocessing.StandardScaler()
    # ts_cv = TimeSeriesSplit(n_splits=cv)

    # if use_entire_session, use the entire session, otherwise use pre task post
    if use_entire_session:
        ep_epochs = session_epoch

    states_ = ["NREM", "WAKE"]

    # iterate over all epochs
    for ep_i, ep in enumerate(ep_epochs):

        for state_i, state in enumerate([nrem_epochs, wake_epochs]):
            # get the ripple epochs for this state
            curr_ripples = ripple_epochs[ep][state]

            # continue if there are too few ripples
            if len(curr_ripples.starts) < min_ripples:
                continue

            # get participation for every cell
            # this will continue if there is insufficient data
            try:
                st_par = functions.get_participation(
                    st[ep][state].data,
                    curr_ripples.starts,
                    curr_ripples.stops,
                    par_type="firing_rate",
                )
                # start and end times can be the same for some reason, which
                # causes nan as the denominator above is 0
                st_par[np.isnan(st_par)] = 0
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

                    # get index of ca1 cells
                    ca1_idx = (
                        cm.brainRegion.str.contains("CA1").values
                        & (cm.deepSuperficial == ca1_sub)
                        & (cm.putativeCellType.str.contains(source_cell_type))
                    )

                    # get index of target cells
                    if target_cell_type is not None:
                        target_idx = (
                            cm.brainRegion.str.contains(region).values
                            & cm.putativeCellType.str.contains(target_cell_type).values
                        )
                    else:
                        target_idx = cm.brainRegion.str.contains(region).values

                    if (sum(ca1_idx) < min_cells) | (sum(target_idx) < min_cells):
                        continue

                    X_train, X_test, y_train, y_test = train_test_split(
                        X[ca1_idx, :].T,
                        X[target_idx, :].T,
                        test_size=0.4,
                        random_state=42,
                        shuffle=False,
                    )
                    grid_search = run_grid_search(
                        X_train, y_train, n_grid=n_grid, cv=cv, max_rank=max_rank
                    )

                    regressor = (
                        kernel_reduced_rank_ridge_regression.ReducedRankRegressor()
                    )
                    regressor.rank = int(grid_search.best_params_["rank"])
                    regressor.reg = reg

                    # evaluate(regressor, X, y, cv, verbose=False)

                    regressor.fit(X_train, y_train)

                    mdl_cca = CCA().fit(X_train, y_train)
                    r2_cca.append(mdl_cca.score(X_test, y_test))
                    mse_cca.append(mean_squared_error(y_test, mdl_cca.predict(X_test)))

                    mdl_plsc = PLSCanonical().fit(X_train, y_train)
                    r2_plsc.append(mdl_plsc.score(X_test, y_test))
                    mse_plsc.append(mean_squared_error(y_test, mdl_plsc.predict(X_test)))

                    mdl_plsr = PLSRegression().fit(X_train, y_train)
                    r2_plsr.append(mdl_plsr.score(X_test, y_test))
                    mse_plsr.append(mean_squared_error(y_test, mdl_plsr.predict(X_test)))

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
                    states.append(states_[state_i])
                    targ_reg.append(region)
                    ca1_sub_layer.append(ca1_sub)
                    n_ca1.append(sum(ca1_idx))
                    n_target_cells.append(
                        sum(cm.brainRegion.str.contains(region).values)
                    )

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

                    # multi-cell_performance
                    train_error_units.append(
                        mean_squared_error(
                            y_train,
                            regressor.predict(X_train),
                            multioutput="raw_values",
                        )
                    )
                    test_error_units.append(
                        mean_squared_error(
                            y_test, regressor.predict(X_test), multioutput="raw_values"
                        )
                    )

                    mse_cca_units.append(
                        mean_squared_error(
                            y_test, mdl_cca.predict(X_test), multioutput="raw_values"
                        )
                    )
                    mse_plsc_units.append(
                        mean_squared_error(
                            y_test, mdl_plsc.predict(X_test), multioutput="raw_values"
                        )
                    )
                    mse_plsr_units.append(
                        mean_squared_error(
                            y_test, mdl_plsr.predict(X_test), multioutput="raw_values"
                        )
                    )

                    # get metadata
                    n_cells = len(test_error_units[-1])
                    rrr_rank_units.append(np.tile(regressor.rank,n_cells))
                    rrr_reg_units.append(np.tile(regressor.reg,n_cells))
                    n_x_components_units.append(np.tile(X.shape[1],n_cells))
                    epoch_units.append(np.tile(ep_df.environment.iloc[ep_i],n_cells))
                    epoch_i_units.append(np.tile(ep_i,n_cells))
                    states_units.append(np.tile(states_[state_i],n_cells))
                    targ_reg_units.append(np.tile(region,n_cells))
                    ca1_sub_layer_units.append(np.tile(ca1_sub,n_cells))
                    n_ca1_units.append(np.tile(sum(ca1_idx),n_cells))
                    n_target_cells_units.append(
                        np.tile(sum(cm.brainRegion.str.contains(region).values),n_cells)
                    )
                    target_uid.append(cm[target_idx].UID.values)

    if len(epoch) == 0:
        return pd.DataFrame()

    # create a dataframe
    df = pd.DataFrame()
    df["epoch"] = np.hstack(epoch)
    df["epoch_i"] = np.hstack(epoch_i)
    df["state"] = np.hstack(states)
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

    df_unit = pd.DataFrame()
    df_unit["epoch"] = np.hstack(epoch_units)
    df_unit["epoch_i"] = np.hstack(epoch_i_units)
    df_unit["state"] = np.hstack(states_units)
    df_unit["targ_reg"] = np.hstack(targ_reg_units)
    df_unit["ca1_sub_layer"] = np.hstack(ca1_sub_layer_units)
    df_unit["n_x_components"] = np.hstack(n_x_components_units)
    df_unit["training_error"] = np.hstack(train_error_units)
    df_unit["testing_error"] = np.hstack(test_error_units)
    df_unit["mse_cca"] = np.hstack(mse_cca_units)
    df_unit["mse_plsc"] = np.hstack(mse_plsc_units)
    df_unit["mse_plsr"] = np.hstack(mse_plsr_units)
    # df_unit["r2_rrr_train"] = np.hstack(r2_rrr_train_units)
    # df_unit["r2_rrr_test"] = np.hstack(r2_rrr_test_units)
    df_unit["rrr_rank"] = np.hstack(rrr_rank_units)
    df_unit["rrr_reg"] = np.hstack(rrr_reg_units)
    # df_unit["r2_cca"] = np.hstack(r2_cca_units)
    # df_unit["r2_plsc"] = np.hstack(r2_plsc_units)
    # df_unit["r2_plsr"] = np.hstack(r2_plsr_units)
    df_unit["n_ca1"] = np.hstack(n_ca1_units)
    df_unit["n_target_cells"] = np.hstack(n_target_cells_units)
    df_unit["basepath"] = basepath
    df_unit["use_entire_session"] = use_entire_session
    df_unit["target_uid"] = np.hstack(target_uid)

    results = {"df": df, "df_unit": df_unit}
    return results


def load_results(save_path, verbose=False):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    df = pd.DataFrame()
    df_unit = pd.DataFrame()
    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if (results is None):
            continue
        if isinstance(results, dict):
            df = pd.concat([df, results['df']], ignore_index=True)
            df_unit = pd.concat([df, results['df_unit']], ignore_index=True)

    return df,df_unit
