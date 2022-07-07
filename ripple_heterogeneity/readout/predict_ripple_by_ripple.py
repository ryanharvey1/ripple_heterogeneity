import glob
import multiprocessing
import os
import pickle
from joblib import Parallel, delayed
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
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import ExtraTreesRegressor


def run_grid_search(X_train, y_train, n_grid=10, cv=5, max_rank=30):
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
        n_jobs=-1,
    )
    return grid_search.fit(X_train, y_train)


def get_data(basepath, target_regions, reference_region, rip_exp=0.5):
    """
    get_data: get the data for the analysis

    """
    st, cm = loading.load_spikes(
        basepath, brainRegion=[*target_regions, *reference_region]
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.peaks, ripples.peaks]).T]).expand(
        rip_exp
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
    rem_epochs = nel.EpochArray(state_dict["REMstate"])
    theta_epochs = nel.EpochArray(state_dict["THETA"])
    nontheta_epochs = nel.EpochArray(state_dict["nonTHETA"])

    return (
        st,
        cm,
        ripple_epochs,
        ep_epochs,
        ep_df,
        session_epoch,
        nrem_epochs,
        wake_epochs,
        rem_epochs,
        theta_epochs,
        nontheta_epochs,
        ripples,
    )


def mse_axis(y_test, y_pred, axis=1):
    """
    mse_axis: calculate the mean squared error along an axis
    """
    return ((y_test - y_pred) ** 2).mean(axis=axis)


def main_analysis(
    rip,
    ca1_sub,
    region,
    st,
    cm,
    min_cells,
    source_cell_type,
    target_cell_type,
    nrem_epochs,
    wake_epochs,
    rem_epochs,
    theta_epochs,
    nontheta_epochs,
    ep_epochs,
    ep_df,
    ripples,
    rip_i,
):

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
        return

    bst = st[rip].bin(ds=0.001).smooth(sigma=0.015)
    scaler = preprocessing.StandardScaler()
    X = scaler.fit_transform(bst.data)

    x = X[ca1_idx, :].T
    y = X[target_idx, :].T
    X_train, X_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=0.4,
        random_state=42,
        shuffle=True,
    )

    # simple linear regression
    reg = LinearRegression().fit(X_train, y_train)
    r2_train_lr = reg.score(X_train, y_train)
    r2_test_lr = reg.score(X_test, y_test)
    mse_train_lr = mean_squared_error(y_train, reg.predict(X_train))
    mse_test_lr = mean_squared_error(y_test, reg.predict(X_test))
    # predict across whole ripple and estimate error
    mse_time_lr = mse_axis(y, reg.predict(x), axis=1)

    reg = ExtraTreesRegressor().fit(X_train, y_train)
    r2_train_et = reg.score(X_train, y_train)
    r2_test_et = reg.score(X_test, y_test)
    mse_train_et = mean_squared_error(y_train, reg.predict(X_train))
    mse_test_et = mean_squared_error(y_test, reg.predict(X_test))
    # predict across whole ripple and estimate error
    mse_time_et = mse_axis(y, reg.predict(x), axis=1)

    # reduced rank regression
    # reg = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()
    # if grid_search:
    #     grid_search_result = run_grid_search(
    #         X_train, y_train, n_grid=n_grid, cv=cv, max_rank=max_rank
    #     )
    # if grid_search:
    #     reg.rank = int(grid_search_result.best_params_["rank"])
    # else:
    #     reg.rank = rank
    # reg.reg = rrrr_reg
    # reg.fit(X_train, y_train)
    # r2_train_rrrr.append(reg.score(X_train, y_train))
    # r2_test_rrrr.append(reg.score(X_test, y_test))
    # mse_train_rrrr.append(mean_squared_error(y_train, reg.predict(X_train)))
    # mse_test_rrrr.append(mean_squared_error(y_test, reg.predict(X_test)))
    # # predict across whole ripple and estimate error
    # mse_time_rrrr.append(mse_axis(y, reg.predict(x), axis=1))

    # get metadata
    if not (rip & nrem_epochs).isempty:
        state = "NREM"
    elif not (rip & wake_epochs).isempty:
        state = "Wake"
    elif not (rip & rem_epochs).isempty:
        state = "REM"
    elif not (rip & theta_epochs).isempty:
        state = "Theta"
    elif not (rip & nontheta_epochs).isempty:
        state = "NonTheta"

    for ep_i, ep in enumerate(ep_epochs):
        if not (rip & ep).isempty:
            environment = ep_df.loc[ep_i].environment
            name = ep_df.loc[ep_i].name
            epoch_i = ep_i
            break

    # rrr_rank.append(reg.rank)
    n_target = sum(target_idx)
    n_ca1 = sum(ca1_idx)
    targ_reg = region
    ca1_sub_layer = ca1_sub

    rip_n = rip_i
    ripple_duation = ripples.duration.loc[rip_i]
    ripple_start = ripples.start.loc[rip_i]
    ripple_stop = ripples.stop.loc[rip_i]
    ripple_peak = ripples.peaks.loc[rip_i]
    ripple_amp = ripples.amplitude.loc[rip_i]
    ripple_freq = ripples.frequency.loc[rip_i]


    df = {
        "state": state,
        "epoch_i": epoch_i,
        "environment": environment,
        "name": name,
        "targ_reg": targ_reg,
        "ca1_sub_layer": ca1_sub_layer,
        "rip_n": rip_n,
        "ripple_duation": ripple_duation,
        "ripple_start": ripple_start,
        "ripple_stop": ripple_stop,
        "ripple_peak": ripple_peak,
        "ripple_amp": ripple_amp,
        "ripple_freq": ripple_freq,
        "r2_train_lr": r2_train_lr,
        "r2_test_lr": r2_test_lr,
        "mse_train_lr": mse_train_lr,
        "mse_test_lr": mse_test_lr,
        # "mse_time_lr": mse_time_lr,
        "r2_train_et": r2_train_et,
        "r2_test_et": r2_test_et,
        "mse_train_et": mse_train_et,
        "mse_test_et": mse_test_et,
        # "mse_time_et": mse_time_et,
        "n_target": n_target,
        "n_ca1": n_ca1,
    }

    return df


def run(
    basepath,  # path to data folder
    reference_region=["CA1"],  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    rank=10,  # rank of the reduced rank regressor (not used)
    rrrr_reg=1,  # regularization parameter
    source_cell_type="Pyr",  # source cell type
    target_cell_type=None,  # cell type to use for target cells
    n_grid=20,  # number of grid search parameters to use
    cv=5,  # number of cross validation folds
    max_rank=30,  # maximum rank to use in the reduced rank regressor
    use_entire_session=False,  # use entire session or just pre task post
    grid_search=True,  # use grid search to find the best rank
    rip_exp=0.5,  # expansion for ripples (center plus/minus this amount)
    parallel=False,  # use parallel processing
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
        rem_epochs,
        theta_epochs,
        nontheta_epochs,
        ripples,
    ) = get_data(basepath, target_regions, reference_region, rip_exp=rip_exp)
    if st is None:
        return None


    # if parallel:
    # get number of cores
    num_cores = multiprocessing.cpu_count()

    # iterate over ca1 sublayers regions
    for ca1_sub in ["Deep", "Superficial"]:
        # iterate over target regions
        for region in target_regions:    
            # run in parallel
            processed_list = Parallel(n_jobs=num_cores)(
                delayed(main_analysis)(
                        rip,
                        ca1_sub,
                        region,
                        st,
                        cm,
                        min_cells,
                        source_cell_type,
                        target_cell_type,
                        nrem_epochs,
                        wake_epochs,
                        rem_epochs,
                        theta_epochs,
                        nontheta_epochs,
                        ep_epochs,
                        ep_df,
                        ripples,
                        rip_i,
                )
                for rip_i, rip in enumerate(ripple_epochs)
            )

    return processed_list
    # df = []
    # # iterate over ca1 sublayers regions
    # for ca1_sub in ["Deep", "Superficial"]:
    #     # iterate over target regions
    #     for region in target_regions:
            # for rip_i, rip in enumerate(ripple_epochs):
            #     df.append(main_analysis(
            #         rip,
            #         ca1_sub,
            #         region,
            #         st,
            #         cm,
            #         min_cells,
            #         source_cell_type,
            #         target_cell_type,
            #         nrem_epochs,
            #         wake_epochs,
            #         rem_epochs,
            #         theta_epochs,
            #         nontheta_epochs,
            #         ep_epochs,
            #         ep_df,
            #         ripples,
            #         rip_i,
            #     ))

    # pd.DataFrame.from_dict(df[0], orient="index")


    # initialize output vars (not the best way to do this, long format is better)
    # r2_train_lr = []
    # r2_test_lr = []
    # mse_train_lr = []
    # mse_test_lr = []
    # mse_time_lr = []

    # r2_train_et = []
    # r2_test_et = []
    # mse_train_et = []
    # mse_test_et = []
    # mse_time_et = []

    # r2_train_rrrr = []
    # r2_test_rrrr = []
    # mse_train_rrrr = []
    # mse_test_rrrr = []
    # mse_time_rrrr = []

    # rrr_rank = []
    # state = []
    # n_target = []
    # n_ca1 = []
    # targ_reg = []
    # ca1_sub_layer = []
    # environment = []
    # name = []
    # epoch_i = []

    # rip_n = []
    # ripple_duation = []

    # iterate over ripples
    # for rip_i, rip in enumerate(ripple_epochs):
    #     # iterate over ca1 sublayers regions
    #     for ca1_sub in ["Deep", "Superficial"]:
    #         # iterate over target regions
    #         for region in target_regions:
    #             if sum(cm.brainRegion.str.contains(region).values) < min_cells:
    #                 continue

    #             # get index of ca1 cells
    #             ca1_idx = (
    #                 cm.brainRegion.str.contains("CA1").values
    #                 & (cm.deepSuperficial == ca1_sub)
    #                 & (cm.putativeCellType.str.contains(source_cell_type))
    #             )

    # # get index of target cells
    # if target_cell_type is not None:
    #     target_idx = (
    #         cm.brainRegion.str.contains(region).values
    #         & cm.putativeCellType.str.contains(target_cell_type).values
    #     )
    # else:
    #     target_idx = cm.brainRegion.str.contains(region).values

    # if (sum(ca1_idx) < min_cells) | (sum(target_idx) < min_cells):
    #     continue

    # bst = st[rip].bin(ds=0.001).smooth(sigma=0.015)
    # scaler = preprocessing.StandardScaler()
    # X = scaler.fit_transform(bst.data)

    # x = X[ca1_idx, :].T
    # y = X[target_idx, :].T
    # X_train, X_test, y_train, y_test = train_test_split(
    #     x,
    #     y,
    #     test_size=0.4,
    #     random_state=42,
    #     shuffle=True,
    # )

    # # simple linear regression
    # reg = LinearRegression().fit(X_train, y_train)
    # r2_train_lr.append(reg.score(X_train, y_train))
    # r2_test_lr.append(reg.score(X_test, y_test))
    # mse_train_lr.append(mean_squared_error(y_train, reg.predict(X_train)))
    # mse_test_lr.append(mean_squared_error(y_test, reg.predict(X_test)))
    # # predict across whole ripple and estimate error
    # mse_time_lr.append(mse_axis(y, reg.predict(x), axis=1))

    # reg = ExtraTreesRegressor().fit(X_train, y_train)
    # r2_train_et.append(reg.score(X_train, y_train))
    # r2_test_et.append(reg.score(X_test, y_test))
    # mse_train_et.append(mean_squared_error(y_train, reg.predict(X_train)))
    # mse_test_et.append(mean_squared_error(y_test, reg.predict(X_test)))
    # # predict across whole ripple and estimate error
    # mse_time_et.append(mse_axis(y, reg.predict(x), axis=1))

    # reduced rank regression
    # reg = kernel_reduced_rank_ridge_regression.ReducedRankRegressor()
    # if grid_search:
    #     grid_search_result = run_grid_search(
    #         X_train, y_train, n_grid=n_grid, cv=cv, max_rank=max_rank
    #     )
    # if grid_search:
    #     reg.rank = int(grid_search_result.best_params_["rank"])
    # else:
    #     reg.rank = rank
    # reg.reg = rrrr_reg
    # reg.fit(X_train, y_train)
    # r2_train_rrrr.append(reg.score(X_train, y_train))
    # r2_test_rrrr.append(reg.score(X_test, y_test))
    # mse_train_rrrr.append(mean_squared_error(y_train, reg.predict(X_train)))
    # mse_test_rrrr.append(mean_squared_error(y_test, reg.predict(X_test)))
    # # predict across whole ripple and estimate error
    # mse_time_rrrr.append(mse_axis(y, reg.predict(x), axis=1))

    # # get metadata
    # if not (rip & nrem_epochs).isempty:
    #     state.append("NREM")
    # elif not (rip & wake_epochs).isempty:
    #     state.append("Wake")
    # elif not (rip & rem_epochs).isempty:
    #     state.append("REM")
    # elif not (rip & theta_epochs).isempty:
    #     state.append("Theta")
    # elif not (rip & nontheta_epochs).isempty:
    #     state.append("NonTheta")

    # for ep_i, ep in enumerate(ep_epochs):
    #     if not (rip & ep).isempty:
    #         environment.append(ep_df.loc[ep_i].environment)
    #         name.append(ep_df.loc[ep_i].name)
    #         epoch_i.append(ep_i)
    #         break

    # # rrr_rank.append(reg.rank)
    # n_target.append(sum(target_idx))
    # n_ca1.append(sum(ca1_idx))
    # targ_reg.append(region)
    # ca1_sub_layer.append(ca1_sub)

    # rip_n.append(rip_i)
    # ripple_duation.append(ripples.duration.loc[rip_i])

    # create a dataframe
    # df = pd.DataFrame()
    # df["state"] = np.hstack(state)
    # df["environment"] = np.hstack(environment)
    # df["name"] = np.hstack(name)
    # df["epoch_i"] = np.hstack(epoch_i)
    # df["n_target"] = np.hstack(n_target)
    # df["n_ca1"] = np.hstack(n_ca1)
    # df["targ_reg"] = np.hstack(targ_reg)
    # df["ca1_sub_layer"] = np.hstack(ca1_sub_layer)
    # df["rrr_rank"] = np.hstack(rrr_rank)
    # df["r2_train_lr"] = np.hstack(r2_train_lr)
    # df["r2_test_lr"] = np.hstack(r2_test_lr)
    # df["mse_train_lr"] = np.hstack(mse_train_lr)
    # df["mse_test_lr"] = np.hstack(mse_test_lr)
    # df["r2_train_et"] = np.hstack(r2_train_et)
    # df["r2_test_et"] = np.hstack(r2_test_et)
    # df["mse_train_et"] = np.hstack(mse_train_et)
    # df["mse_test_et"] = np.hstack(mse_test_et)
    # df["r2_train_rrrr"] = np.hstack(r2_train_rrrr)
    # df["r2_test_rrrr"] = np.hstack(r2_test_rrrr)
    # df["mse_train_rrrr"] = np.hstack(mse_train_rrrr)
    # df["mse_test_rrrr"] = np.hstack(mse_test_rrrr)
    # df["n_epochs"] = np.hstack(ep_epochs.n_intervals)
    # df["n_ripples"] = np.hstack(ripple_epochs.n_intervals)
    # df["basepath"] = basepath

    # df["n_samples"] = np.hstack(n_samples)
    # df["n_features"] = np.hstack(n_features)
    # df["n_grid"] = np.hstack(n_grid)
    # df["max_rank"] = np.hstack(max_rank)
    # df["cv"] = np.hstack(cv)

    # results = {
    #     "df": df,
    #     "mse_time_lr": mse_time_lr,
    #     "mse_time_et": mse_time_et,
    #     # "mse_time_rrrr": mse_time_rrrr,
    # }
    # return results
