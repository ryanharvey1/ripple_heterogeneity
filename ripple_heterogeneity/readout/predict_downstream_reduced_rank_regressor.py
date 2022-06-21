import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel
from ripple_heterogeneity.utils import compress_repeated_epochs
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from ripple_heterogeneity.utils import reduced_rank_regressor
from scipy import around
from scipy import size
from scipy.linalg import norm


def sqerr(matrix1, matrix2):
    """Squared error (frobenius norm of diff) between two matrices."""
    return around(pow(norm(matrix1 - matrix2, "fro"), 2) / size(matrix2, 0), 5)


def shuffle_data(X, y, rank, reg, n_shuff=1000):
    testing_error = []
    for i in range(n_shuff):
        idx = np.random.permutation(X.shape[0])
        X_train, X_test, y_train, y_test = train_test_split(
            X[idx, :], y, test_size=0.4, random_state=42
        )
        regressor = reduced_rank_regressor.ReducedRankRegressor(
            X_train, y_train, rank, reg
        )
        # get model performance
        testing_error.append(sqerr(regressor.predict(X_test), y_test))

    return testing_error


def run(
    basepath,  # path to data folder
    reference_region=["CA1"],  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    ripple_expand=0.1,  # in seconds, how much to expand ripples
    min_ripples=10,  # minimum number of ripples per epoch
    n_shuff=1000,  # number of shuffles to do
    rank=10,  # rank of the reduced rank regressor
    reg=1e-6,  # regularization parameter
    target_cell_type=None,  # cell type to use for target cells
):

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
    # locate pre task post structure
    idx, _ = functions.find_pre_task_post(ep_df.environment)
    if idx is None:
        return None

    ep_df = ep_df[idx]
    ep_epochs = nel.EpochArray([np.array([ep_df.startTime, ep_df.stopTime]).T])

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

    scaler = preprocessing.StandardScaler()

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
                    & (cm.putativeCellType.str.contains("Pyr"))
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
                regressor = reduced_rank_regressor.ReducedRankRegressor(
                    X_train, y_train, rank, reg
                )
                # get model performance
                training_error.append(sqerr(regressor.predict(X_train), y_train))
                testing_error.append(sqerr(regressor.predict(X_test), y_test))

                # get metadata
                n_x_components.append(X.shape[1])
                epoch.append(ep_df.environment.iloc[ep_i])
                epoch_i.append(ep_i)
                targ_reg.append(region)
                ca1_sub_layer.append(ca1_sub)
                n_ca1.append(sum(ca1_idx))
                n_target_cells.append(sum(cm.brainRegion.str.contains(region).values))
                # get vars for prediction gain
                error_shuff = shuffle_data(
                    X[ca1_idx, :].T, X[target_idx, :].T, rank, reg, n_shuff=n_shuff
                )
                median_error_shuff.append(np.median(error_shuff))
                mean_error_shuff.append(np.mean(error_shuff))

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
    df["mean_error_shuff"] = np.hstack(mean_error_shuff)
    df["median_error_shuff"] = np.hstack(median_error_shuff)
    df["n_ca1"] = np.hstack(n_ca1)
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
