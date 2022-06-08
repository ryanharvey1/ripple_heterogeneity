import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import nelpy as nel
from quantities import s
import quantities as pq
from neo.core import SpikeTrain
from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import correlation_coefficient
from ripple_heterogeneity.utils import compress_repeated_epochs
import itertools
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

def get_downstream_data(st,cell_metrics,target_regions,ripple_epochs):
    target_idx = cell_metrics.brainRegion.str.contains(target_regions).values
    target_par = functions.get_participation(
        st.iloc[:, target_idx].data,
        ripple_epochs.starts,
        ripple_epochs.stops,
        par_type="counts",
    )
    return target_par


def run(
    basepath,  # path to data folder
    reference_region="CA1",  # reference region
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],  # regions to compare ref to
    min_cells=5,  # minimum number of cells per region
    restrict_task=False,  # restrict restriction_type to task epochs
    restriction_type="ripples",  # "ripples" or "NREMstate"
    ripple_expand=0.05,  # in seconds, how much to expand ripples
):

    st, cell_metrics = loading.load_spikes(basepath, brainRegion=["CA1", "PFC", "MEC"])
    cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cell_metrics)

    ripples = loading.load_ripples_events(basepath)
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T]).expand(
        0.1
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

    ca1_deep_par = functions.get_participation(
        st.iloc[:, ca1_deep_idx].data,
        ripple_epochs.starts,
        ripple_epochs.stops,
        par_type="binary",
    )
    ca1_sup_par = functions.get_participation(
        st.iloc[:, ca1_sup_idx].data,
        ripple_epochs.starts,
        ripple_epochs.stops,
        par_type="binary",
    )
    deep_sup_ratio = ca1_deep_par.sum(axis=0) / ca1_sup_par.sum(axis=0)

    y = np.log(deep_sup_ratio+1).reshape(-1,1)

    bad_idx = np.isinf(deep_sup_ratio) | np.isnan(deep_sup_ratio)
    y = y[~bad_idx]
    X = X[~bad_idx]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)

    reg = LinearRegression().fit(X_train, y_train)
    pred = reg.predict(X_test)

    mse = mean_squared_error(y_test, pred)

    # for region in target_regions:
