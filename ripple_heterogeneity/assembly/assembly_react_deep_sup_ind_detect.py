from ripple_heterogeneity.assembly import assembly_reactivation
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import numpy as np
import nelpy as nel
import pandas as pd
import os
import pickle
import glob


def run(basepath: str):

    # initiate model
    assembly_react = assembly_reactivation.AssemblyReact(basepath=basepath)
    assembly_react.load_data()

    if assembly_react.isempty:
        return None

    # if this fails, there is no pre task post structure in session
    try:
        assembly_react.restrict_epochs_to_pre_task_post()
    except:
        return None

    assembly_react.cell_metrics = (
        add_new_deep_sup.deep_sup_from_deepSuperficialDistance(
            assembly_react.cell_metrics
        )
    )

    # assign copy so I can subset and add to assembly_react
    st = assembly_react.st.copy()

    state_dict = loading.load_SleepState_states(basepath)
    theta_epochs = nel.EpochArray(state_dict["THETA"])
    nrem_epochs = nel.EpochArray(state_dict["NREMstate"])

    # define pre / task / post epochs
    #   restrict to first hour of sleep
    #   restrict task to theta epochs and sleep to nrem
    pre_start = assembly_react.epochs[0].start
    pre = assembly_react.epochs[0][nel.EpochArray([pre_start, pre_start + 3600])][
        nrem_epochs
    ]

    task = assembly_react.epochs[1][theta_epochs]

    post_start = assembly_react.epochs[2].start
    post = assembly_react.epochs[2][nel.EpochArray([post_start, post_start + 3600])][
        nrem_epochs
    ]

    if pre.isempty | task.isempty | post.isempty:
        return None

    response_df = pd.DataFrame()
    peth = pd.DataFrame()

    # iter over deep and superficial cells
    for deepSuperficial in ["Deep", "Superficial"]:

        # subset cell label (deep or superficial)
        idx = (assembly_react.cell_metrics.deepSuperficial == deepSuperficial).values

        if sum(idx) == 0:
            continue

        # add only single cell class
        assembly_react.add_st(st.iloc[:, idx])

        # detect assembies in task theta epoch
        assembly_react.get_weights(epoch=task)

        if assembly_react.n_assemblies() == 0:
            continue

        # to save time, calculate only during pre and post epochs
        assembly_act = assembly_react.get_assembly_act()

        peth_avg_pre, time_lags = functions.event_triggered_average(
            assembly_act.abscissa_vals,
            assembly_act.data.T,
            assembly_react.ripples[pre].starts,
            sampling_rate=assembly_act.fs,
            window=[-0.5, 0.5],
        )
        peth_avg_post, time_lags = functions.event_triggered_average(
            assembly_act.abscissa_vals,
            assembly_act.data.T,
            assembly_react.ripples[post].starts,
            sampling_rate=assembly_act.fs,
            window=[-0.5, 0.5],
        )

        response_df_temp = pd.DataFrame()
        # get mean of peth 0 to 100ms
        response_df_temp["response"] = np.hstack(
            [
                np.nanmean(
                    peth_avg_pre[(time_lags > 0) & (time_lags < 0.1), :], axis=0
                ),
                np.nanmean(
                    peth_avg_post[(time_lags > 0) & (time_lags < 0.1), :], axis=0
                ),
            ]
        )
        response_df_temp["assembly_n"] = np.hstack(
            [np.arange(peth_avg_pre.shape[1]), np.arange(peth_avg_pre.shape[1])]
        )

        # add labels for epochs
        response_df_temp["epoch"] = np.hstack(
            [["pre"] * peth_avg_pre.shape[1], ["post"] * peth_avg_pre.shape[1]]
        )

        response_df_temp["deepSuperficial"] = deepSuperficial
        response_df_temp["n_cells"] = assembly_react.st.n_active

        response_df = pd.concat([response_df, response_df_temp], ignore_index=True)

        peth_avg_pre = pd.DataFrame(
            index=time_lags, columns=np.arange(peth_avg_pre.shape[1]), data=peth_avg_pre
        )
        peth_avg_post = pd.DataFrame(
            index=time_lags,
            columns=np.arange(peth_avg_post.shape[1]),
            data=peth_avg_post,
        )

        peth = pd.concat(
            [peth, pd.concat([peth_avg_pre, peth_avg_post], axis=1, ignore_index=True)],
            axis=1,
            ignore_index=True,
        )

    response_df["basepath"] = basepath

    results = {"results_df": response_df, "peth": peth}

    return results


def load_results(save_path: str, verbose: bool = False):
    """
    load_results: load results (pandas dataframe) from a pickle file

    This is the most basic results loader and
        **will only work if your output was a pandas dataframe (long format)**

    This will have to be adapted if your output was more complicated, but you can
        use this function as an example.
    """

    if not os.path.exists(save_path):
        raise ValueError(f"folder {save_path} does not exist")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    results = pd.DataFrame()
    peth = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results_ = pickle.load(f)
        if results_ is None:
            continue
        if results_["results_df"].shape[0] == 0:
            continue
        if (results_["peth"].shape[0] == 501):
            test = 0
        results = pd.concat([results, results_["results_df"]], ignore_index=True)
        peth = pd.concat([peth, results_["peth"]], axis=1)

    return results, peth
