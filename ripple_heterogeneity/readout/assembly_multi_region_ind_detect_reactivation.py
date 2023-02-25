from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import numpy as np
import nelpy as nel
import pandas as pd
import os
import pickle
import glob


def run(
    basepath: str,
    regions:str="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # brain regions to load
    cross_regions:tuple=(("CA1", "PFC"), ("CA1", "EC1|EC2|EC3|EC4|EC5|MEC")),
    ca1_layers:list=["Deep", "Superficial"],
    putativeCellType:str="Pyr",  # type of cells to load (can be multi ex. Pyr|Int)
    weight_dt:float=0.05,  # dt in seconds for binning st to get weights for each assembly
    restrict_to_nrem:bool=True, # to restrict sleep to nrem
    z_mat_dt:float=0.005 # AssemblyReact default is 2ms, but 5ms is faster to run
):
    """
    Run assembly reactivation for a single session
    Created for multi region analysis
    Inputs:
        basepath: str, path to session
        regions: str, brain regions to load
        cross_regions: tuple of tuples, regions to compare
        ca1_layers: list of str, deep or superficial
        putativeCellType: str, type of cells to load (can be multi ex. Pyr|Int)
        weight_dt: float, dt in seconds for binning st to get weights for each assembly
        restrict_to_nrem: bool, to restrict sleep to nrem
        z_mat_dt: float, default is 5ms
    Returns:
        response_df: pd.DataFrame, response of each assembly to pre post ripples
        peth: pd.DataFrame, peth of each assembly to pre post ripples

    """

    # initiate model and load data 
    assembly_react = assembly_reactivation.AssemblyReact(
        basepath=basepath,
        brainRegion=regions,
        weight_dt=weight_dt,
        putativeCellType=putativeCellType,
        z_mat_dt=z_mat_dt
    )
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
    assembly_react.cell_metrics.deepSuperficial = (
        assembly_react.cell_metrics.deepSuperficial.replace(np.nan, "unknown")
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
    pre = assembly_react.epochs[0][nel.EpochArray([pre_start, pre_start + 3600])]

    task = assembly_react.epochs[1][theta_epochs]

    post_start = assembly_react.epochs[2].start
    post = assembly_react.epochs[2][nel.EpochArray([post_start, post_start + 3600])]

    if restrict_to_nrem:
        pre = pre[nrem_epochs]
        post = post[nrem_epochs]

    if pre.isempty | task.isempty | post.isempty:
        return None

    response_df = pd.DataFrame()
    peth = pd.DataFrame()
    # iter over cross region types
    for cross_region in cross_regions:
        # check for current cross regions
        if (
            not assembly_react.cell_metrics.brainRegion.str.contains(
                cross_region[0]
            ).any()
            & assembly_react.cell_metrics.brainRegion.str.contains(
                cross_region[1]
            ).any()
        ):
            continue

        # iter over deep and superficial cells
        for deepSuperficial in ca1_layers:
            # check for deep sup cells
            if not assembly_react.cell_metrics.deepSuperficial.str.contains(
                deepSuperficial
            ).any():
                continue

            # subset cell label (deep or superficial and brain regions)
            query_region = cross_region[0] + "|" + cross_region[1]
            query_deep_sup = deepSuperficial + "|" + "unknown"
            idx = (
                assembly_react.cell_metrics.deepSuperficial.str.contains(query_deep_sup)
                & assembly_react.cell_metrics.brainRegion.str.contains(query_region)
            ).values

            if sum(idx) == 0:
                continue
            cell_metrics = assembly_react.cell_metrics[idx]

            # add only single cell class
            assembly_react.add_st(st.iloc[:, idx])

            # detect assembies in task theta epoch
            assembly_react.get_weights(epoch=task)

            if assembly_react.n_assemblies() == 0:
                continue

            # make sure assembly members represent cross region label
            keep_assembly = []
            n_ca1 = []
            n_cortex = []
            _, _, keep_assembly_, is_member = find_sig_assembly.main(
                assembly_react.patterns
            )
            for assembly_i in range(assembly_react.n_assemblies()):

                member_idx = is_member[assembly_i, :]

                cortex_check = (
                    cell_metrics[member_idx]
                    .brainRegion.str.contains(cross_region[1])
                    .sum()
                )
                ca1_check = (
                    cell_metrics[member_idx]
                    .brainRegion.str.contains(cross_region[0])
                    .sum()
                )
                ca1_layer_check = (
                    cell_metrics[member_idx].deepSuperficial == deepSuperficial
                ).sum()
                n_ca1.append(ca1_check)
                n_cortex.append(cortex_check)
                keep_assembly.append(
                    (cortex_check > 0) & (ca1_check > 0) & (ca1_layer_check > 0)
                )

            # only keep assemblies with cross-region members
            assembly_react.patterns = assembly_react.patterns[keep_assembly, :]
            n_ca1 = np.array(n_ca1)[keep_assembly]
            n_cortex = np.array(n_cortex)[keep_assembly]

            if assembly_react.n_assemblies() == 0:
                continue

            # to save time, calculate only during pre and post epochs
            assembly_act = assembly_react.get_assembly_act(epoch=pre + post)

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
            response_df_temp["cross_region_label"] = (
                deepSuperficial + "_" + cross_region[1]
            )
            response_df_temp["n_ca1"] = np.tile(n_ca1, 2)
            response_df_temp["n_cortex"] = np.tile(n_cortex, 2)
            response_df_temp["n_ca1_total"] = cell_metrics.brainRegion.str.contains(
                cross_region[0]
            ).sum()
            response_df_temp["n_cortex_total"] = cell_metrics.brainRegion.str.contains(
                cross_region[1]
            ).sum()
            response_df_temp["n_cells_total"] = assembly_react.st.n_active

            response_df = pd.concat([response_df, response_df_temp], ignore_index=True)

            peth_avg_pre = pd.DataFrame(
                index=time_lags,
                columns=np.arange(peth_avg_pre.shape[1]),
                data=peth_avg_pre,
            )
            peth_avg_post = pd.DataFrame(
                index=time_lags,
                columns=np.arange(peth_avg_post.shape[1]),
                data=peth_avg_post,
            )

            peth = pd.concat(
                [
                    peth,
                    pd.concat([peth_avg_pre, peth_avg_post], axis=1, ignore_index=True),
                ],
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

        results = pd.concat([results, results_["results_df"]], ignore_index=True)
        peth = pd.concat([peth, results_["peth"]], axis=1, ignore_index=True)

    return results, peth
