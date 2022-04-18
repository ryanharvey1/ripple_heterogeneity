import numpy as np
import pandas as pd
import itertools
from scipy import stats
import nelpy as nel
import os
import multiprocessing
from joblib import Parallel, delayed
from ripple_heterogeneity.utils import functions,loading
from ripple_heterogeneity.utils import add_new_deep_sup
from ripple_heterogeneity.assembly import assembly_run

def pairwise_corr(unit_mat):
    x = np.arange(0, unit_mat.shape[0])
    c = np.array(list(itertools.combinations(x, 2)))
    rho = []
    pval = []
    for i, s in enumerate(c):
        rho_, pval_ = stats.spearmanr(unit_mat[s[0], :], unit_mat[s[1], :])
        rho.append(rho_)
        pval.append(pval_)
    return rho, pval, c


def pairwise_cross_corr(spks, binsize=0.001, nbins=100):
    # Get unique combo without repeats
    x = np.arange(0, spks.shape[0])
    c = np.array(list(itertools.combinations(x, 2)))
    # prepare a pandas dataframe to receive the data
    times = np.arange(0, binsize * (nbins + 1), binsize) - (nbins * binsize) / 2
    crosscorrs = pd.DataFrame(index=times, columns=np.arange(len(c)))

    # Now we can iterate over spikes
    for i, s in enumerate(c):
        # Calling the crossCorr function
        crosscorrs[i] = functions.crossCorr(spks[s[0]], spks[s[1]], binsize, nbins)
    return crosscorrs


def get_group_cor_vectors(df, temp_df, basepath):
    deep = []
    sup = []
    member_rho = []
    non_member_rho = []
    member_deep = []
    member_sup = []
    non_member_deep = []
    non_member_sup = []
    member_deep_sup = []
    non_member_deep_sup = []

    assembly_id = {}
    assembly_id["deep"] = []
    assembly_id["sup"] = []
    assembly_id["member"] = []
    assembly_id["non_member"] = []
    assembly_id["member_deep"] = []
    assembly_id["member_sup"] = []
    assembly_id["non_member_deep"] = []
    assembly_id["non_member_sup"] = []
    assembly_id["member_deep_sup"] = []
    assembly_id["non_member_deep_sup"] = []

    for assembly_n in df[df.basepath == basepath].assembly_n.unique():
        current_df = df[(df.basepath == basepath) & (df.assembly_n == assembly_n)]
        member_uid = current_df[current_df.membership == True].UID.values
        nonmember_uid = current_df[current_df.membership == False].UID.values

        deep.append(
            temp_df[
                (temp_df.layer_ref == "Deep") & (temp_df.layer_target == "Deep")
            ].rho.values
        )
        sup.append(
            temp_df[
                (temp_df.layer_ref == "Superficial")
                & (temp_df.layer_target == "Superficial")
            ].rho.values
        )
        member_rho.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, member_uid))
                & (np.in1d(temp_df.UID_target, member_uid))
            ].rho.values
        )
        non_member_rho.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, nonmember_uid))
                & (np.in1d(temp_df.UID_target, nonmember_uid))
            ].rho.values
        )
        member_deep.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, member_uid))
                & (np.in1d(temp_df.UID_target, member_uid))
                & (temp_df.layer_ref == "Deep")
                & (temp_df.layer_target == "Deep")
            ].rho.values
        )
        member_sup.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, member_uid))
                & (np.in1d(temp_df.UID_target, member_uid))
                & (temp_df.layer_ref == "Superficial")
                & (temp_df.layer_target == "Superficial")
            ].rho.values
        )
        non_member_deep.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, nonmember_uid))
                & (np.in1d(temp_df.UID_target, nonmember_uid))
                & (temp_df.layer_ref == "Deep")
                & (temp_df.layer_target == "Deep")
            ].rho.values
        )
        non_member_sup.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, nonmember_uid))
                & (np.in1d(temp_df.UID_target, nonmember_uid))
                & (temp_df.layer_ref == "Superficial")
                & (temp_df.layer_target == "Superficial")
            ].rho.values
        )
        member_deep_sup.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, member_uid))
                & (np.in1d(temp_df.UID_target, member_uid))
                & (
                    (temp_df.layer_ref == "Deep")
                    & (temp_df.layer_target == "Superficial")
                    | (temp_df.layer_ref == "Superficial")
                    & (temp_df.layer_target == "Deep")
                )
            ].rho.values
        )
        non_member_deep_sup.append(
            temp_df[
                (np.in1d(temp_df.UID_ref, nonmember_uid))
                & (np.in1d(temp_df.UID_target, nonmember_uid))
                & (
                    (temp_df.layer_ref == "Deep")
                    & (temp_df.layer_target == "Superficial")
                    | (temp_df.layer_ref == "Superficial")
                    & (temp_df.layer_target == "Deep")
                )
            ].rho.values
        )

        assembly_id["deep"].append([assembly_n] * len(deep[-1]))
        assembly_id["sup"].append([assembly_n] * len(sup[-1]))
        assembly_id["member"].append([assembly_n] * len(member_rho[-1]))
        assembly_id["non_member"].append([assembly_n] * len(non_member_rho[-1]))
        assembly_id["member_deep"].append([assembly_n] * len(member_deep[-1]))
        assembly_id["member_sup"].append([assembly_n] * len(member_sup[-1]))
        assembly_id["non_member_deep"].append([assembly_n] * len(non_member_deep[-1]))
        assembly_id["non_member_sup"].append([assembly_n] * len(non_member_sup[-1]))
        assembly_id["member_deep_sup"].append([assembly_n] * len(member_deep_sup[-1]))
        assembly_id["non_member_deep_sup"].append(
            [assembly_n] * len(non_member_deep_sup[-1])
        )

    deep = np.hstack(deep)
    sup = np.hstack(sup)
    member = np.hstack(member_rho)
    non_member = np.hstack(non_member_rho)
    member_deep = np.hstack(member_deep)
    member_sup = np.hstack(member_sup)
    non_member_deep = np.hstack(non_member_deep)
    non_member_sup = np.hstack(non_member_sup)
    member_deep_sup = np.hstack(member_deep_sup)
    non_member_deep_sup = np.hstack(non_member_deep_sup)

    for key_ in assembly_id.keys():
        assembly_id[key_] = np.hstack(assembly_id[key_])

    return (
        deep,
        sup,
        member,
        non_member,
        member_deep,
        member_sup,
        non_member_deep,
        non_member_sup,
        member_deep_sup,
        non_member_deep_sup,
        assembly_id,
    )


def get_pairwise_corrs(basepath, epoch):
    """
    Get pairwise correlations for all assemblies.
    Input:
        basepath: path to the directory containing the data
        epoch: epoch of interest
    Output:
        Data frame with pairwise correlations for all assemblies:
            deep: deep layer correlations
            sup: superficial layer correlations
            member: member-member correlations
            non_member: non-member-non-member correlations
            member_deep: member-member deep layer correlations
            member_sup: member-member superficial layer correlations
            non_member_deep: non-member-non-member deep layer correlations
            non_member_sup: non-member-non-member superficial layer correlations
            member_deep_sup: member-member deep-superficial layer correlations
            non_member_deep_sup: non-member-non-member deep-superficial layer correlations
            assembly_id: assembly id for each correlation
    """
    cell_metrics, data, ripples, fs_dat = assembly_run.load_basic_data(basepath)
    restrict_idx = (
        (cell_metrics.putativeCellType == "Pyramidal Cell")
        & (cell_metrics.brainRegion.str.contains("CA1"))
        & (cell_metrics.bad_unit == False)
    )

    # restrict cell metrics
    cell_metrics = cell_metrics[restrict_idx]

    # add deep sup classification
    cell_metrics = add_new_deep_sup.add_new_deep_sup_class(cell_metrics)

    st_unit = nel.SpikeTrainArray(
        timestamps=np.array(data["spikes"], dtype=object)[restrict_idx], fs=fs_dat
    )
    ripple_epochs = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])

    if epoch is not None:
        epochs = loading.load_epoch(basepath)
        idx = functions.find_pre_task_post(epochs.environment)
        epochs = epochs[idx[0]]
        behav_epochs = nel.EpochArray([np.array([epochs.startTime, epochs.stopTime]).T])
        if epoch == "pre":
            ripple_epochs = ripple_epochs[behav_epochs[0]]
            st_unit = st_unit[behav_epochs[0]]
        elif epoch == "task":
            ripple_epochs = ripple_epochs[behav_epochs[1]]
            st_unit = st_unit[behav_epochs[1]]
        elif epoch == "post":
            ripple_epochs = ripple_epochs[behav_epochs[2]]
            st_unit = st_unit[behav_epochs[2]]
        else:
            raise ValueError("Invalid epoch")

    spk_count_rip = functions.get_participation(
        st_unit.data, ripple_epochs.starts, ripple_epochs.stops
    )

    rho, pval, c = pairwise_corr(spk_count_rip)

    temp_df = pd.DataFrame()
    temp_df["ref"] = c[:, 0]
    temp_df["target"] = c[:, 1]
    temp_df["UID_ref"] = cell_metrics.UID.iloc[temp_df["ref"]].values
    temp_df["UID_target"] = cell_metrics.UID.iloc[temp_df["target"]].values
    temp_df["rho"] = rho
    temp_df["pval"] = pval
    temp_df["layer_ref"] = cell_metrics.deepSuperficial.iloc[temp_df["ref"]].values
    temp_df["layer_target"] = cell_metrics.deepSuperficial.iloc[
        temp_df["target"]
    ].values
    temp_df["basepath"] = basepath

    return temp_df


def get_and_organize_pairwise_corrs(basepath, df, epoch):

    temp_df = get_pairwise_corrs(basepath, epoch)
    (
        deep,
        sup,
        member,
        non_member,
        member_deep,
        member_sup,
        non_member_deep,
        non_member_sup,
        member_deep_sup,
        non_member_deep_sup,
        assembly_id,
    ) = get_group_cor_vectors(df, temp_df, basepath)

    df_save = pd.DataFrame()

    df_save["rho"] = np.hstack(
        [
            deep,
            sup,
            member,
            non_member,
            member_deep,
            member_sup,
            non_member_deep,
            non_member_sup,
            member_deep_sup,
            non_member_deep_sup,
        ]
    )

    df_save["label"] = np.hstack(
        [
            ["deep"] * len(deep),
            ["sup"] * len(sup),
            ["member"] * len(member),
            ["non_member"] * len(non_member),
            ["member_deep"] * len(member_deep),
            ["member_sup"] * len(member_sup),
            ["non_member_deep"] * len(non_member_deep),
            ["non_member_sup"] * len(non_member_sup),
            ["member_deep_sup"] * len(member_deep_sup),
            ["non_member_deep_sup"] * len(non_member_deep_sup),
        ]
    )

    df_save["assembly_id"] = np.hstack(
        [
            assembly_id["deep"],
            assembly_id["sup"],
            assembly_id["member"],
            assembly_id["non_member"],
            assembly_id["member_deep"],
            assembly_id["member_sup"],
            assembly_id["non_member_deep"],
            assembly_id["non_member_sup"],
            assembly_id["member_deep_sup"],
            assembly_id["non_member_deep_sup"],
        ]
    )
    df_save["basepath"] = basepath

    return df_save


def session_loop(basepath, df, save_path, epoch):

    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".csv"
    )
    if os.path.exists(save_file):
        return

    df_save = get_and_organize_pairwise_corrs(basepath, df, epoch)
    df_save.to_csv(save_file)


def assembly_corr_run(df, save_path, parallel=True, epoch=None):
    """
    Run the pairwise correlation analysis on all assemblies in the dataframe.
    Input:
        df: dataframe with assembly information
        save_path: path to save the results
        parallel: whether to run in parallel
        epoch: which epoch to use [pre, task, post] (assumes sleep,linear,sleep)
    """
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(session_loop)(basepath, df, save_path, epoch)
            for basepath in basepaths
        )
    else:
        for basepath in basepaths:
            print(basepath)
            session_loop(basepath, df, save_path, epoch)
