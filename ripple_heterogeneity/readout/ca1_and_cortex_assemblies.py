from ripple_heterogeneity.utils import functions, loading, compress_repeated_epochs
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.readout import ca1_assembly_downstream_psth
import pandas as pd
import numpy as np
from itertools import combinations
from scipy import signal
import glob
import os
import pickle


def compute_cross_correlogram(X, dt=1, window=0.5):
    """
    Cross-correlate two N-dimensional arrays.
    Input:
        X: N-dimensional array of shape  (n_signals, n_timepoints)
        dt: time step in seconds
        window: window size in seconds, output will be +- window
    Output:
        cross_correlogram: pandas dataframe with pairwise cross-correlogram
    """

    crosscorrs = {}
    pairs = list(combinations(np.arange(X.shape[0]), 2))
    for i, j in pairs:
        auc = signal.correlate(X[i], X[j])
        times = signal.correlation_lags(len(X[i]), len(X[j])) * dt
        # normalize by coeff
        normalizer = np.sqrt((X[i] ** 2).sum(axis=0) * (X[j] ** 2).sum(axis=0))
        auc /= normalizer

        crosscorrs[(i, j)] = pd.Series(index=times, data=auc, dtype="float32")
    crosscorrs = pd.DataFrame.from_dict(crosscorrs)

    if window is None:
        return crosscorrs
    else:
        return crosscorrs[(crosscorrs.index >= -window) & (crosscorrs.index <= window)]


def construct_assembly_df(react):
    # construct dataframe with assembly data
    assembly_df = pd.DataFrame()
    for r in react:
        # find significant assemblies and significant members
        _, _, keep_assembly, is_member = find_sig_assembly.main(r.patterns)
        # get assembly dataframe using code from previous analysis
        results = {"react": r}
        assembly_df_temp = ca1_assembly_downstream_psth.get_assembly_df(
            results, r.patterns, is_member
        )

        # add if the assembly is significant
        keep_assembly_df = pd.DataFrame(
            index=np.arange(r.patterns.shape[0]), data=keep_assembly
        )
        keep_assembly_dict = keep_assembly_df.to_dict()
        assembly_df_temp["assembly_keep"] = assembly_df_temp["assembly_n"]
        assembly_df_temp["assembly_keep"] = assembly_df_temp["assembly_keep"].map(
            keep_assembly_dict[0]
        )

        # add brain region
        assembly_df_temp["region"] = r.brainRegion
        assembly_df = pd.concat([assembly_df, assembly_df_temp], ignore_index=True)

    # add offset to assembly_n
    assem_offset = 0
    for region in assembly_df.region.unique():
        assembly_df.loc[assembly_df.region == region, "assembly_n"] = (
            assembly_df[assembly_df.region == region].assembly_n.values + assem_offset
        )
        assem_offset += assembly_df[assembly_df.region == region].assembly_n.max() + 1
    return assembly_df


def add_affiliation(assembly_df):

    labels = ["Superficial", "Deep"]
    assembly_df["affiliation"] = "unknown"

    # iterate over potential assembly labels
    for n in assembly_df.assembly_n.unique():
        # pull out current assembly and units that assembly members
        cur_assembly = assembly_df[
            (assembly_df.assembly_n == n)
            & (assembly_df.assembly_keep == True)
            & (assembly_df.is_member == True)
        ]
        if cur_assembly.shape[0] == 0:
            continue
        # if there no ca1 in the assembly, assign it to the cortex
        if not (cur_assembly.brainRegion.str.contains("CA1")).all():
            assembly_df.loc[
                (assembly_df.assembly_n == n),
                "affiliation",
            ] = "Cortex"
            continue
        # if there is a ca1 in the assembly, assign it to the deep or superficial
        for label in labels:
            if (cur_assembly.deepSuperficial == label).all():
                assembly_df.loc[
                    (assembly_df.assembly_n == n),
                    "affiliation",
                ] = label
    return assembly_df


def run(
    basepath,
    regions="CA1",  # not used
    target_regions=["CA1", "PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],
    putativeCellType="Pyr",
    weight_dt=0.02,  # dt in seconds for binning st to get weights for each assembly
    verbose=False,  # print out progress
    rip_exp_start=0.05,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.2,  # ripple expansion stop, in seconds, how much to expand ripples
    min_cells=5,  # minimum number of cells in analysis (n deep, n superficial, n cortex)
):

    epoch_df = loading.load_epoch(basepath)
    epoch_df = compress_repeated_epochs.main(epoch_df, epoch_name="sleep")
    idx, _ = functions.find_pre_task_post(epoch_df.environment)
    if idx is None:
        return None

    cm = pd.DataFrame()
    assembly_act = []
    assembly_act_all = []
    react = []
    for regions in target_regions:

        # initialize session
        m1 = assembly_reactivation.AssemblyReact(
            basepath,
            brainRegion=regions,
            putativeCellType=putativeCellType,
            weight_dt=weight_dt,
        )

        # load data
        m1.load_data()

        # restrict to pre/task/post epochs
        m1.restrict_epochs_to_pre_task_post()

        if m1.st.n_active < min_cells:
            continue

        # detect assemblies during the task
        m1.get_weights(m1.epochs[1])

        if len(m1.patterns) == 0:
            print("No patterns found")
            continue

        # extend ripples to include some extra time
        m1.ripples = m1.ripples.expand(rip_exp_start, direction="start")
        m1.ripples = m1.ripples.expand(rip_exp_stop, direction="stop")
        ac = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[2]])

        assembly_act.append(ac)
        assembly_act_all.append(ac.data)
        cm = pd.concat([cm, m1.cell_metrics], ignore_index=True)
        react.append(m1)

    if len(assembly_act_all) == 0:
        return None

    assembly_act_all = np.vstack(assembly_act_all)

    crosscorrs = compute_cross_correlogram(assembly_act_all, dt=m1.z_mat_dt)

    assembly_df = construct_assembly_df(react)

    assembly_df = add_affiliation(assembly_df)

    # deep_idx = assembly_df[assembly_df.affiliation == "Deep"].assembly_n.unique()
    # cortex_idx = assembly_df[assembly_df.affiliation == "Cortex"].assembly_n.unique()
    # # temp_crosscorrs = crosscorrs.loc[deep_idx,cortex_idx]
    # temp_crosscorrs = pd.DataFrame()
    # for val in deep_idx:
    #     temp_crosscorrs = pd.concat([temp_crosscorrs,crosscorrs[val][cortex_idx]],axis=1)
    results = {
        "assembly_df": assembly_df,
        "crosscorrs": crosscorrs,
        "assembly_act": assembly_act,
        "assembly_act_all": assembly_act_all,
        "cell_metrics": cm,
        "react": react,
    }
    return results


def load_results(save_path):
    sessions = glob.glob(os.path.join(save_path, "*.pkl"))

    deep_mec_cc = pd.DataFrame()
    deep_pfc_cc = pd.DataFrame()
    sup_mec_cc = pd.DataFrame()
    sup_pfc_cc = pd.DataFrame()

    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)

        if results is None:
            continue

        assembly_df = results.get("assembly_df")
        crosscorrs = results.get("crosscorrs")

        deep_idx = assembly_df[assembly_df.affiliation == "Deep"].assembly_n.unique()
        sup_idx = assembly_df[
            assembly_df.affiliation == "Superficial"
        ].assembly_n.unique()

        mec_idx = assembly_df[
            (assembly_df.affiliation == "Cortex")
            & assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC")
        ].assembly_n.unique()

        pfc_idx = assembly_df[
            (assembly_df.affiliation == "Cortex")
            & assembly_df.brainRegion.str.contains("PFC")
        ].assembly_n.unique()

        for val in deep_idx:
            try:
                deep_mec_cc = pd.concat([deep_mec_cc, crosscorrs[val][mec_idx]], axis=1,ignore_index=True)
            except:
                pass

        for val in deep_idx:
            try:
                deep_pfc_cc = pd.concat([deep_pfc_cc, crosscorrs[val][pfc_idx]], axis=1,ignore_index=True)
            except:
                pass
        for val in sup_idx:
            try:
                sup_mec_cc = pd.concat([sup_mec_cc, crosscorrs[val][mec_idx]], axis=1,ignore_index=True)
            except:
                pass
        for val in sup_idx:
            try:
                sup_pfc_cc = pd.concat([sup_pfc_cc, crosscorrs[val][pfc_idx]], axis=1,ignore_index=True)
            except:
                pass
    results = {
        "deep_mec": deep_mec_cc,
        "deep_pfc": deep_pfc_cc,
        "sup_mec": sup_mec_cc,
        "sup_pfc": sup_pfc_cc,
    }
    return results