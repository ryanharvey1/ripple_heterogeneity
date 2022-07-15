import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.utils import add_new_deep_sup
import warnings


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",
    putativeCellType="Pyr", 
    weight_dt=0.1, # dt in seconds for binning st to get weights for each assembly
    verbose=False, # print out progress
    rip_expand=0.05, # expand the ripple region by this amount (not used)
    rip_exp_start=0.05,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.2,  # ripple expansion stop, in seconds, how much to expand ripples
):
    """
    Gets the pre and post assembly strengths
    """
    # initialize session
    m1 = assembly_reactivation.AssemblyReact(
        basepath,
        brainRegion=regions,
        putativeCellType=putativeCellType,
        weight_dt=weight_dt,
    )

    # load data
    m1.load_data()

    # extend ripples to include some extra time
    m1.ripples = m1.ripples.expand(rip_exp_start, direction="start")
    m1.ripples = m1.ripples.expand(rip_exp_stop, direction="stop")

    # check if no cells were found
    if m1.cell_metrics.shape[0] == 0:
        return None

    # restrict to pre/task/post epochs
    try:
        m1.restrict_epochs_to_pre_task_post()
    except:
        print("No pre/task/post epochs found")
        return None
    # get weights for task outside ripples
    # % (TODO: use more robust method to locate epochs than index)
    if verbose:
        print("Getting weights...")
    m1.get_weights(m1.epochs[1][~m1.ripples])

    if len(m1.patterns) == 0:
        print("No patterns found")
        return None

    # get assembly activity
    if verbose:
        print("Getting assembly activity...")
    assembly_act_pre = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[0]])
    assembly_act_task = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[1]])
    assembly_act_post = m1.get_assembly_act(epoch=m1.ripples[m1.epochs[2]])
    results = {
        "assembly_act_pre": assembly_act_pre,
        "assembly_act_task": assembly_act_task,
        "assembly_act_post": assembly_act_post,
        "react": m1,
    }

    return results


def compile_results_df(results):

    patterns, is_member_sig, keep_assembly, is_member = find_sig_assembly.main(
        results.get("react").patterns
    )

    assembly_df = pd.DataFrame()
    assembly_df["patterns"] = patterns.ravel()
    assembly_df["is_member_sig"] = is_member_sig.ravel()
    assembly_df["assembly_n"] = (
        (np.ones_like(patterns).T * np.arange(patterns.shape[0])).T.astype(int).ravel()
    )
    assembly_df["UID"] = np.tile(
        results.get("react").cell_metrics.UID.values, patterns.shape[0]
    )
    assembly_df["putativeCellType"] = np.tile(
        results.get("react").cell_metrics.putativeCellType.values, patterns.shape[0]
    )
    assembly_df["brainRegion"] = np.tile(
        results.get("react").cell_metrics.brainRegion.values, patterns.shape[0]
    )
    assembly_df["deepSuperficial"] = np.tile(
        results.get("react").cell_metrics.deepSuperficial.values, patterns.shape[0]
    )
    assembly_df["deepSuperficialDistance"] = np.tile(
        results.get("react").cell_metrics.deepSuperficialDistance.values,
        patterns.shape[0],
    )
    assembly_df = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(assembly_df)

    deep_mec = []
    deep_pfc = []
    superficial_mec = []
    superficial_pfc = []
    n_deep = []
    n_sup = []
    n_mec = []
    n_pfc = []
    for n in assembly_df.assembly_n.unique():
        temp_assembly_df = assembly_df[
            (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
        ]
        n_deep.append(np.sum(temp_assembly_df.deepSuperficial == "Deep"))
        n_sup.append(np.sum(temp_assembly_df.deepSuperficial == "Superficial"))
        n_mec.append(np.sum(temp_assembly_df.deepSuperficial == "EC1|EC2|EC3|EC4|EC5|MEC"))
        n_pfc.append(np.sum(temp_assembly_df.deepSuperficial == "PFC"))

        deep_mec.append(
            any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
            & any((temp_assembly_df.deepSuperficial == "Deep"))
        )
        deep_pfc.append(
            any(temp_assembly_df.brainRegion.str.contains("PFC"))
            & any((temp_assembly_df.deepSuperficial == "Deep"))
        )
        superficial_mec.append(
            any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
            & any((temp_assembly_df.deepSuperficial == "Superficial"))
        )
        superficial_pfc.append(
            any(temp_assembly_df.brainRegion.str.contains("PFC"))
            & any((temp_assembly_df.deepSuperficial == "Superficial"))
        )

    prop_df = pd.DataFrame()
    prop_df["prop_cross_region"] = [
        np.mean(np.array(superficial_pfc) > 0),
        np.mean(np.array(superficial_mec) > 0),
        np.mean(np.array(deep_pfc) > 0),
        np.mean(np.array(deep_mec) > 0),
    ]

    prop_df["n_cross_region"] = [
        np.sum(np.array(superficial_pfc) > 0),
        np.sum(np.array(superficial_mec) > 0),
        np.sum(np.array(deep_pfc) > 0),
        np.sum(np.array(deep_mec) > 0),
    ]
    
    prop_df["labels"] = ["Superficial PFC", "Superficial MEC", "Deep PFC", "Deep MEC"]

    results.get(
        "react"
    ).cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(
        results.get("react").cell_metrics
    )

    prop_df["n_deep"] = sum(results.get("react").cell_metrics.deepSuperficial == "Deep")
    prop_df["n_sup"] = sum(
        results.get("react").cell_metrics.deepSuperficial == "Superficial"
    )
    prop_df["n_mec"] = sum(
        results.get("react").cell_metrics.brainRegion.str.contains(
            "EC1|EC2|EC3|EC4|EC5|MEC"
        )
    )
    prop_df["n_pfc"] = sum(
        results.get("react").cell_metrics.brainRegion.str.contains("PFC")
    )
    prop_df["n_assemblies"] = len(assembly_df.assembly_n.unique())
    prop_df["basepath"] = results.get("react").basepath

    return prop_df, assembly_df


def load_results(save_path, verbose=False):
    """
    load_results: load results from a directory
    """
    warnings.filterwarnings("ignore")
    print("Loading results...")
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    prop_df = pd.DataFrame()
    assembly_df = pd.DataFrame()
    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        prop_df_, assembly_df_ = compile_results_df(results)

        prop_df = pd.concat([prop_df, prop_df_], ignore_index=True)
        assembly_df = pd.concat(
            [assembly_df, assembly_df_], ignore_index=True
        )
    return prop_df, assembly_df


# def run(basepath, regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC", putativeCellType="Pyr"):

#     ar = assembly_reactivation.AssemblyReact(
#         basepath, brainRegion=regions, putativeCellType=putativeCellType
#     )
#     # use built in function to load needed data
#     ar.load_data()
#     # retain pre task post task structure
#     ar.restrict_epochs_to_pre_task_post()
#     # get weights for the task
#     ar.get_weights([ar.epochs[1]])
#     # locate sig assemblies and sig members
#     patterns, is_member_sig, keep_assembly, is_member = find_sig_assembly.main(
#         ar.patterns
#     )

#     assembly_df = pd.DataFrame()
#     assembly_df["patterns"] = patterns.ravel()
#     assembly_df["is_member_sig"] = is_member_sig.ravel()
#     assembly_df["assembly_n"] = (
#         (np.ones_like(patterns).T * np.arange(patterns.shape[0])).T.astype(int).ravel()
#     )
#     assembly_df["UID"] = np.tile(ar.cell_metrics.UID.values, patterns.shape[0])
#     assembly_df["putativeCellType"] = np.tile(
#         ar.cell_metrics.putativeCellType.values, patterns.shape[0]
#     )
#     assembly_df["brainRegion"] = np.tile(
#         ar.cell_metrics.brainRegion.values, patterns.shape[0]
#     )
#     assembly_df["deepSuperficial"] = np.tile(
#         ar.cell_metrics.deepSuperficial.values, patterns.shape[0]
#     )
#     assembly_df["deepSuperficialDistance"] = np.tile(
#         ar.cell_metrics.deepSuperficialDistance.values, patterns.shape[0]
#     )


# deep_mec = []
# deep_pfc = []
# superficial_mec = []
# superficial_pfc = []

# for n in assembly_df.assembly_n.unique():
#     temp_assembly_df = assembly_df[
#         (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
#     ]
#     deep_mec.append(
#         any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
#         & any((temp_assembly_df.deepSuperficial == "Deep"))
#     )
#     deep_pfc.append(
#         any(temp_assembly_df.brainRegion.str.contains("PFC"))
#         & any((temp_assembly_df.deepSuperficial == "Deep"))
#     )
#     superficial_mec.append(
#         any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
#         & any((temp_assembly_df.deepSuperficial == "Superficial"))
#     )
#     superficial_pfc.append(
#         any(temp_assembly_df.brainRegion.str.contains("PFC"))
#         & any((temp_assembly_df.deepSuperficial == "Superficial"))
#     )


# prop_df = pd.DataFrame()
# prop_df["prop_cross_region"] = [
#     np.mean(np.array(superficial_pfc) > 0),
#     np.mean(np.array(superficial_mec) > 0),
#     np.mean(np.array(deep_pfc) > 0),
#     np.mean(np.array(deep_mec) > 0),
# ]
# prop_df["labels"] = ["Superficial PFC", "Superficial MEC", "Deep PFC", "Deep MEC"]
# prop_df["n_assemblies"] = len(superficial_pfc)
# prop_df['basepath'] = basepath

# results = {"assembly_df": assembly_df, "prop_df": prop_df}
# return results

# target_regions = ["PFC","EC1|EC2|EC3|EC4|EC5|MEC"]
# for ca1_sub in ["Deep", "Superficial"]:
#     # iterate over target regions
#     for region in target_regions:
#         for n in assembly_df.assembly_n.unique():
#             temp_assembly_df = assembly_df[
#                 (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
#             ]
