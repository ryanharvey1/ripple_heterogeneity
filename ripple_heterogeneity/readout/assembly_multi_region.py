import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.utils import add_new_deep_sup, functions
import warnings
import logging
import nelpy as nel


def run(
    basepath,
    regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",
    putativeCellType="Pyr",
    weight_dt=0.1,  # dt in seconds for binning st to get weights for each assembly
    z_mat_dt=0.03, # dt in seconds for binning st to get activation strength
    verbose=False,  # print out progress
    rip_expand=0.05,  # expand the ripple region by this amount (not used)
    rip_exp_start=0.05,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.2,  # ripple expansion stop, in seconds, how much to expand ripples
):
    """
    Gets the pre and post assembly strengths
    """
    logging.getLogger().setLevel(logging.ERROR)

    # initialize session
    m1 = assembly_reactivation.AssemblyReact(
        basepath,
        brainRegion=regions,
        putativeCellType=putativeCellType,
        weight_dt=weight_dt,
        z_mat_dt=z_mat_dt
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
    if "FujisawaS" in basepath:
        idx = functions.find_env_paradigm_pre_task_post(m1.epoch_df)
        if len(idx) == 1:
            print("No pre/task/post epochs found")
            no_pre_task_post = True
        else:
            m1.epoch_df = m1.epoch_df[idx]
            # convert to epoch array and add to object
            m1.epochs = nel.EpochArray(
                [np.array([m1.epoch_df.startTime, m1.epoch_df.stopTime]).T]
            )
            no_pre_task_post = False
    else:
        try:
            m1.restrict_epochs_to_pre_task_post()
            no_pre_task_post = False
        except:
            print("No pre/task/post epochs found")
            no_pre_task_post = True

    # get weights for task outside ripples
    # % (TODO: use more robust method to locate epochs than index)
    if verbose:
        print("Getting weights...")

    # if no pre task post, take longest task
    if no_pre_task_post:
        epoch_df = m1.epoch_df.reset_index()
        epoch_df = epoch_df.query("environment != 'sleep'")
        epoch_df["duration"] = epoch_df.stopTime.values - epoch_df.startTime.values
        task_idx = int(epoch_df.sort_values("duration", ascending=False).index[0])
        m1.get_weights(m1.epochs[task_idx][~m1.ripples])
    else:
        m1.get_weights(m1.epochs[1][~m1.ripples])

    if len(m1.patterns) == 0:
        print("No patterns found")
        return None

    # get assembly activity
    if not no_pre_task_post:
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
    else:
        results = {
            "assembly_act_pre": None,
            "assembly_act_task": m1.get_assembly_act(
                epoch=m1.ripples[m1.epochs[task_idx]]
            ),
            "assembly_act_post": None,
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

    cm = results.get("react").cell_metrics
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

    assembly_df["UID"] = np.tile(cm.UID.values, patterns.shape[0])
    assembly_df["putativeCellType"] = np.tile(
        cm.putativeCellType.values, patterns.shape[0]
    )
    assembly_df["brainRegion"] = np.tile(cm.brainRegion.values, patterns.shape[0])
    assembly_df["deepSuperficial"] = np.tile(
        cm.deepSuperficial.values, patterns.shape[0]
    )
    assembly_df["deepSuperficialDistance"] = np.tile(
        cm.deepSuperficialDistance.values, patterns.shape[0]
    )

    deep_mec = []
    deep_pfc = []
    superficial_mec = []
    superficial_pfc = []
    n_deep = []
    n_sup = []
    n_mec = []
    n_pfc = []

    # make sure non ca1 is unknown deep/sup
    assembly_df.loc[
        ~assembly_df.brainRegion.str.contains("CA1"), "deepSuperficial"
    ] = "unknown"

    for n in assembly_df.assembly_n.unique():
        temp_assembly_df = assembly_df[
            (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
        ]
        n_deep.append(np.sum(temp_assembly_df.deepSuperficial == "Deep"))
        n_sup.append(np.sum(temp_assembly_df.deepSuperficial == "Superficial"))
        n_mec.append(
            np.sum(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
        )
        n_pfc.append(np.sum(temp_assembly_df.brainRegion == "PFC"))

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

    # make sure non ca1 is unknown deep/sup
    cm.loc[~cm.brainRegion.str.contains("CA1"), "deepSuperficial"] = "unknown"

    prop_df["n_deep"] = sum(cm.deepSuperficial == "Deep")
    prop_df["n_sup"] = sum(cm.deepSuperficial == "Superficial")
    prop_df["n_mec"] = sum(cm.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
    prop_df["n_pfc"] = sum(cm.brainRegion.str.contains("PFC"))
    prop_df["n_assemblies"] = len(assembly_df.assembly_n.unique())
    prop_df["basepath"] = results.get("react").basepath

    return prop_df, assembly_df, keep_assembly


def load_reactivation(results):
    """
    Loads reactivation data from results in pickle file.
    """

    # suppress root warnings unless severity level ERROR
    logging.getLogger().setLevel(logging.ERROR)

    def get_strength_matrix(analog_signal):
        """
        Returns a matrix of the strength of the signal at each timepoint.
        input:
            analog_signal: a nelpy AnalogSignalArray of the assembly strengths, must have segments
        output:
            strength_matrix: a numpy array of the strength of the signal at each timepoint
        """
        strength = np.zeros([analog_signal.n_signals, analog_signal.n_intervals])

        for rip_i, rip in enumerate(analog_signal):
            strength[:, rip_i] = rip.max()
        return strength

    def compile_strength_df(strength, keep_assembly):
        """
        Compiles a dataframe of the strength of the signal at each timepoint.
        input:
            strength: a numpy array of the strength of the signal at each timepoint
            keep_assembly: a list of assemblies to keep marked as significant
        output:
            strength_df: a dataframe of the strength of the signal at each timepoint
        """

        df_strength = pd.DataFrame()
        df_strength["strength"] = strength.flatten()
        # assembly index
        df_strength["assembly_n"] = (
            np.ones_like(strength).T * np.arange(strength.shape[0])
        ).T.flatten()
        # ripple index
        df_strength["ripple_n"] = (
            (np.ones_like(strength) * np.arange(strength.shape[1]))
        ).flatten()
        # if the assembly is significant, mark it as 1
        df_strength["sig"] = (np.ones_like(strength).T * keep_assembly).T.flatten()
        return df_strength

    def get_assembly_df(results, is_member):
        """
        Returns a dataframe of the assembly properties.
        """
        assembly_df = pd.DataFrame()

        patterns = results.get("react").patterns

        assembly_df["patterns"] = patterns.ravel()
        assembly_df["is_member"] = is_member.ravel()
        assembly_df["assembly_n"] = (
            (np.ones_like(patterns).T * np.arange(patterns.shape[0]))
            .T.astype(int)
            .ravel()
        )
        cm = results.get("react").cell_metrics
        cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

        assembly_df["UID"] = np.tile(cm.UID.values, patterns.shape[0])
        assembly_df["putativeCellType"] = np.tile(
            cm.putativeCellType.values, patterns.shape[0]
        )
        assembly_df["brainRegion"] = np.tile(cm.brainRegion.values, patterns.shape[0])
        assembly_df["deepSuperficial"] = np.tile(
            cm.deepSuperficial.values, patterns.shape[0]
        )
        assembly_df["deepSuperficialDistance"] = np.tile(
            cm.deepSuperficialDistance.values, patterns.shape[0]
        )

        return assembly_df

    def get_assembly_cross_members(assembly_df):
        """
        Returns lists of the number of assemblies that are cross-members of each brain region.
        """
        deep_mec = []
        deep_pfc = []
        superficial_mec = []
        superficial_pfc = []
        for n in assembly_df.assembly_n.unique():
            temp_assembly_df = assembly_df[
                (assembly_df.assembly_n == n) & (assembly_df.is_member)
            ]
            deep_mec.append(
                any(
                    temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC")
                )
                & any((temp_assembly_df.deepSuperficial == "Deep"))
            )
            deep_pfc.append(
                any(temp_assembly_df.brainRegion.str.contains("PFC"))
                & any((temp_assembly_df.deepSuperficial == "Deep"))
            )
            superficial_mec.append(
                any(
                    temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC")
                )
                & any((temp_assembly_df.deepSuperficial == "Superficial"))
            )
            superficial_pfc.append(
                any(temp_assembly_df.brainRegion.str.contains("PFC"))
                & any((temp_assembly_df.deepSuperficial == "Superficial"))
            )
        return deep_mec, deep_pfc, superficial_mec, superficial_pfc

    patterns, is_member_sig, keep_assembly, is_member = find_sig_assembly.main(
        results.get("react").patterns
    )

    strength = get_strength_matrix(results["assembly_act_pre"])
    df_strength_pre = compile_strength_df(strength, keep_assembly)
    df_strength_pre["epoch"] = "pre"

    strength = get_strength_matrix(results["assembly_act_task"])
    df_strength_task = compile_strength_df(strength, keep_assembly)
    df_strength_task["epoch"] = "task"

    strength = get_strength_matrix(results["assembly_act_post"])
    df_strength_post = compile_strength_df(strength, keep_assembly)
    df_strength_post["epoch"] = "post"

    df_strength = pd.DataFrame()
    df_strength = pd.concat(
        [df_strength_pre, df_strength_task, df_strength_post], ignore_index=True
    )

    assembly_df = get_assembly_df(results, is_member)
    deep_mec, deep_pfc, superficial_mec, superficial_pfc = get_assembly_cross_members(
        assembly_df
    )

    labels = ["deep_mec", "deep_pfc", "superficial_mec", "superficial_pfc"]
    for label in labels:
        df_strength[label] = np.zeros_like(df_strength.assembly_n)

    df_strength.loc[
        np.in1d(df_strength.assembly_n, np.where(deep_mec)[0]), "deep_mec"
    ] = 1
    df_strength.loc[
        np.in1d(df_strength.assembly_n, np.where(deep_pfc)[0]), "deep_pfc"
    ] = 1
    df_strength.loc[
        np.in1d(df_strength.assembly_n, np.where(superficial_mec)[0]), "superficial_mec"
    ] = 1
    df_strength.loc[
        np.in1d(df_strength.assembly_n, np.where(superficial_pfc)[0]), "superficial_pfc"
    ] = 1

    df_strength["basepath"] = results["react"].basepath

    return df_strength


def load_results(save_path, verbose=False):
    """
    load_results: load results from a pickle file
    """
    warnings.filterwarnings("ignore")
    print("Loading results...")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    prop_df = pd.DataFrame()
    assembly_df = pd.DataFrame()
    df_strength = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        prop_df_, assembly_df_, keep_assembly_ = compile_results_df(results)
        prop_df = pd.concat([prop_df, prop_df_], ignore_index=True)
        assembly_df = pd.concat([assembly_df, assembly_df_], ignore_index=True)

        if results["assembly_act_pre"] is not None:
            df_strength_ = load_reactivation(results)
            df_strength = pd.concat([df_strength, df_strength_], ignore_index=True)

    return prop_df, assembly_df, df_strength
