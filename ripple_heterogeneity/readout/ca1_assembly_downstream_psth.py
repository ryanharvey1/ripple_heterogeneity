import glob
import os
import pickle
import numpy as np
import pandas as pd
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly
from ripple_heterogeneity.utils import add_new_deep_sup
from ripple_heterogeneity.utils import functions, loading, add_new_deep_sup
import warnings

warnings.filterwarnings("ignore")


def get_assembly_df(results, patterns, is_member):
    """
    Returns a dataframe of the assembly properties.
    """
    assembly_df = pd.DataFrame()

    # patterns = results.get("react").patterns

    assembly_df["patterns"] = patterns.ravel()
    assembly_df["is_member"] = is_member.ravel()
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
    return assembly_df


def get_strength_matrix(analog_signal):
    """
    Returns a matrix of the strength of the signal at each timepoint.
    input:
        analog_signal: a nelpy AnalogSignalArray of the assembly strengths, must have segments
    output:
        strength_matrix: a numpy array of the strength of the signal at each timepoint
    """
    strength = np.zeros([analog_signal.n_signals, analog_signal.n_intervals])
    time = strength.copy()
    for rip_i, rip in enumerate(analog_signal):
        strength[:, rip_i] = rip.max()
        time[:, rip_i] = rip.abscissa_vals[np.argmax(rip.data, axis=1)]
    return strength, time


def get_psths(
    st,
    assembly_strengths,
    assembly_df,
    labels=["Deep", "Superficial", "mixed"],
    assembly_strength_thres=5,
    bin_width=0.005,
    n_bins=200,
):
    """
    Returns a dict of the PSTH of st relative to when the assembly is active.
    input:
        st: a nelpy SpikeTrainArray of the st of the assembly
        assembly_strengths: a nelpy analogsignal array of the strength of the signal at each ripple
        assembly_df: a dataframe of the assembly properties
        labels: a list of the labels of the assemblies
        assembly_strength_thres: the threshold for the strength of the assembly
        bin_width: the width of the bins in seconds
        n_bins: the number of bins
    output:
        psths: a dict of the PSTH of st relative to when the assembly is active
    """
    psths = {}
    for label in labels:
        psth = pd.DataFrame()
        # get the max value and time of max value for each assembly over ripples
        strength, time = get_strength_matrix(assembly_strengths)

        # keep the assemblies that are specific to the label
        assembly_n = assembly_df[
            (assembly_df.assembly_label == label) & (assembly_df.assembly_keep == True)
        ].assembly_n.unique()

        strength = strength[assembly_n, :]
        time = time[assembly_n, :]

        # iterate over the assemblies and calculate the psth relative to the time of max strength
        for i in range(strength.shape[0]):
            psth = pd.concat(
                [
                    psth,
                    functions.compute_psth(
                        st.data,
                        time[i, (strength[i, :] > assembly_strength_thres)],
                        bin_width=bin_width,
                        n_bins=n_bins,
                    ),
                ],
                ignore_index=True,
                axis=1,
            )
        # add the psth to the dict
        psths[label] = psth
    return psths


def run(
    basepath,
    regions="CA1",
    target_regions=["PFC", "EC1|EC2|EC3|EC4|EC5|MEC"],
    putativeCellType="Pyr",
    weight_dt=0.02,  # dt in seconds for binning st to get weights for each assembly
    verbose=False,  # print out progress
    rip_exp_start=0.05,  # ripple expansion start, in seconds, how much to expand ripples
    rip_exp_stop=0.05,  # ripple expansion stop, in seconds, how much to expand ripples
    min_cells=5,  # minimum number of cells in analysis (n deep, n superficial, n cortex)
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

    # check if enough data is available
    m1.cell_metrics = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(
        m1.cell_metrics
    )

    if not ((m1.cell_metrics.deepSuperficial == "Deep").sum() >= min_cells) | (
        (m1.cell_metrics.deepSuperficial == "Superficial").sum() >= min_cells
    ):
        print("Not enough cells")
        return None
    
    # if not ((m1.cell_metrics.deepSuperficial == "Deep").sum() >= min_cells) & (
    #     (m1.cell_metrics.deepSuperficial == "Superficial").sum() >= min_cells
    # ):
    #     print("Not enough cells")
    #     return None

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

    # get activity for target regions
    st, cell_metrics = loading.load_spikes(
        basepath, brainRegion=target_regions, putativeCellType=putativeCellType
    )

    target_regions_1,target_regions_2 = target_regions 
    if not cell_metrics.brainRegion.str.contains(target_regions_1).any() | cell_metrics.brainRegion.str.contains(target_regions_2).any():
        print("Not enough cells")
        return None
    
    _, _, keep_assembly, is_member = find_sig_assembly.main(
        results.get("react").patterns
    )
    assembly_df = get_assembly_df(results, results.get("react").patterns, is_member)

    # add keep assembly labels to the assembly df
    keep_assembly_df = pd.DataFrame(
        index=np.arange(results.get("react").patterns.shape[0]), data=keep_assembly
    )
    keep_assembly_dict = keep_assembly_df.to_dict()
    assembly_df["assembly_keep"] = assembly_df["assembly_n"]
    assembly_df["assembly_keep"] = assembly_df["assembly_keep"].map(
        keep_assembly_dict[0]
    )

    # iterate over unique assemblies
    for assembly_n in assembly_df.assembly_n.unique():
        assem_idx = assembly_df.assembly_n == assembly_n
        sig_mem_idx = assembly_df.is_member == True

        # check if the significant members are deep or superficial or mixed
        if (assembly_df[assem_idx & sig_mem_idx].deepSuperficial == "Deep").all():
            assembly_df.loc[assem_idx, "assembly_label"] = "Deep"
        elif (
            assembly_df[assem_idx & sig_mem_idx].deepSuperficial == "Superficial"
        ).all():
            assembly_df.loc[assem_idx, "assembly_label"] = "Superficial"
        else:
            assembly_df.loc[assem_idx, "assembly_label"] = "mixed"

    result = {}
    result["pre"] = get_psths(st, results.get("assembly_act_pre"), assembly_df)

    result["task"] = get_psths(st, results.get("assembly_act_task"), assembly_df)

    result["post"] = get_psths(st, results.get("assembly_act_post"), assembly_df)

    result["cell_metrics"] = cell_metrics

    return result


def load_results(save_path, verbose=False):
    """
    load_results: load results from a pickle file
    """
    warnings.filterwarnings("ignore")
    print("Loading results...")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    pre_deep = pd.DataFrame()
    task_deep = pd.DataFrame()
    post_deep = pd.DataFrame()
    pre_sup = pd.DataFrame()
    task_sup = pd.DataFrame()
    post_sup = pd.DataFrame()
    meta_data_deep = pd.DataFrame()
    meta_data_sup = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue

        if not (
            results["pre"]["Superficial"].shape[1]
            == results["task"]["Superficial"].shape[1]
            == results["post"]["Superficial"].shape[1]
        ):
            print("Error: unequal number of columns")
            continue

        if not (
            results["pre"]["Deep"].shape[1]
            == results["task"]["Deep"].shape[1]
            == results["post"]["Deep"].shape[1]
        ):
            print("Error: unequal number of columns")
            continue

        pre_deep = pd.concat(
            [pre_deep, results["pre"]["Deep"]], axis=1, ignore_index=True
        )
        task_deep = pd.concat(
            [task_deep, results["task"]["Deep"]], axis=1, ignore_index=True
        )
        post_deep = pd.concat(
            [post_deep, results["post"]["Deep"]], axis=1, ignore_index=True
        )

        pre_sup = pd.concat(
            [pre_sup, results["pre"]["Superficial"]], axis=1, ignore_index=True
        )
        task_sup = pd.concat(
            [task_sup, results["task"]["Superficial"]], axis=1, ignore_index=True
        )
        post_sup = pd.concat(
            [post_sup, results["post"]["Superficial"]], axis=1, ignore_index=True
        )

        fields_to_keep = ["UID", "brainRegion", "putativeCellType", "basepath"]

        if results["cell_metrics"].shape[0] > 0:
            n_assembly = int(
                results["pre"]["Deep"].shape[1] / results["cell_metrics"].shape[0]
            )
            for n_assem in range(n_assembly):
                temp_df = results["cell_metrics"][fields_to_keep].copy()
                temp_df["assembly_n"] = n_assem
                meta_data_deep = pd.concat([meta_data_deep, temp_df], ignore_index=True)

            if not (meta_data_deep.shape[0] == post_deep.shape[1]):
                print("Error: unequal number of columns")

            n_assembly = int(
                results["pre"]["Superficial"].shape[1]
                / results["cell_metrics"].shape[0]
            )
            for n_assem in range(n_assembly):
                temp_df = results["cell_metrics"][fields_to_keep].copy()
                temp_df["assembly_n"] = n_assem
                meta_data_sup = pd.concat([meta_data_sup, temp_df], ignore_index=True)
        else:
            print("No cell metrics found")
    return (
        pre_deep,
        task_deep,
        post_deep,
        pre_sup,
        task_sup,
        post_sup,
        meta_data_deep,
        meta_data_sup,
    )
