import glob
import pickle
from ripple_heterogeneity.utils import (
    functions,
    loading,
    add_new_deep_sup,
)
import pandas as pd
import numpy as np
import nelpy as nel
import os


def get_ripple_info_df(rip_par_mat, cm):

    n_deep = []
    n_sup = []
    n_middle = []
    n_pfc = []
    n_mec = []
    n_spikes_deep = []
    n_spikes_sup = []
    n_spikes_middle = []
    n_spikes_pfc = []
    n_spikes_mec = []

    for rip in rip_par_mat.T:
        n_deep.append((cm.deepSuperficial[rip > 0] == "Deep").sum())
        n_sup.append((cm.deepSuperficial[rip > 0] == "Superficial").sum())
        n_middle.append((cm.deepSuperficial[rip > 0] == "middle").sum())
        n_pfc.append((cm.brainRegion[rip > 0] == "PFC").sum())
        n_mec.append((cm.brainRegion[rip > 0] == "MEC").sum())

        n_spikes_deep.append(rip[cm.deepSuperficial == "Deep"].sum())
        n_spikes_sup.append(rip[cm.deepSuperficial == "Superficial"].sum())
        n_spikes_middle.append(rip[cm.deepSuperficial == "middle"].sum())
        n_spikes_pfc.append(rip[cm.brainRegion == "PFC"].sum())
        n_spikes_mec.append(rip[cm.brainRegion == "MEC"].sum())

        rip_resp_df = pd.DataFrame(
            {
                "n_deep": n_deep,
                "n_sup": n_sup,
                "n_middle": n_middle,
                "n_pfc": n_pfc,
                "n_mec": n_mec,
                "n_spikes_deep": n_spikes_deep,
                "n_spikes_sup": n_spikes_sup,
                "n_spikes_middle": n_spikes_middle,
                "n_spikes_pfc": n_spikes_pfc,
                "n_spikes_mec": n_spikes_mec,
            }
        )
    # convert to int16 to save space
    columns = [
        "n_deep",
        "n_sup",
        "n_middle",
        "n_pfc",
        "n_mec",
        "n_spikes_deep",
        "n_spikes_sup",
        "n_spikes_middle",
        "n_spikes_pfc",
        "n_spikes_mec",
    ]
    rip_resp_df[columns] = rip_resp_df[columns].astype("int16")

    # calculate deep sup ratios
    rip_resp_df["deep_sup_spike_ratio"] = (
        rip_resp_df.n_spikes_deep / rip_resp_df.n_spikes_sup
    )
    rip_resp_df["deep_sup_cell_count_ratio"] = (
        rip_resp_df.n_deep - rip_resp_df.n_sup
    ) / (rip_resp_df.n_deep + rip_resp_df.n_sup)
    return rip_resp_df


def load_and_format_data(
    basepath,
    putativeCellType,
    brainRegion,
    convert_regions,
    ripple_expand,
    min_cell_per_group,
):

    st, cm = loading.load_spikes(
        basepath, putativeCellType=putativeCellType, brainRegion=brainRegion
    )
    cm = add_new_deep_sup.deep_sup_from_deepSuperficialDistance(cm)

    if ((cm.deepSuperficial == "Superficial").sum() < min_cell_per_group) | (
        (cm.deepSuperficial == "Deep").sum() < min_cell_per_group
    ):
        return None, None

    for key in convert_regions.keys():
        cm.loc[cm["brainRegion"].str.contains(key), "brainRegion"] = convert_regions[
            key
        ]

    ripples_df = loading.load_ripples_events(basepath)

    ripple_epochs = nel.EpochArray(np.array([ripples_df.start, ripples_df.stop]).T)
    # ripples were already extended in the analysis, but we want to
    # extend them again here to get xxms after the ripple
    ripple_epochs = ripple_epochs.expand(ripple_expand, direction="stop")

    rip_par_mat = functions.get_participation(
        st.data, ripple_epochs.starts, ripple_epochs.stops, par_type="counts"
    )
    return rip_par_mat, cm


def run(
    basepath,
    putativeCellType="Pyr",  # neuron type to use for analysis
    brainRegion="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC",  # regions to include
    convert_regions={
        "CA1": "CA1",
        "PFC": "PFC",
        "EC1|EC2|EC3|EC4|EC5|MEC": "MEC",
    },  # value pair to re-label regions
    ripple_expand=0.150,
    min_cell_per_ripple=1,
    min_cell_per_group=5,
):

    rip_par_mat, cm = load_and_format_data(
        basepath,
        putativeCellType,
        brainRegion,
        convert_regions,
        ripple_expand,
        min_cell_per_group,
    )
    if rip_par_mat is None:
        return None

    rip_resp_df = get_ripple_info_df(rip_par_mat, cm)

    corr_df = rip_resp_df.query(
        "n_deep>@min_cell_per_ripple & n_sup>@min_cell_per_ripple"
    ).corr()
    correlation = []
    correlation.append(corr_df["deep_sup_cell_count_ratio"]["n_pfc"])
    correlation.append(corr_df["deep_sup_cell_count_ratio"]["n_mec"])
    correlation.append(corr_df["deep_sup_cell_count_ratio"]["n_spikes_pfc"])
    correlation.append(corr_df["deep_sup_cell_count_ratio"]["n_spikes_mec"])

    correlation.append(corr_df["deep_sup_spike_ratio"]["n_pfc"])
    correlation.append(corr_df["deep_sup_spike_ratio"]["n_mec"])
    correlation.append(corr_df["deep_sup_spike_ratio"]["n_spikes_pfc"])
    correlation.append(corr_df["deep_sup_spike_ratio"]["n_spikes_mec"])

    correlation.append(corr_df["n_deep"]["n_pfc"])
    correlation.append(corr_df["n_deep"]["n_mec"])
    correlation.append(corr_df["n_deep"]["n_spikes_pfc"])
    correlation.append(corr_df["n_deep"]["n_spikes_mec"])

    correlation.append(corr_df["n_sup"]["n_pfc"])
    correlation.append(corr_df["n_sup"]["n_mec"])
    correlation.append(corr_df["n_sup"]["n_spikes_pfc"])
    correlation.append(corr_df["n_sup"]["n_spikes_mec"])

    correlation.append(corr_df["n_spikes_deep"]["n_pfc"])
    correlation.append(corr_df["n_spikes_deep"]["n_mec"])
    correlation.append(corr_df["n_spikes_deep"]["n_spikes_pfc"])
    correlation.append(corr_df["n_spikes_deep"]["n_spikes_mec"])

    correlation.append(corr_df["n_spikes_sup"]["n_pfc"])
    correlation.append(corr_df["n_spikes_sup"]["n_mec"])
    correlation.append(corr_df["n_spikes_sup"]["n_spikes_pfc"])
    correlation.append(corr_df["n_spikes_sup"]["n_spikes_mec"])


    label = (
        "cell_count_n_pfc",
        "cell_count_n_mec",
        "cell_count_n_spikes_pfc",
        "cell_count_n_spikes_mec",
        "spike_count_n_pfc",
        "spike_count_n_mec",
        "spike_count_n_spikes_pfc",
        "spike_count_n_spikes_mec",
        "n_deep_n_pfc",
        "n_deep_n_mec",
        "n_deep_n_spikes_pfc",
        "n_deep_n_spikes_mec",
        "n_sup_n_pfc",
        "n_sup_n_mec",
        "n_sup_n_spikes_pfc",
        "n_sup_n_spikes_mec",
        "n_spikes_deep_n_pfc",
        "n_spikes_deep_n_mec",
        "n_spikes_deep_n_spikes_pfc",
        "n_spikes_deep_n_spikes_mec",
        "n_spikes_sup_n_pfc",
        "n_spikes_sup_n_mec",
        "n_spikes_sup_n_spikes_pfc",
        "n_spikes_sup_n_spikes_mec",
    )
    results_df = pd.DataFrame(
        {
            "correlation": correlation,
            "label": label,
        }
    )
    results_df["n_deep"] = (cm.deepSuperficial == "Deep").sum()
    results_df["n_sup"] = (cm.deepSuperficial == "Superficial").sum()
    results_df["n_middle"] = (cm.deepSuperficial == "middle").sum()
    results_df["n_pfc"] = (cm.brainRegion == "PFC").sum()
    results_df["n_mec"] = (cm.brainRegion == "MEC").sum()
    results_df["basepath"] = basepath

    results = {"results_df": results_df, "rip_resp_df": rip_resp_df}

    return results


def load_results(save_path):
    """
    load_results: load results from a directory
    """
    sessions = glob.glob(save_path + os.sep + "*.pkl")
    results_df = pd.DataFrame()
    rip_resp_df = pd.DataFrame()
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
        if results is None:
            continue
        results_df = pd.concat([results_df, results["results_df"]], ignore_index=True)
        rip_resp_df = pd.concat([rip_resp_df, results["rip_resp_df"]], ignore_index=True)
    return results_df, rip_resp_df
