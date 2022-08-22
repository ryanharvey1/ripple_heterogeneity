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


def get_ripple_info_df(rip_par_mat, cm, ripple_epochs):

    n_deep = []
    n_sup = []
    n_pfc = []
    n_mec = []

    n_spikes_deep = []
    n_spikes_sup = []
    n_spikes_pfc = []
    n_spikes_mec = []

    pop_rate_deep = []
    pop_rate_sup = []
    pop_rate_pfc = []
    pop_rate_mec = []

    avg_pop_rate_deep = []
    avg_pop_rate_sup = []
    avg_pop_rate_pfc = []
    avg_pop_rate_mec = []

    for rip, duration in zip(rip_par_mat.T, ripple_epochs.lengths):
        n_deep.append((cm.deepSuperficial[rip > 0] == "Deep").sum())
        n_sup.append((cm.deepSuperficial[rip > 0] == "Superficial").sum())
        n_pfc.append((cm.brainRegion[rip > 0] == "PFC").sum())
        n_mec.append((cm.brainRegion[rip > 0] == "MEC").sum())

        n_spikes_deep.append(rip[cm.deepSuperficial == "Deep"].sum())
        n_spikes_sup.append(rip[cm.deepSuperficial == "Superficial"].sum())
        n_spikes_pfc.append(rip[cm.brainRegion == "PFC"].sum())
        n_spikes_mec.append(rip[cm.brainRegion == "MEC"].sum())

        pop_rate_deep.append(rip[cm.deepSuperficial == "Deep"].sum() / duration)
        pop_rate_sup.append(rip[cm.deepSuperficial == "Superficial"].sum() / duration)
        pop_rate_pfc.append(rip[cm.brainRegion == "PFC"].sum() / duration)
        pop_rate_mec.append(rip[cm.brainRegion == "MEC"].sum() / duration)

        avg_pop_rate_deep.append((rip[cm.deepSuperficial == "Deep"] / duration).mean())
        avg_pop_rate_sup.append((rip[cm.deepSuperficial == "Superficial"] / duration).mean())
        avg_pop_rate_pfc.append((rip[cm.deepSuperficial == "PFC"] / duration).mean())
        avg_pop_rate_mec.append((rip[cm.deepSuperficial == "MEC"] / duration).mean())

    rip_resp_df = pd.DataFrame(
        {
            "n_deep": n_deep,
            "n_sup": n_sup,
            "n_pfc": n_pfc,
            "n_mec": n_mec,
            "n_spikes_deep": n_spikes_deep,
            "n_spikes_sup": n_spikes_sup,
            "n_spikes_pfc": n_spikes_pfc,
            "n_spikes_mec": n_spikes_mec,
            "pop_rate_deep": pop_rate_deep,
            "pop_rate_sup": pop_rate_sup,
            "pop_rate_pfc": pop_rate_pfc,
            "pop_rate_mec": pop_rate_mec,
            'avg_pop_rate_deep':avg_pop_rate_deep,
            'avg_pop_rate_sup':avg_pop_rate_sup,
            'avg_pop_rate_pfc':avg_pop_rate_pfc,
            'avg_pop_rate_mec':avg_pop_rate_mec
        }
    )
    # convert to int16 to save space
    columns = [
        "n_deep",
        "n_sup",
        "n_pfc",
        "n_mec",
        "n_spikes_deep",
        "n_spikes_sup",
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
        return None, None, None

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
    return rip_par_mat, cm, ripple_epochs


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

    rip_par_mat, cm, ripple_epochs = load_and_format_data(
        basepath,
        putativeCellType,
        brainRegion,
        convert_regions,
        ripple_expand,
        min_cell_per_group,
    )
    if rip_par_mat is None:
        return None

    rip_resp_df = get_ripple_info_df(rip_par_mat, cm, ripple_epochs)

    corr_df = rip_resp_df.query(
        "n_deep>@min_cell_per_ripple & n_sup>@min_cell_per_ripple"
    ).corr()

    # take upper triangle of corr matrix
    corr_df = (
        corr_df.where(np.triu(np.ones(corr_df.shape)).astype(np.bool))
        .stack()
        .reset_index()
        .rename(columns={0: "correlation"})
    )
    # remove diag
    corr_df = corr_df.drop(corr_df[corr_df.level_0 == corr_df.level_1].index)

    # re-label comparisons
    corr_df["label"] = corr_df.level_0 + "_" + corr_df.level_1
    corr_df = corr_df.drop(['level_0', 'level_1'], axis=1)

    corr_df["n_deep"] = (cm.deepSuperficial == "Deep").sum()
    corr_df["n_sup"] = (cm.deepSuperficial == "Superficial").sum()
    corr_df["n_pfc"] = (cm.brainRegion == "PFC").sum()
    corr_df["n_mec"] = (cm.brainRegion == "MEC").sum()
    corr_df["basepath"] = basepath

    results = {"results_df": corr_df, "rip_resp_df": rip_resp_df}

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
        results["rip_resp_df"]["basepath"] = results["results_df"]["basepath"].iloc[0]
        rip_resp_df = pd.concat(
            [rip_resp_df, results["rip_resp_df"]], ignore_index=True
        )
    return results_df, rip_resp_df
