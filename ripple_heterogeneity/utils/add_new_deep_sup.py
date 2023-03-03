import pandas as pd
import numpy as np
from ripple_heterogeneity.utils import loading


def add_new_deep_sup_class(df: pd.DataFrame, layer_dist: int = 30) -> pd.DataFrame:
    """
    Take df dataframe and update deepSuperficial classification
    Inputs:
        df: dataframe with at least basepath and UID
        layer_dist: distance from pyramidal layer
    Outputs:
        df: dataframe with deepSuperficial classification
    """

    cell_metrics = loading.load_all_cell_metrics(df.basepath.unique())

    df_out = df.merge(
        cell_metrics[["basepath", "UID", "deepSuperficialDistance", "brainRegion"]],
        on=["basepath", "UID"],
        how="left",
    )

    deep = -layer_dist
    middle = [-layer_dist, layer_dist]
    sup = layer_dist
    df_out.loc[df_out.deepSuperficialDistance <= deep, "deepSuperficial"] = "Deep"

    df_out.loc[
        (df_out.deepSuperficialDistance > middle[0])
        & (df_out.deepSuperficialDistance < middle[1]),
        "deepSuperficial",
    ] = "middle"

    df_out.loc[df_out.deepSuperficialDistance >= sup, "deepSuperficial"] = "Superficial"

    # make sure these are ca1
    df_out.loc[
        ~df_out.brainRegion.str.contains("CA1|CA2"), "deepSuperficial"
    ] = "unknown"

    return df_out


def deep_sup_from_deepSuperficialDistance(
    cell_metrics: pd.DataFrame, layer_dist: int = 30, hpc_assumption: bool = True
) -> pd.DataFrame:
    """
    Assign deepSuperficial classification based on distance from the pyramidal layer
    Will work if you already have the (up to date) deepSuperficialDistance in the cell_metrics dataframe

    Input:
        cell_metrics: dataframe with deepSuperficialDistance
    Output:
        deepSuperficial: dataframe with deepSuperficial classification
    """

    deep = -layer_dist
    middle = [-layer_dist, layer_dist]
    sup = layer_dist

    cell_metrics.loc[
        cell_metrics.deepSuperficialDistance <= deep, "deepSuperficial"
    ] = "Deep"

    cell_metrics.loc[
        (cell_metrics.deepSuperficialDistance > middle[0])
        & (cell_metrics.deepSuperficialDistance < middle[1]),
        "deepSuperficial",
    ] = "middle"

    cell_metrics.loc[
        cell_metrics.deepSuperficialDistance >= sup, "deepSuperficial"
    ] = "Superficial"

    # make sure these are ca1
    if hpc_assumption:
        cell_metrics.loc[
            ~cell_metrics.brainRegion.str.contains("CA1|CA2"), "deepSuperficial"
        ] = "unknown"

    return cell_metrics


def deep_sup_from_distance(
    cell_metrics: pd.DataFrame, layer_dist: int = 30, hpc_assumption: bool = True
) -> pd.DataFrame:
    """calls deep_sup_from_deepSuperficialDistance"""
    return deep_sup_from_deepSuperficialDistance(
        cell_metrics, layer_dist=layer_dist, hpc_assumption=hpc_assumption
    )
