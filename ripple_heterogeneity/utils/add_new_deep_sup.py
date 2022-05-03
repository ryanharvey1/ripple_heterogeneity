import pandas as pd
import numpy as np
from ripple_heterogeneity.utils import loading


def add_new_deep_sup_class(df, layer_dist=30):
    """
    Take df dataframe and update deepSuperficial classification
    Inputs:
        df: dataframe with at least basepath and UID
        layer_dist: distance from pyramidal layer
    Outputs:
        df: dataframe with deepSuperficial classification
    """

    def assign_region(df, cell_metrics_):
        """
        Assign deepSuperficial classification based on distance from the pyramidal layer
        """
        basepath = df.basepath.iloc[0]
        cell_metrics = cell_metrics_[cell_metrics_.basepath == basepath]

        for uid in df.UID.unique():
            deepSuperficialDistance = cell_metrics[
                cell_metrics.UID == uid
            ].deepSuperficialDistance
            df.loc[df.UID == uid, "deepSuperficialDistance"] = np.tile(
                deepSuperficialDistance, sum(df.UID == uid)
            )
        return df

    cell_metrics = loading.load_all_cell_metrics(df.basepath.unique())

    df_out = pd.DataFrame()
    # iter over each unique basepath
    for basepath in df.basepath.unique():
        df_out = pd.concat(
            [df_out, assign_region(df[df.basepath == basepath].copy(), cell_metrics)],
            ignore_index=True,
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

    return df_out


def deep_sup_from_deepSuperficialDistance(cell_metrics, layer_dist=30):
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

    return cell_metrics
