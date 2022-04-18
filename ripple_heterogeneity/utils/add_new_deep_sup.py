import pandas as pd
import numpy as np
from ripple_heterogeneity.utils import loading

def add_new_deep_sup_class(df, layer_dist=30):
    """
    Take df dataframe and update deepSuperficial classification

    """
    LAYERDIST = layer_dist

    def assign_region(df,cell_metrics_):
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

        deep = -LAYERDIST
        middle = [-LAYERDIST, LAYERDIST]
        sup = LAYERDIST
        df.loc[df.deepSuperficialDistance <= deep, "deepSuperficial"] = "Deep"
        df.loc[
            (df.deepSuperficialDistance > middle[0])
            & (df.deepSuperficialDistance < middle[1]),
            "deepSuperficial",
        ] = "middle"
        df.loc[df.deepSuperficialDistance >= sup, "deepSuperficial"] = "Superficial"
        return df

    cell_metrics = loading.load_all_cell_metrics(df.basepath.unique())

    df_temp = pd.DataFrame()
    # iter over each unique basepath
    for basepath in df.basepath.unique():
        df_temp = pd.concat(
            [df_temp, assign_region(df[df.basepath == basepath].copy(),cell_metrics)],
            ignore_index=True,
        )

    return df_temp