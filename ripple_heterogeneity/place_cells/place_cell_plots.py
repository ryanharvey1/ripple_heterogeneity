import os
import pickle
import statistics
from matplotlib import pyplot as plt
import numpy as np
from ripple_heterogeneity.utils import loading, functions

def run(basepath, unit_id, save_path):
    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )

    beh_df = loading.load_animal_behavior(basepath)
    
    # interp over nans
    beh_df.x = beh_df.x.interpolate(
        method="linear",
        limit=int(1 / statistics.mode(np.diff(beh_df.time))) * 5,
    )
    beh_df.y = beh_df.y.interpolate(
        method="linear",
        limit=int(1 / statistics.mode(np.diff(beh_df.time))) * 5,
    )

    epoch_df = loading.load_epoch(basepath)
    # remove sleep and wheel running
    epoch_df = epoch_df[
        (epoch_df.environment != "sleep") & (epoch_df.environment != "wheel")
    ]
    # remove sessions < 5 minutes
    epoch_df = epoch_df[(epoch_df.stopTime - epoch_df.startTime) / 60 > 5]

    with open(save_file, "rb") as f:
        result = pickle.load(f)

    ratemaps = np.array(result["ratemaps"])[result["df"]["UID"] == unit_id]
    occ = np.array(result["occupancies"])[result["df"]["UID"] == unit_id]

    # x = np.array(result['x'])[result['df']['UID'] == unit_id]
    # y = np.array(result['y'])[result['df']['UID'] == unit_id]
    name = result["df"].name.values[result["df"]["UID"] == unit_id]
    st = np.array(result["st"])[result["df"]["UID"] == unit_id]

    # n_panels = int(np.ceil(len(ratemaps)/2))
    n_panels = len(ratemaps)
    fig, axs = plt.subplots(
        2,
        n_panels,
        figsize=functions.set_size("thesis", fraction=1.25, subplots=(3, n_panels)),
        edgecolor="k",
    )
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    axs = axs.ravel()

    max_rate = [np.max(r) for r in ratemaps]
    v_max = np.min(max_rate)

    for i in range(len(ratemaps)):
        # axs[i].imshow(ratemaps[i])
        # plt.plot(x[i],y[i])

        ts = beh_df[
            beh_df["time"].between(
                epoch_df.iloc[i].startTime, epoch_df.iloc[i].stopTime
            )
        ].time
        x1 = beh_df[
            beh_df["time"].between(
                epoch_df.iloc[i].startTime, epoch_df.iloc[i].stopTime
            )
        ].x
        y1 = beh_df[
            beh_df["time"].between(
                epoch_df.iloc[i].startTime, epoch_df.iloc[i].stopTime
            )
        ].y

        if (len(x1) == 0) | (len(st[i]) == 0):
            test = 1

        axs[i].plot(x1, y1, color="grey", alpha=0.5)
        axs[i].plot(
            np.interp(st[i], ts, x1), np.interp(st[i], ts, y1), ".k", alpha=0.25
        )
        axs[i].axis("equal")
        axs[i].axis("off")

        axs[i].set_title("epoch " + str(i) + "\n" + name[i], fontsize=7)
        # axs[i].show()
        # ratemap_,_ = get_ratemap(ts,x1,y1,st[i],bin_width=3,smooth_sigma=1,add_nan_back=True)

        ratemap_ = ratemaps[i].copy()
        ratemap_[occ[i] < 0.01] = np.nan
        axs[i + len(ratemaps)].imshow(
            ratemap_,
            interpolation="nearest",
            origin="lower",
            vmax=np.nanmax(ratemap_) * 0.8,
        )
        axs[i + len(ratemaps)].axis("off")

        # sns.heatmap(ratemaps[i],ax=axs[i+len(ratemaps)])
        axs[i + len(ratemaps)].axis("equal")

    return axs,fig