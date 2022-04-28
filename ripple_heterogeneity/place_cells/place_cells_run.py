import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
from ripple_heterogeneity.utils import functions, loading
import statistics
import nelpy as nel
from scipy.ndimage import gaussian_filter
from scipy.ndimage import rotate


def get_ratemap(ts, x, y, st, bin_width=3, smooth_sigma=1, add_nan_back=False):

    fs = 1 / statistics.mode(np.diff(ts))

    x_edges = np.arange(np.nanmin(x), np.nanmax(x), bin_width)
    y_edges = np.arange(np.nanmin(y), np.nanmax(y), bin_width)

    if len(y_edges) == 0:
        y_edges = 1

    occ, _, _ = np.histogram2d(x, y, bins=(x_edges, y_edges))
    occ = occ / fs

    spk_mat, _, _ = np.histogram2d(
        np.interp(st, ts, x), np.interp(st, ts, y), bins=(x_edges, y_edges)
    )

    ratemap = spk_mat / occ
    bad_idx = np.isnan(ratemap) | np.isinf(ratemap)
    ratemap[bad_idx] = 0

    ratemap = gaussian_filter(ratemap, sigma=smooth_sigma)

    if add_nan_back:
        ratemap[bad_idx] = np.nan

    ratemap = rotate(ratemap, angle=-90)
    occ = rotate(occ, angle=-90)

    ratemap = np.fliplr(ratemap)
    occ = np.fliplr(occ)

    return ratemap, occ


def surrogate_test_spatial_info(ts, x, y, s, n_shuff=500, bin_width=3):
    def shuff_coords(x, y, n_shuff):
        range_ = len(x)

        if range_ * 2 < n_shuff:
            n_shuff = range_ * 2

        surrogate = np.random.choice(
            np.arange(-range_, range_), size=n_shuff, replace=False
        )
        x_temp = []
        y_temp = []
        for n in surrogate:
            x_temp.append(np.roll(x, n))
            y_temp.append(np.roll(y, n))
        return x_temp, y_temp

    def pvalue(shuff_dist, score):
        # DOI: 10.2202/1544-6115.1585
        return (sum(np.abs(shuff_dist) > np.abs(score)) + 1) / (len(shuff_dist) + 1)

    ratemap, occupancy = get_ratemap(ts, x, y, s, bin_width=bin_width)
    obs_ic = functions.spatial_information(ratemap, occupancy)

    x_temp, y_temp = shuff_coords(x, y, n_shuff)

    null_ic = []

    for new_xy in zip(x_temp, y_temp):
        ratemap, occupancy = get_ratemap(
            ts, new_xy[0], new_xy[1], s, bin_width=bin_width
        )
        null_ic.append(functions.spatial_information(ratemap, occupancy))

    return pvalue(null_ic, obs_ic), null_ic, obs_ic


def run(basepath, n_shuff=500, min_session_duration=5, bin_width=3, speed_thres=4):

    st_unit, cell_metrics = loading.load_spikes(
        basepath, brainRegion="CA1", putativeCellType="Pyr"
    )

    beh_df = loading.load_animal_behavior(basepath)

    epoch_df = loading.load_epoch(basepath)
    # remove sleep and wheel running
    epoch_df = epoch_df[
        (epoch_df.environment != "sleep") & (epoch_df.environment != "wheel")
    ]
    # remove sessions < 5 minutes
    epoch_df = epoch_df[
        (epoch_df.stopTime - epoch_df.startTime) / 60 > min_session_duration
    ]

    if len(beh_df) == 0:
        print("no beh data")
        return

    if np.isnan(beh_df.x).all():
        beh_df.x = beh_df.linearized
        beh_df.y = np.zeros_like(beh_df.x)

    # make linear track linear
    for ep in epoch_df.itertuples():
        if "linear" in ep.environment:
            x = beh_df[beh_df["time"].between(ep.startTime, ep.stopTime)].x
            y = beh_df[beh_df["time"].between(ep.startTime, ep.stopTime)].y

            if np.isnan(x).sum() == len(x):
                continue

            x, y = functions.linearize_position(x, y)

            beh_df.loc[beh_df["time"].between(ep.startTime, ep.stopTime), "x"] = x
            beh_df.loc[beh_df["time"].between(ep.startTime, ep.stopTime), "y"] = y

    beh_epochs = nel.EpochArray([np.array([epoch_df.startTime, epoch_df.stopTime]).T])

    # interp over nans
    beh_df.x = beh_df.x.interpolate(
        method="linear",
        limit=int(1 / statistics.mode(np.diff(beh_df.time))) * 5,
    )
    beh_df.y = beh_df.y.interpolate(
        method="linear",
        limit=int(1 / statistics.mode(np.diff(beh_df.time))) * 5,
    )

    pos = nel.AnalogSignalArray(
        data=np.array([beh_df.x, beh_df.y]),
        timestamps=beh_df.time,
        fs=1 / statistics.mode(np.diff(beh_df.time)),
    )

    # iter over epochs
    spatial_infos = []
    pvals = []
    ratemaps = []
    occupancies = []
    null_ics = []
    UID = []
    epoch = pd.DataFrame()
    xs = []
    ys = []
    tss = []
    st = []
    for ep_i, ep in enumerate(beh_epochs):

        x = pos[ep].data[0, :]

        if len(x) == 0:
            print("no beh data in epoch")
            continue

        if (np.nanmax(x) - np.nanmin(x)) < 2:
            bin_width = 0.03
            speed_thres = 0.04
        else:
            bin_width = 3
            speed_thres = 4

        # get speed
        speed = nel.utils.ddt_asa(pos[ep], smooth=True, sigma=0.1, norm=True)
        run_epochs = nel.utils.get_run_epochs(speed, v1=speed_thres, v2=speed_thres)
        # limit units/pos by epoch and by speed
        st_run = st_unit[ep][run_epochs]
        pos_run = pos[ep][run_epochs]

        ts = pos_run.abscissa_vals
        x = pos_run.data[0, :]
        y = pos_run.data[1, :]

        if st_run.n_active == 0:
            print("no spk data in epoch")
            continue

        # iter over cells
        for cell_id in range(st_run.data.shape[0]):

            ratemap, occupancy = get_ratemap(
                ts, x, y, st_run.data[cell_id], bin_width=bin_width
            )
            ratemaps.append(ratemap)
            occupancies.append(occupancy)

            pval, null_ic, spatial_info = surrogate_test_spatial_info(
                ts, x, y, st_run.data[cell_id], bin_width=bin_width, n_shuff=n_shuff
            )
            pvals.append(pval)
            null_ics.append(null_ic)
            spatial_infos.append(spatial_info)

            UID.append(cell_metrics.UID.iloc[cell_id])
            epoch = epoch.append(epoch_df.iloc[ep_i], ignore_index=True)
            xs.append(x)
            ys.append(y)
            tss.append(ts)
            st.append(st_run.data[cell_id])

    epoch["UID"] = UID
    epoch["spatial_infos"] = spatial_infos
    epoch["pvals"] = pvals
    epoch["basepath"] = basepath

    results = {}
    results["df"] = epoch
    results["ratemaps"] = ratemaps
    results["occupancies"] = occupancies
    results["x"] = xs
    results["y"] = ys
    results["ts"] = tss
    results["st"] = st

    return results
