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


def get_maps_and_score(cell_id, st_run, x, y, ts, bin_width, n_shuff):

    ratemaps, occupancies = get_ratemap(
        ts, x, y, st_run.data[cell_id], bin_width=bin_width
    )

    pvals, null_ics, spatial_infos = surrogate_test_spatial_info(
        ts, x, y, st_run.data[cell_id], bin_width=bin_width, n_shuff=n_shuff
    )

    return ratemaps, occupancies, pvals, null_ics, spatial_infos


def run(
    basepath,
    n_shuff=500,
    min_session_duration=5,
    bin_width=3,
    speed_thres=4,
    epochs_to_skip=["sleep", "wheel", "cheeseboard"],
    brainRegion="CA1",
    cell_type="Pyr",
):

    st_unit, cell_metrics = loading.load_spikes(
        basepath, brainRegion=brainRegion, putativeCellType=cell_type
    )
    if cell_metrics.shape[0] == 0:
        return

    epoch_df = loading.load_epoch(basepath)

    # remove sleep and wheel running and others
    for epoch in epochs_to_skip:
        epoch_df = epoch_df[epoch_df.environment != epoch]

    # remove sessions < 5 minutes
    epoch_df = epoch_df[
        (epoch_df.stopTime - epoch_df.startTime) / 60 > min_session_duration
    ]
    # exit if no sessions left
    if epoch_df.shape[0] == 0:
        return

    beh_df = loading.load_animal_behavior(basepath)

    # exist if no behavior data
    if len(beh_df) == 0:
        return
    if (
        np.isnan(beh_df.x).all()
        & np.isnan(beh_df.y).all()
        & np.isnan(beh_df.linearized).all()
    ):
        return

    # if no x, use linearized
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
    env = []
    name = []
    startTime = []
    stopTime = []
    xs = []
    ys = []
    tss = []
    st = []
    linear_dir = []
    for ep_i, ep in enumerate(beh_epochs):

        x = pos[ep].data[0, :]

        if (len(x) == 0) | np.isnan(x).all():
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

        if st_run.n_active == 0:
            print("no spk data in epoch")
            continue

        if "linear" in epoch_df.iloc[ep_i].environment:
            (outbound_epochs, inbound_epochs) = functions.get_linear_track_lap_epochs(
                pos_run[ep].abscissa_vals, pos_run.data[0], newLapThreshold=20
            )
            for dir_epoch_i, dir_epoch in enumerate([outbound_epochs, inbound_epochs]):
                # check if no laps in this direction
                if dir_epoch.n_intervals == 0:
                    continue

                ts = pos_run[dir_epoch].abscissa_vals
                x = pos_run[dir_epoch].data[0, :]
                y = np.ones_like(x)

                for cell_id in range(st_run.data.shape[0]):

                    (
                        ratemap_,
                        occupancy_,
                        pval,
                        null_ic,
                        spatial_info,
                    ) = get_maps_and_score(
                        cell_id, st_run[dir_epoch], x, y, ts, bin_width, n_shuff
                    )

                    pvals.append(pval)
                    null_ics.append(null_ic)
                    spatial_infos.append(spatial_info)
                    UID.append(cell_metrics.iloc[cell_id].UID)
                    env.append(epoch_df.iloc[ep_i].environment)
                    name.append(epoch_df.name.iloc[ep_i])
                    startTime.append(epoch_df.iloc[ep_i].startTime)
                    stopTime.append(epoch_df.iloc[ep_i].stopTime)
                    xs.append(x)
                    ys.append(y)
                    tss.append(ts)
                    st.append(st_run[dir_epoch].data[cell_id])
                    ratemaps.append(ratemap_)
                    occupancies.append(occupancy_)
                    linear_dir.append(dir_epoch_i)
        else:
            ts = pos_run.abscissa_vals
            x = pos_run.data[0, :]
            y = pos_run.data[1, :]

            for cell_id in range(st_run.data.shape[0]):
                (
                    ratemap_,
                    occupancy_,
                    pval,
                    null_ic,
                    spatial_info,
                ) = get_maps_and_score(cell_id, st_run, x, y, ts, bin_width, n_shuff)

                pvals.append(pval)
                null_ics.append(null_ic)
                spatial_infos.append(spatial_info)
                UID.append(cell_metrics.iloc[cell_id].UID)
                env.append(epoch_df.iloc[ep_i].environment)
                name.append(epoch_df.name.iloc[ep_i])
                startTime.append(epoch_df.iloc[ep_i].startTime)
                stopTime.append(epoch_df.iloc[ep_i].stopTime)
                xs.append(x)
                ys.append(y)
                tss.append(ts)
                st.append(st_run.data[cell_id])
                ratemaps.append(ratemap_)
                occupancies.append(occupancy_)
                linear_dir.append(np.nan)

    if len(UID) == 0:
        print("no data in epoch")
        return None

    epoch["UID"] = np.hstack(UID)
    epoch["environment"] = np.hstack(env)
    epoch["name"] = np.hstack(name)
    epoch["startTime"] = np.hstack(startTime)
    epoch["stopTime"] = np.hstack(stopTime)
    epoch["spatial_infos"] = np.hstack(spatial_infos)
    epoch["pvals"] = np.hstack(pvals)
    epoch["linear_dir"] = np.hstack(linear_dir)
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
