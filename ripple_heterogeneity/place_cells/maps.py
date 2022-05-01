import numpy as np
import nelpy as nel


class SpatialMap(object):
    def __init__(
        self,
        pos,
        st,
        dim=None,
        dir_epoch=None,
        speed_thres=4,
        ds_bst=0.05,
        s_binsize=3,
        tuning_curve_sigma=3,
    ):
        self.pos = pos
        self.st = st
        self.dim = dim
        self.dir_epoch = dir_epoch
        self.speed_thres = speed_thres
        self.ds_bst = ds_bst
        self.s_binsize = s_binsize
        self.tuning_curve_sigma = tuning_curve_sigma
        self.speed = nel.utils.ddt_asa(self.pos, smooth=True, sigma=0.1, norm=True)
        self.run_epochs = nel.utils.get_run_epochs(
            self.speed, v1=self.speed_thres, v2=self.speed_thres
        )

        if dim == 2:
            self.map_2d()
        elif dim == 1:
            self.map_1d()
        else:
            raise ValueError("dim must be 1 or 2")

    def map_1d(self):

        if self.dir_epoch is None:
            raise ValueError("dir_epoch must be specified")

        # restrict spike trains to those epochs during which the animal was running
        st_run = self.st[self.dir_epoch][self.run_epochs]

        # smooth and re-bin:
        bst_run = st_run.bin(ds=self.ds_bst)

        x_max = np.ceil(np.nanmax(self.pos[self.dir_epoch].data))
        x_min = np.floor(np.nanmin(self.pos[self.dir_epoch].data))

        n_bins = int((x_max - x_min) / self.s_binsize)

        tc = nel.TuningCurve1D(
            bst=bst_run,
            extern=self.pos[self.dir_epoch][self.run_epochs],
            n_extern=n_bins,
            extmin=x_min,
            extmax=x_max,
            sigma=self.tuning_curve_sigma,
            min_duration=0,
        )
        return tc, st_run, bst_run

    def map_2d(self):

        # restrict spike trains to those epochs during which the animal was running
        st_run = self.st[self.run_epochs]

        # smooth and re-bin:
        bst_run = st_run.bin(ds=self.ds_bst)

        x_max = np.ceil(np.nanmax(self.pos.data[0,:]))
        x_min = np.floor(np.nanmin(self.pos.data[0,:]))
        y_max = np.ceil(np.nanmax(self.pos.data[1,:]))
        y_min = np.floor(np.nanmin(self.pos.data[1,:]))

        n_bins = int((x_max - x_min) / self.s_binsize)

        tc = nel.TuningCurve2D(
            bst=bst_run,
            extern=self.pos[self.run_epochs],
            n_extern=n_bins,
            extmin=x_min,
            extmax=x_max,
            ext_xmin=x_min,
            ext_ymin=y_min,
            ext_xmax=x_max,
            ext_ymax=y_max,
            sigma=self.tuning_curve_sigma,
            min_duration=0,
        )
        return tc, st_run, bst_run
