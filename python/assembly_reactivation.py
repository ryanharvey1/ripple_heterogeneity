import multiprocessing
from joblib import Parallel, delayed
import os
import pandas as pd
import numpy as np
import nelpy as nel
import pickle
import assembly_run
import loading
import sys

sys.path.append(r"D:\github\neurocode\reactivation\assemblies")
import assembly


class AssemblyReact(object):
    """
    Class for running assembly reactivation analysis

    Parameters:
    -----------
    basepath: str
        Path to the session folder
    brainRegion: str
        Brain region to restrict to
    putativeCellType: str
        Cell type to restrict to
    weight_dt: float
        Time resolution of the weight matrix
    z_mat_dt: float
        Time resolution of the z matrix

    attributes:
    -----------
    st: spike train
    cell_metrics: cell metrics
    ripples: ripples
    zactmat: z matrix
    ts: timestamps
    patterns: patterns
    assembly_act: assembly activity

    methods:
    --------
    load_data: load data
    restrict_to_epoch: restrict to epoch
    restrict_to_events: restrict to events
    get_z_mat: get z matrix
    get_weights: get weights
    get_assembly_act: get assembly activity
    session_loop_activation: session loop activation

    """

    def __init__(
        self,
        basepath,
        brainRegion="CA1",
        putativeCellType="Pyramidal Cell",
        weight_dt=0.01,
        z_mat_dt=0.002,
    ):
        self.basepath = basepath
        self.brainRegion = brainRegion
        self.putativeCellType = putativeCellType
        self.weight_dt = weight_dt
        self.z_mat_dt = z_mat_dt

    def load_data(self):
        self.st, self.cell_metrics = loading.load_spikes(
            self.basepath,
            brainRegion=self.brainRegion,
            putativeCellType=self.putativeCellType,
        )
        ripples = loading.load_ripples_events(self.basepath)
        self.ripples = nel.EpochArray([np.array([ripples.start, ripples.stop]).T])

        epoch_df = loading.load_epoch(self.basepath)
        self.epoch_df = nel.EpochArray(
            [np.array([epoch_df.startTime, epoch_df.stopTime]).T],
            label=epoch_df.environment.values,
        )

    def add_st(self, st):
        self.st = st
    
    def add_ripples(self, ripples):
        self.ripples = ripples
    
    def add_epoch_df(self, epoch_df):
        self.epoch_df = epoch_df

    def restrict_to_epoch(self, epoch):
        self.st = self.st[epoch]

    def restrict_to_events(self, events):
        self.st = self.st[events]

    def get_z_mat(self):
        self.zactmat, self.ts = assembly_run.get_z_t(self.st, ds=self.z_mat_dt)

    def get_weights(self):
        (
            self.patterns,
            _,
            _,
        ) = assembly.runPatterns(self.st.bin(ds=self.weight_dt).data)

    def get_assembly_act(self):
        self.assembly_act = nel.AnalogSignalArray(
            data=assembly.computeAssemblyActivity(self.patterns, self.zactmat),
            timestamps=self.ts,
            fs=1 / self.z_mat_dt,
        )

    def session_loop_activation(self):
        pass
