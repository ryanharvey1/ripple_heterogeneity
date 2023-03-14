from itertools import chain
from scipy import stats
import numpy as np
import nelpy as nel
from ripple_heterogeneity.assembly import assembly_reactivation

def test_assembly_reactivation():

    def lif_neuron(n_steps=1000, alpha=0.01, rate=10):
        """ Simulate a linear integrate-and-fire neuron.

        Args:
        n_steps (int): The number of time steps to simulate the neuron's activity.
        alpha (float): The input scaling factor
        rate (int): The mean rate of incoming spikes

        """
        # Precompute Poisson samples for speed
        exc = stats.poisson(rate).rvs(n_steps)

        # Initialize voltage and spike storage
        v = np.zeros(n_steps)
        spike_times = []

        # Loop over steps
        for i in range(1, n_steps):

            # Update v
            dv = alpha * exc[i]
            v[i] = v[i-1] + dv

            # If spike happens, reset voltage and record
            if v[i] > 1:
                spike_times.append(i)
                v[i] = 0

        return v, spike_times

    def jitter(spike_times, rate=10):
        return [
            spike_times + stats.poisson(rate).rvs(len(spike_times)) * 0.01
            for _ in np.arange(0, 0.1, 0.01)
        ]

    def generate_data():
        v, spike_times_ = lif_neuron()
        spike_times_1 = jitter(spike_times_, rate=10)

        v, spike_times_ = lif_neuron()
        spike_times_2 = jitter(spike_times_, rate=10)

        v, spike_times_ = lif_neuron()
        spike_times_3 = jitter(spike_times_, rate=10)

        v, spike_times_ = lif_neuron()
        spike_times_4 = jitter(spike_times_, rate=10)

        v, spike_times_ = lif_neuron()
        spike_times_5 = jitter(spike_times_, rate=10)

        spike_times = np.array(
            list(
                chain.from_iterable(
                    [spike_times_1, spike_times_2, spike_times_3, spike_times_4, spike_times_5]
                )
            ),
            dtype=object,
        )
        return spike_times
    
    # generate spike trains    
    st = nel.SpikeTrainArray(timestamps=generate_data())

    # create assembly reactivation object
    assembly_react = assembly_reactivation.AssemblyReact(
        )
    
    # test is empty
    assert assembly_react.isempty == True

    # load spike trains
    assembly_react.add_st(st)

    # test is empty
    assert assembly_react.isempty == False

    # detect assemblies
    assembly_react.get_weights()

    # test number of assemblies
    assert assembly_react.n_assemblies() == 5

    assembly_react.find_members()
    n_members_per_assembly = [np.sum(assembly_react.assembly_members[i]) for i in range(0,assembly_react.n_assemblies())]
    assert n_members_per_assembly == [10, 10, 10, 10, 10]