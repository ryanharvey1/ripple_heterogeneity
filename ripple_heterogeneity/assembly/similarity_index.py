import itertools
import numpy as np
from ripple_heterogeneity.utils import functions


def similarity_index(patterns, n_shuffles=1000):
    """
    Calculate the similarity index of a set of patterns.

    Based on Almeida-Filho et al., 2014 to detect similar assemblies.

    To use a quantitative criterion to compare assembly composition,
    a Similarity Index (SI) was defined as the absolute value of the
    inner product between the assembly patterns (unitary vectors) of
    two given assemblies, varying from 0 to 1. Thus, if two assemblies
    attribute large weights to the same neurons, SI will be large;
    if assemblies are orthogonal, SI will be zero.

    Input:
        patterns: list of patterns (n patterns x n neurons)
        n_shuffles: number of shuffles to calculate the similarity index
    Output:
        si: similarity index: float (0-1)
        combos: list of all possible combinations of patterns
        pvalues: list of p-values for each pattern combination
    """
    # check to see if patterns are numpy arrays
    if not isinstance(patterns, np.ndarray):
        patterns = np.array(patterns)

    # check to see if patterns have more rows than columns and transpose if necessary
    # should have fewer patterns than neurons
    if len(patterns) > len(patterns[0]):
        patterns = patterns.T

    # check if all values in matrix are less than 1
    if not all(i <= 1 for i in patterns.flatten()):
        raise ValueError("All values in matrix must be less than 1")

    # shuffle patterns over neurons
    def shuffle_patterns(patterns):
        shuffled_patterns = []
        for pattern in patterns:
            shuffled_patterns.append(np.random.permutation(pattern))
        return np.array(shuffled_patterns)

    # calculate absolute inner product between patterns
    def get_si(patterns):
        x = np.arange(0, patterns.shape[0])
        # use itertools to get all combinations of patterns
        combos = np.array(list(itertools.combinations(x, 2)))
        si = []
        for s in combos:
            si.append(np.abs(np.inner(patterns[s[0], :], patterns[s[1], :])))
        return si, combos

    # calculate observed si
    si, combos = get_si(patterns)

    # shuffle patterns and calculate si
    si_shuffles = []
    for _ in range(n_shuffles):
        si_shuffled, _ = get_si(shuffle_patterns(patterns))
        si_shuffles.append(si_shuffled)

    # calculate p-values for each pattern combination
    _, pvalues = functions.get_significant_events(np.array(si), np.array(si_shuffles))

    return si, combos, pvalues