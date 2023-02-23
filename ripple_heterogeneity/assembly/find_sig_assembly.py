import numpy as np


def Otsu(vector):
    """
    The Otsu method for splitting data into two groups.
    This is somewhat equivalent to kmeans(vector,2), but while the kmeans implementation
    finds a local minimum and may therefore produce different results each time,
    the Otsu implementation is guaranteed to find the best division every time.

    input:
        vector: arbitrary vector
    output:
        group: binary class
        threshold: threshold used for classification
        em: effectiveness metric

    From Raly
    """
    sorted = np.sort(vector)
    n = len(vector)
    intraClassVariance = [np.nan] * n
    for i in np.arange(n):
        p = (i + 1) / n
        p0 = 1 - p
        if i + 1 == n:
            intraClassVariance[i] = np.nan
        else:
            intraClassVariance[i] = p * np.var(sorted[0 : i + 1]) + p0 * np.var(
                sorted[i + 1 :]
            )

    minIntraVariance = np.nanmin(intraClassVariance)
    idx = np.nanargmin(intraClassVariance)
    threshold = sorted[idx]
    group = vector > threshold

    em = 1 - (minIntraVariance / np.var(vector))

    return group, threshold, em


def main(patterns):
    """
    Finds significant assembly patterns and signficant assembly members

    Input:
        patterns: a list of patterns
    Output:
        patterns[keep_assembly]: a list of patterns that are significant
        is_member[keep_assembly]: a list of booleans indicating whether each pattern is a significant assembly
        keep_assembly: a list of indices of the significant assemblies
        is_member: a list of booleans indicating whether each unit is a significant member of an assembly
    """
    is_member = []
    keep_assembly = []
    for pat in patterns:
        isMember, _, _ = Otsu(np.abs(pat))
        is_member.append(isMember)

        if np.any(pat[isMember] < 0) & np.any(pat[isMember] > 0):
            keep_assembly.append(False)
        elif sum(isMember) == 0:
            keep_assembly.append(False)
        else:
            keep_assembly.append(True)

    is_member = np.array(is_member)

    return patterns[keep_assembly], is_member[keep_assembly], keep_assembly, is_member
