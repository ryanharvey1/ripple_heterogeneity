import numpy as np
from ripple_heterogeneity.assembly import similarity_index
import pytest

def test_similarityindex():
    patterns = [np.random.normal(0, .1, 100) for _ in range(10)]

    si, combos, pvalues = similarity_index.similarity_index(patterns)
    assert all(i <= 1 for i in si)
    assert combos[0,0] == 0
    assert all(i <= 1 for i in pvalues)