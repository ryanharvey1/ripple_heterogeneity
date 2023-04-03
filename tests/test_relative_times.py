import numpy as np
from ripple_heterogeneity.utils import functions

def test_relative_times():
    # Test 1: basic test case
    t = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    intervals = np.array([[1, 3], [4, 6], [7, 9]])
    assert np.allclose(
        functions.relative_times(t, intervals),
        (
            np.array([np.nan, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]),
            np.array([np.nan, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]),
            np.array([np.nan, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]),
        ),
        atol=1e-6,
        equal_nan=True,
    )

    # Test 2: with values assigned to each interval
    t = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    intervals = np.array([[1, 3], [4, 6], [7, 9]])
    values = np.array([10, 20, 30])
    assert np.allclose(
        functions.relative_times(t, intervals, values),
        (
            np.array([np.nan, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0]),
            np.array([np.nan, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0]),
            np.array([np.nan, 10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 30.0, 30.0, 30.0]),
        ),
        atol=1e-6,
        equal_nan=True,
    )

    # Test 3: when t is outside of all intervals
    t = np.array([-2, -1])
    intervals = np.array([[1, 3], [4, 6], [7, 9]])
    assert np.allclose(
        functions.relative_times(t, intervals),
        (
            np.array([np.nan, np.nan]),
            np.array([np.nan, np.nan]),
            np.array([np.nan, np.nan]),
        ),
        atol=1e-6,
        equal_nan=True,
    )

    # Test 4: when t is at the start of the first interval
    t = np.array([1])
    intervals = np.array([[1, 3], [4, 6], [7, 9]])
    assert np.allclose(
        functions.relative_times(t, intervals),
        (np.array([0.0]), np.array([0.0]), np.array([0.0])),
        atol=1e-6,
    )

    # Test 5: when t is at the start or end of an interval
    t = np.array([1, 4, 7])
    intervals = np.array([[1, 3], [4, 6], [7, 9]])
    expected_rt = np.array([0.0, 0.0, 0.0])
    expected_intervalID = np.array([0.0, 1.0, 2.0])
    expected_rt_values = np.array([10.0, 20.0, 30.0])
    assert np.array_equal(
        functions.relative_times(t, intervals, np.array([10, 20, 30]))[0], expected_rt
    )
    assert np.array_equal(
        functions.relative_times(t, intervals, np.array([10, 20, 30]))[1],
        expected_intervalID,
    )
    assert np.array_equal(
        functions.relative_times(t, intervals, np.array([10, 20, 30]))[2],
        expected_rt_values,
    )
