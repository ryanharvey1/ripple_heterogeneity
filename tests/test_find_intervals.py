from ripple_heterogeneity.utils import functions

def test_find_interval():
    # Test empty list
    logical = []
    expected_output = []
    assert functions.find_interval(logical) == expected_output

    # Test list with length 1
    logical = [True]
    expected_output = [(0, 0)]
    assert functions.find_interval(logical) == expected_output

    logical = [False]
    expected_output = []
    assert functions.find_interval(logical) == expected_output

    # Test list with alternating True and False values
    logical = [True, False, True, False]
    expected_output = [(0, 0), (2, 2)]
    assert functions.find_interval(logical) == expected_output

    # Test list with multiple intervals of True values
    logical = [0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1]
    expected_output = [(2, 4), (6, 7), (10, 11)]
    assert functions.find_interval(logical) == expected_output

    logical = [1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1]
    expected_output = [(0, 2), (4, 5), (9, 10)]
    assert functions.find_interval(logical) == expected_output

    # Test list with all True values
    logical = [1, 1, 1, 1, 1]
    expected_output = [(0, 4)]
    assert functions.find_interval(logical) == expected_output