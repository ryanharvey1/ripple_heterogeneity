import glob
import os
import pickle
import pandas as pd
from ripple_heterogeneity.utils import batch_analysis


def test_batchanalysis():
    def test_analysis(basepath, a="Real", b="Python", c="Is", d="Great", e="!"):
        results = {
            "basepath": os.path.exists(basepath),
            "exists?": basepath,
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "e": e,
        }
        return results

    df = pd.read_csv(r"Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv")
    df = df.iloc[0:10]
    
    save_path = r"Z:\home\ryanh\projects\ripple_heterogeneity\batch_analysis_test"

    batch_analysis.run(
        df,
        save_path,
        test_analysis,
        parallel=False,
        overwrite=True,
        a="fake",
        b="julksdjflm",
    )

    sessions = glob.glob(save_path + os.sep + "*.pkl")
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)

        assert results['a'] == 'fake'
        assert results['b'] == 'julksdjflm'
        assert results['c'] == 'Is'
        assert results['d'] == 'Great'
        assert results['e'] == '!'

    # test parallel
    batch_analysis.run(
        df,
        save_path,
        test_analysis,
        parallel=True,
        overwrite=True,
        a="fake",
        b="julksdjflm",
    )

    sessions = glob.glob(save_path + os.sep + "*.pkl")
    for session in sessions:
        with open(session, "rb") as f:
            results = pickle.load(f)
            
        assert results['a'] == 'fake'
        assert results['b'] == 'julksdjflm'
        assert results['c'] == 'Is'
        assert results['d'] == 'Great'
        assert results['e'] == '!'