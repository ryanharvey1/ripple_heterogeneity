import multiprocessing
import os
import pickle
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm


def main_loop(basepath, save_path, func, overwrite, **kwargs):
    """
    main_loop: file management & run function
    Inputs:
        basepath: str
        save_path: str
        func: function
        overwrite: bool
        kwargs: dict
    """
    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file) & ~overwrite:
        return

    # calc some features
    results = func(basepath, **kwargs)

    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def run(df, save_path, func, parallel=True, verbose=False, overwrite=False, **kwargs):
    """
    Inputs:
        df: pandas dataframe
        save_path: str
        func: function
        parallel: bool
        verbose: bool
        overwrite: bool
        kwargs: dict
    """
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(main_loop)(basepath, save_path, func, overwrite, **kwargs)
            for basepath in tqdm(basepaths)
        )
    else:
        for basepath in tqdm(basepaths):
            if verbose:
                print(basepath)
            main_loop(basepath, save_path, func, overwrite, **kwargs)
