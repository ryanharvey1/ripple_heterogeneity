import multiprocessing
import os
import pickle
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm
import traceback


def main_loop(basepath, save_path, func, overwrite, skip_if_error, **kwargs):
    """
    main_loop: file management & run function
    Inputs:
        basepath: str path to session
        save_path: str path to save results to (will be created if it doesn't exist)
        func: function to run on each basepath in df (see run)
        overwrite: bool whether to overwrite existing files in save_path
        skip_if_error: bool whether to skip if an error occurs
        kwargs: dict of keyword arguments to pass to func (see run)
    """
    # get file name from basepath
    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    # if file exists and overwrite is False, skip
    if os.path.exists(save_file) & ~overwrite:
        return

    # calc some features
    if skip_if_error:
        try:
            results = func(basepath, **kwargs)
        except Exception:
            traceback.print_exc()
            print(f"Error in {basepath}")
            return
    else:
        results = func(basepath, **kwargs)

    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def run(
    df,
    save_path,
    func,
    parallel=True,
    verbose=False,
    overwrite=False,
    skip_if_error=False,
    num_cores=None,
    **kwargs,
):
    """
    Inputs:
        df: pandas dataframe with basepath column
        save_path: str path to save results to (will be created if it doesn't exist)
        func: function to run on each basepath in df (see main_loop)
        parallel: bool whether to run in parallel or not
        verbose: bool whether to print progress
        overwrite: bool whether to overwrite existing files in save_path
        skip_if_error: bool whether to skip if an error occurs
        kwargs: dict of keyword arguments to pass to func
    """
    # find sessions to run
    basepaths = pd.unique(df.basepath)
    # create save_path if it doesn't exist
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    # run in parallel if parallel is True
    if parallel:
        # get number of cores
        if num_cores is None:
            num_cores = multiprocessing.cpu_count()
        # run in parallel
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(main_loop)(
                basepath, save_path, func, overwrite, skip_if_error, **kwargs
            )
            for basepath in tqdm(basepaths)
        )
    else:
        # run in serial
        for basepath in tqdm(basepaths):
            if verbose:
                print(basepath)
            # run main_loop on each basepath in df
            main_loop(basepath, save_path, func, overwrite, skip_if_error, **kwargs)
