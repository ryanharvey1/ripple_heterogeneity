import glob
import multiprocessing
import os
import pickle
from joblib import Parallel, delayed
import pandas as pd
from tqdm import tqdm
import traceback


def main_loop(
    basepath: str,
    save_path: str,
    func,
    overwrite: bool,
    skip_if_error: bool,
    **kwargs,
) -> None:
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
    basepath = os.path.normpath(basepath)
    save_path = os.path.normpath(save_path)
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
    df: pd.core.frame.DataFrame,
    save_path: str,
    func,
    parallel: bool = True,
    verbose: bool = False,
    overwrite: bool = False,
    skip_if_error: bool = False,
    num_cores: int = None,
    **kwargs,
) -> None:
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


def load_results(
    save_path: str, verbose: bool = False, add_save_file_name: bool = False
) -> pd.core.frame.DataFrame:
    """
    load_results: load results (pandas dataframe) from a pickle file

    This is the most basic results loader and
        **will only work if your output was a pandas dataframe (long format)**

    This will have to be adapted if your output was more complicated, but you can
        use this function as an example.

    Inputs:
        save_path: str path to save results
        verbose: bool whether to print progress
        add_save_file_name: bool whether to add a column with the save file name

    Returns:
        results: pandas dataframe with results
    """

    if not os.path.exists(save_path):
        raise ValueError(f"folder {save_path} does not exist")

    sessions = glob.glob(save_path + os.sep + "*.pkl")

    results = pd.DataFrame()

    for session in sessions:
        if verbose:
            print(session)
        with open(session, "rb") as f:
            results_ = pickle.load(f)
        if results_ is None:
            continue

        if add_save_file_name:
            results_["save_file_name"] = os.path.basename(session)

        results = pd.concat([results, results_], ignore_index=True)

    return results
