import multiprocessing
import os
import pickle
import pandas as pd
from joblib import Parallel, delayed
from ripple_heterogeneity.assembly import assembly_reactivation


def get_assembly_patterns(basepath):
    # initialize session
    m1 = assembly_reactivation.AssemblyReact(basepath, weight_dt=0.02)
    # load data
    m1.load_data()
    # check if no cells were found
    if m1.cell_metrics.shape[0] == 0:
        return None
    # restrict to pre/task/post epochs
    m1.restrict_epochs_to_pre_task_post()
    # get weights for task outside ripples
    # % (TODO: use more robust method to locate epochs than index)
    m1.get_weights(m1.epochs[1][~m1.ripples])
    return m1

def session_loop(basepath, save_path):

    save_file = os.path.join(
        save_path, basepath.replace(os.sep, "_").replace(":", "_") + ".pkl"
    )
    if os.path.exists(save_file):
        return
    results = get_assembly_patterns(basepath)
    # save file
    with open(save_file, "wb") as f:
        pickle.dump(results, f)


def run(df, save_path, parallel=True):
    # find sessions to run
    basepaths = pd.unique(df.basepath)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if parallel:
        num_cores = multiprocessing.cpu_count()
        processed_list = Parallel(n_jobs=num_cores)(
            delayed(session_loop)(basepath, save_path) for basepath in basepaths
        )
    else:
        for basepath in basepaths:
            print(basepath)
            session_loop(basepath, save_path)
