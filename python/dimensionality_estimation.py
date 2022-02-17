import sys
sys.path.append(r'D:\github\ripple_heterogeneity\python')
import functions,loading

import pandas as pd
import numpy as np 
import glob
import os
import pickle
import seaborn as sns
import nelpy as nel
from sklearn.decomposition import PCA
# from numba import jit


# @jit(nopython=True)
def SVCA_(X):
    """
    computes a cross-validated form of PCA: Shared Variance Component
    Analysis. Components are extracted from the covariance between neuron
    sets "ntrain" and "ntest" on training timepts "itrain". The variance of
    these covariance components is computed on testing timepts "itest". 
    This variance is the amount of reliable variance of that component 
    (because it's consistent across timepts).
    Inputs:
        X: (neurons x timepts)
    Outputs:
        scov (shared variance of each covariance component)
        varcov (total variance of each covariance component)

    From: https://github.com/MouseLand/stringer-et-al-2019/blob/master/utils.py    
    """
    
    NN,NT = X.shape

    # split cells into test and train
    norder = np.random.permutation(NN)
    nhalf = int(norder.size/2)
    ntrain = norder[:nhalf]
    ntest = norder[nhalf:]

    # split time into test and train
    torder = np.random.permutation(NT)
    thalf = int(torder.size/2)
    ttrain = torder[:thalf]
    ttest = torder[thalf:]

    cov = X[np.ix_(ntrain, ttrain)] @ X[np.ix_(ntest, ttrain)].T
    u = PCA(n_components=min(1024, nhalf-1), svd_solver='randomized').fit_transform(cov)
    u /= (u**2).sum(axis=0)**0.5
    v = cov.T @ u
    v /= (v**2).sum(axis=0)**0.5

    strain = u.T @ X[np.ix_(ntrain,ttest)]
    stest = v.T @ X[np.ix_(ntest,ttest)]

    # covariance k is uk.T * F * G.T * vk / npts
    scov = (strain * stest).mean(axis=1)
    varcov = (strain**2 + stest**2).mean(axis=1) / 2

    return scov,varcov,strain,stest


def SVCA(X,folds=10,verbose=False):
    """ calls SVCA_ """

    #  do multiple subsets and avg
    scov = []
    varcov = []
    strain = []
    stest = []
    svc_neur = []
    for i in range(folds):
        if verbose:
            print(str(i)+' of ',str(folds-1),'...')
        scov_,varcov_,strain_,stest_ = SVCA_(X)
        scov.append(scov_)
        varcov.append(varcov_)
        strain.append(strain_)
        stest.append(stest_)
        svc_neur.append(scov_ / varcov_)

    return svc_neur,scov,varcov,strain,stest