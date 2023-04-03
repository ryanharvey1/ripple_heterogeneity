import numpy as np
from scipy.interpolate import interp1d
from scipy.signal.windows import dpss
import pandas as pd
from numba import jit


def getfgrid(Fs, nfft, fpass):
    """
    get frequency grid for evaluation
    Inputs:
        Fs: sampling frequency
        nfft: number of points for fft
        fpass: frequency range to evaluate
    Outputs:
        f: frequency vector
        findx: frequency index
    """
    df = Fs / nfft
    f = np.arange(0, Fs + df, df)
    f = f[0:nfft]
    findx = np.logical_and(f >= fpass[0], f <= fpass[-1])
    f = f[findx]
    return f, findx


def dpsschk(tapers, N, Fs):
    """
    check tapers
    Inputs:
        tapers: [NW, K] or [tapers, eigenvalues]
        N: number of points for fft
        Fs: sampling frequency
    Outputs:
        tapers: [tapers, eigenvalues]
    """
    tapers, eigs = dpss(N, NW=tapers[0], Kmax=tapers[1], sym=False, return_ratios=True)
    tapers = tapers * np.sqrt(Fs)
    tapers = tapers.T
    return tapers


def mtfftpt(data, tapers, nfft, t, f, findx):
    """
    mt fft for point process times
    Inputs:
        data: 1d array of spike times (in seconds)
        tapers: [NW, K] or [tapers, eigenvalues]
        nfft: number of points for fft
        t: time vector
        f: frequency vector
        findx: frequency index
    Outputs:
        J: fft of data
        Msp: number of spikes in data
        Nsp: number of spikes in data
    """
    C = 1
    K = tapers.shape[1]
    nfreq = len(f)
    H = np.zeros((nfft, tapers.shape[1]), dtype=np.complex128)
    for i in range(tapers.shape[1]):
        H[:, i] = np.fft.fft(tapers[:, i], nfft, axis=0)

    H = H[findx, :]
    w = 2 * np.pi * f
    dtmp = data
    indx = np.logical_and(dtmp >= np.min(t), dtmp <= np.max(t))
    if len(indx):
        dtmp = dtmp[indx]
    Nsp = len(dtmp)
    Msp = Nsp / len(t)

    if Msp != 0:
        data_proj = []
        for i in range(tapers.shape[1]):
            ff = interp1d(t, tapers[:, i])
            data_proj.append(ff(dtmp))
        data_proj = np.array(data_proj).T
        exponential = np.exp(np.atleast_2d(-1j * w).T * (dtmp - t[0]))
        J = np.dot(exponential, data_proj) - H * Msp
    else:
        J = np.zeros((nfreq, tapers.shape[1]))

    return J, Msp, Nsp

# @jit(nopython=True)
def mtspectrumpt(data, Fs, fpass, tapers):
    """
    mtspectrumpt from chronux toolbox
    Inputs:
        data: 1d array of spike times (in seconds)
        Fs: sampling frequency
        fpass: frequency range to evaluate
        tapers: [NW, K] or [tapers, eigenvalues]
    Outputs:
        S: power spectrum
    """
    mintime = np.min(data)
    maxtime = np.max(data)
    dt = 1 / Fs
    t = np.arange(mintime - dt, maxtime + dt, dt)
    N = len(t)
    nfft = np.max(
        [int(2 ** np.ceil(np.log2(N))), N]
    )  # number of points in fft of prolates
    f, findx = getfgrid(Fs, nfft, fpass)
    tapers = dpsschk(tapers, N, Fs)
    J, Msp, Nsp = mtfftpt(data, tapers, nfft, t, f, findx)
    S = np.real(np.mean(np.conj(J) * J, 1))
    return pd.Series(index=f, data=S)


def mtfftc(data, tapers, nfft, Fs):
    """ 
    mt fft for continuous data
    Inputs:
        data: 1d array 
        tapers: [NW, K] or [tapers, eigenvalues]
        nfft: number of points for fft
        Fs: sampling frequency
    Outputs:
        J: fft of data
    """
    NC = len(data)
    NK, K = tapers.shape
    assert NK == NC, "length of tapers is incompatible with length of data"
    tmp = np.repeat(np.atleast_2d(data), K, 0).T
    tmp2 = tmp * tapers
    J = np.fft.fft(tmp2.T, nfft) / float(Fs)
    return J


def mtspectrumc(data, Fs, fpass, tapers):
    """ 
    mtspectrumc from chronux toolbox
    Inputs:
        data: 1d array
        Fs: sampling frequency
        fpass: frequency range to evaluate
        tapers: [NW, K] or [tapers, eigenvalues]
    Outputs:
        S: power spectrum
    """
    N = len(data)
    nfft = np.max(
        [int(2 ** np.ceil(np.log2(N))), N]
    )  # number of points in fft of prolates
    f, findx = getfgrid(Fs, nfft, fpass)
    tapers = dpsschk(tapers, N, Fs)
    J = mtfftc(data, tapers, nfft, Fs)
    J = J.T[findx, :]
    S = np.real(np.mean(np.conj(J) * J, 1))
    return pd.Series(index=f, data=S)
