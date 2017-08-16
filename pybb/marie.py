import numpy as np
import pywt
import matplotlib.pyplot as plt

def farge_denoise(sig, wtname, mode):
    
    # @sig : 1d signal
    # @wtname : wavelet name
    # @mode: signal extension mode for reconstruction
    
    # bulid wavelet, calculate maxlevel coeffs
    w = pywt.Wavelet(wtname)
    lev = pywt.dwt_max_level(len(sig), w)
    coeffs = pywt.wavedec(sig, w, level=lev)
    
    # total number of wavelet coefficients
    N = sum([len(i) for i in coeffs[1:]])
    
    # initial sigma square, eps, Sc, Si
    sigma2 = sum([sum([i*i for i in j]) for j in coeffs[1:]]) / N
    eps = np.sqrt(2.*np.log(N)*sigma2)
    Sc = pywt.wavedec(sig, w, level=lev)
    Si = pywt.wavedec(sig, w, level=lev)
    for i in range(len(coeffs[0])):
        Si[0][i] = 0.
    
    # main loop
    Ni_old = -1
    Ni = 0
    while Ni_old != Ni:
        Ni_old = Ni
        Nc = 0
        Ni = 0
        for i in range(1,len(coeffs)):
            for j in range(len(coeffs[i])):
                isNoise = abs(coeffs[i][j]) <= eps
                Ni += int(isNoise)
                Nc += int(not isNoise)
                Si[i][j] = float(isNoise) * coeffs[i][j]
                Sc[i][j] = float(not isNoise) * coeffs[i][j]
        sigma2 = sigma2 = sum([sum([i*i for i in j]) for j in Si]) / N
        eps = np.sqrt(2.*np.log(N)*sigma2)
        #print Ni, Nc, eps
    
    rec = pywt.waverec(Sc, w, mode)
    return rec
