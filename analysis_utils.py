# analysis_utils.py

import numpy as np
import scipy.constants as const

class analysis_Helper:
    def __init__(self, G=150):
        self.G = G

    def Rs(self, Q):
        return self.G / Q

    def Xs(self, f, f0, X0):
        return -2 * self.G * (f - f0) / f0 + X0

    def sigmaRX(self, Rs, Xs, freq0):
        omega = 2 * np.pi * freq0
        mu0 = const.mu_0
        sigma = omega * mu0 * (2 * Rs * Xs / (Rs**2 + Xs**2)**2 +
                               (Xs**2 - Rs**2) / (Rs**2 + Xs**2)**2 * 1j)
        sigma1 = np.real(sigma)
        sigma2 = np.imag(sigma)
        return sigma1, sigma2, sigma

    def sigmaTrunin(self, Rs, Xs, Rn):
        sigma1 = 4 * Rn**2 * Rs * Xs / (Rs**2 + Xs**2)**2
        sigma2 = 2 * Rn**2 * (Xs**2 - Rs**2) / (Rs**2 + Xs**2)**2
        return sigma1, sigma2, sigma1 + 1j * sigma2

    def deltaLambda(self, freq, temp, G=192):
        f0 = freq[np.where(temp <= 5.0)[0][0]]
        dL = -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)
        return dL

    def deltaLFit(self, temp, Tc, lLondon, l, eps, l0):
        dl = lLondon * np.sqrt(1 + eps / l) / np.sqrt(1 - (temp / Tc)**4) - l0
        return dl
