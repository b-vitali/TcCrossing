# analysis_utils.py

import numpy as np
import matplotlib.pyplot as plt
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

    #! Why a different G?
    def deltaLambda(self, freq, temp, G=192):
        f0 = freq[np.where(temp <= 5.0)[0][0]]
        dL = -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)
        return dL

    def deltaLFit(self, temp, Tc, lLondon, l, eps, l0):
        dl = lLondon * np.sqrt(1 + eps / l) / np.sqrt(1 - (temp / Tc)**4) - l0
        return dl

    # ---------------------- Plotting Methods ----------------------

    def plot_frequency_vs_temp(self, df):
        plt.figure()
        plt.plot(df["Temp"], df["Freq"], '.')
        plt.xlabel('Temperature [K]')
        plt.ylabel('Frequency [Hz]')
        plt.xlim([4, 11])
        plt.ylim([6.4982e8, 6.4984e8])
        plt.grid(True)
        plt.title('Frequency vs Temperature')
        plt.tight_layout()

    def plot_q0_vs_temp(self, df):
        plt.figure()
        plt.plot(df["Temp"], df["Q0"], '.')
        plt.xlabel('Temperature [K]')
        plt.ylabel('$Q_0$')
        plt.xlim([4, 11])
        plt.grid(True)
        plt.title('$Q_0$ vs Temperature')
        plt.tight_layout()

    def plot_dlambda_fit(self, temp, deltaL, fit_result):

        #! Why do i need this here?
        Tc = 9.2

        plt.figure()
        plt.plot(temp, deltaL * 1e10, label="Data")
        plt.plot(temp, fit_result.best_fit, label="Fit")
        plt.plot(temp, self.deltaLFit(temp, Tc-0.2, 1000, 1, 1, 0 ))
        plt.xlabel("Temperature [K]")
        plt.ylabel("$\Delta \lambda\ [\mathring{\mathrm{A}}]$")
        plt.grid(True)
        plt.legend()
        plt.title("Delta Lambda Fit")
        plt.tight_layout()
