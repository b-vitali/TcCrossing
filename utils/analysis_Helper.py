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

        # If there is a 'raw' plot both corrected and raw frequencies
        if "Freq_raw" in df.columns:
            plt.plot(df["Temp"], df["Freq"] , '.', label="Corrected Freq")
            plt.plot(df["Temp"], df["Freq_raw"] , '.', label="Raw Freq", alpha=0.6)
            plt.legend()
        else:
            plt.plot(df["Temp"], df["Freq"] , '.', label="Frequency")

        plt.xlabel('Temperature [K]')
        plt.ylabel('Frequency [Hz]')
        plt.grid(True)
        plt.title('Frequency vs Temperature')
        
        plt.tight_layout()

    def plot_pressure_vs_temp(self, df):
        plt.figure()
        plt.semilogy(df["Temp"], df["MKS1000"] - 1000, '.', label="Pressure")
        plt.xlabel('Temperature [K]')
        plt.ylabel('Pressure [mbar]')
        plt.xlim([4, 13])
        plt.grid(True, which='both')
        plt.title('Pressure vs Temperature')
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

    def plot_freq_q0_dual(self, df, tempS, deltaf, Q):
        fig, ax = plt.subplots(2, 1, figsize=(6, 6))

        # Frequency vs Temp + delta f (twin axis)
        ax[0].plot(df["Temp"], df["Freq"] / 1e6, label="Resonant Freq")
        ax[0].set_xlim([7, 9])
        ax[0].set_ylabel("Frequency [MHz]")

        ax1 = ax[0].twinx()
        ax1.plot(tempS, deltaf / 1e6, '-.', color='orange', label="Δf")
        ax1.set_xlim([7, 9])
        ax1.set_ylabel("Δf [MHz]")

        # Q0 vs Temp
        ax[1].semilogy(df["Temp"], df["Q0"], label="Q₀ Data")
        ax[1].semilogy(tempS, Q, '-', label="Q₀ Model")
        ax[1].set_ylabel("Q$_0$")
        ax[1].set_xlabel("Temperature [K]")

        ax[0].grid(True)
        ax[1].grid(True)

        plt.tight_layout()

    def plot_rs_xs_dual(self, df, tempS, RsData, ZsS, XsData):
        fig, ax = plt.subplots(2, 1, figsize=(6, 6))

        # Rs plot
        ax[0].semilogy(df["Temp"], RsData, label="Rs Data")
        ax[0].semilogy(tempS, np.real(ZsS), label="Re(Zs)")
        ax[0].set_ylabel("R$_s$ [$\Omega$]")
        ax[0].legend()
        ax[0].grid(True)

        # Xs plot
        ax[1].plot(df["Temp"], XsData / 1e-3, label="Xs Data")
        ax1 = ax[1].twinx()
        ax1.plot(tempS, np.imag(ZsS) / 1e-3, 'c', label="Im(Zs)")
        ax[1].set_ylabel("X$_s$ [m$\Omega$]")
        ax[1].set_xlabel("Temperature [K]")
        ax[1].legend(loc="upper left")
        ax1.legend(loc="upper right")
        ax[1].grid(True)

        plt.tight_layout()

    def plot_sigma_dual(self, df, tempS, Tc, sigma1, sigma2, sigma1T, sigma2T, s1S, s2S):
        fig, ax = plt.subplots(2, 1, figsize=(6, 6))

        # Sigma1 plot
        ax[0].plot(df["Temp"] / Tc, sigma1 / 2e8, label="RX")
        ax[0].plot(df["Temp"] / Tc, sigma1T, label="Trunin")
        ax0 = ax[0].twinx()
        ax0.plot(tempS / Tc, s1S, linestyle="--", label="Mattis-Bardeen")
        ax[0].set_ylabel(r"$\sigma_1$")
        ax[0].legend(loc="upper left")
        ax0.legend(loc="upper right")
        ax[0].grid(True)

        # Sigma2 plot
        ax[1].plot(df["Temp"] / Tc, sigma2 / 2e8)
        ax[1].plot(df["Temp"] / Tc, sigma2T)
        ax1 = ax[1].twinx()
        ax1.plot(tempS / Tc, s2S)
        ax[1].set_ylabel(r"$\sigma_2$")
        ax[1].set_xlabel("Temperature [K]")
        ax[1].grid(True)

        plt.tight_layout()
