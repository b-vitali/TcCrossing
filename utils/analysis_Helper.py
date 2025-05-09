# analysis_utils.py

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import os

class analysis_Helper:
    def __init__(self, G=150, save=False, filename='data', save_dir="figs"):
        self.G = G
        self.save = save
        self.filename = filename
        self.save_dir = save_dir+"/"+filename

        if self.save and not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

    def _handle_plot(self, plotname):
        fig = plt.gcf()  # Get the current figure
        fig.text(0.98, 0.02, f"{self.filename}", ha='right', va='center', fontsize=10, color='black')

        if self.save:
            full_path = os.path.join(self.save_dir, plotname)
            plt.savefig(full_path)
            print(f"Saved plot to {full_path}")
            plt.close()
        else:
            plt.show()

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

    def deltaLambda(self, freq, temp):
        f0 = np.max(freq)
        return -self.G * (freq - f0) / (const.pi * const.mu_0 * f0**2)

    def deltaLFit(self, temp, Tc, lLondon, l, eps, l0):
        dl = lLondon * np.sqrt(1 + eps / l) / np.sqrt(1 - (temp / Tc)**4) - l0
        return dl

    def deltaLFit_simple(self, temp, Tc, A, l0):
        dl = A / np.sqrt(1 - (temp / Tc)**4) - l0
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
        self._handle_plot("plot_frequency_vs_temp.png")


    def plot_pressure_vs_temp(self, df):
        plt.figure()
        plt.semilogy(df["Temp"], df["MKS1000"] - 1000, '.', label="Pressure")
        plt.xlabel('Temperature [K]')
        plt.ylabel('Pressure [mbar]')
        plt.xlim([4, 13])
        plt.grid(True, which='both')
        plt.title('Pressure vs Temperature')
        plt.tight_layout()
        self._handle_plot("plot_pressure_vs_temp.png")

    def plot_q0_vs_temp(self, df):
        plt.figure()
        plt.semilogy(df["Temp"], df["Q0"], '.')
        plt.xlabel('Temperature [K]')
        plt.ylabel('$Q_0$')
        plt.xlim([4, 11])
        plt.grid(True)
        plt.title('$Q_0$ vs Temperature')

        plt.tight_layout()
        self._handle_plot("plot_q0_vs_temp.png")

    def plot_dlambda_fit(self, temp, deltaL, fit_result):
        import matplotlib.pyplot as plt

        plt.figure()
        
        # Extract fit parameters, uncertainties, and chi-squared value
        params = fit_result.params
        chi_squared = fit_result.chisqr

        # Create the label for the fit line including parameters and chi-squared
        param_lines = [f"{param}={params[param].value:.2e}±{params[param].stderr:.1e}" for param in params]
        chi_squared_line = f"$\chi^2$={chi_squared:.2f}"
        full_fit_label = "Fit:\n" + "\n".join(param_lines + [chi_squared_line])

        # Plot data and fit with updated legend label
        plt.plot(temp, deltaL * 1e10, label="Data")
        plt.plot(temp, fit_result.best_fit, label=full_fit_label)

        plt.xlabel("Temperature [K]")
        plt.ylabel("$\Delta \lambda\ [\mathring{\mathrm{A}}]$")
        plt.grid(True)
        plt.legend(loc="best", fontsize=8)
        plt.title("Delta Lambda Fit")
        plt.tight_layout()
        self._handle_plot("plot_dlambda_fit.png")

    def plot_freq_q0_dual(self, df, tempS, deltaf, Q):
        fig, ax = plt.subplots(2, 1)

        # === Top Plot: Frequency & Delta f ===
        # Frequency
        line1, = ax[0].plot(df["Temp"], df["Freq"] / 1e1, label="Resonant Frequency", color='blue')
        ax[0].set_xlim([7, 9])
        ax[0].set_ylabel("Frequency [Hz]")
        ax[0].grid(True)

        # Delta f on twin axis
        ax1 = ax[0].twinx()
        line2, = ax1.plot(tempS, deltaf / 1e1, '-.', label="$\Delta f$", color='orange')
        ax1.set_xlim([7, 9])
        ax1.set_ylabel("$\Delta f$ [Hz]")

        # Combine legends from both axes
        lines = [line1, line2]
        labels = [line.get_label() for line in lines]
        ax[0].legend(lines, labels, loc='upper left')

        # === Bottom Plot: Q0 ===
        ax[1].semilogy(df["Temp"], df["Q0"], 'o', label="$Q_0$ Data", color='green', markersize=4)
        ax[1].semilogy(tempS, Q, '-', label="Q₀ Model", color='purple')
        ax[1].set_ylabel("$Q_0$")
        ax[1].set_xlabel("Temperature [K]")
        ax[1].legend(loc='best')
        ax[1].grid(True)

        # === Figure Title ===
        fig.suptitle("Frequency and $Q_0$ vs Temperature")
        plt.tight_layout()
        self._handle_plot("plot_freq_q0_dual.png")

    def plot_rs_xs_dual(self, df, tempS, RsData, ZsS, XsData):
        fig, ax = plt.subplots(2, 1)

        # === Top Plot: Rs ===
        line1, = ax[0].semilogy(df["Temp"], RsData, label="$R_s$ Data", color='tab:red')
        line2, = ax[0].semilogy(tempS, np.real(ZsS), label="Re($Z_s$)", color='tab:blue')
        ax[0].set_ylabel("$R_s$ [$\Omega$]")
        ax[0].grid(True)

        # Combine legends
        lines = [line1, line2]
        labels = [line.get_label() for line in lines]
        ax[0].legend(lines, labels, loc='best')

        # === Bottom Plot: Xs ===
        line3, = ax[1].plot(df["Temp"], XsData / 1e-3, label="$X_s$ Data", color='tab:green')
        ax1 = ax[1].twinx()
        line4, = ax1.plot(tempS, np.imag(ZsS) / 1e-3, '-.', label="Im($Z_s$)", color='tab:orange')
        ax[1].set_ylabel("$X_s$ [m$\Omega$]")
        ax[1].set_xlabel("Temperature [K]")
        ax[1].grid(True)

        # Combine legends
        lines = [line3, line4]
        labels = [line.get_label() for line in lines]
        ax[1].legend(lines, labels, loc='best')

        # === Title ===
        fig.suptitle("$R_s$ and $X_s$ vs Temperature")
        plt.tight_layout()
        self._handle_plot("plot_rs_xs_dual.png")


    def plot_sigma_dual(self, df, tempS, Tc, sigma1, sigma2, sigma1T, sigma2T, s1S, s2S):
        fig, ax = plt.subplots(2, 1)

        # === Top Plot: Sigma1 ===
        line1, = ax[0].plot(df["Temp"] / Tc, sigma1 / 2e8, label="RX", color='tab:blue')
        line2, = ax[0].plot(df["Temp"] / Tc, sigma1T, label="Trunin", color='tab:green')
        ax0 = ax[0].twinx()
        line3, = ax0.plot(tempS / Tc, s1S, '--', label="Mattis-Bardeen", color='tab:orange')

        ax[0].set_ylabel(r"$\sigma_1$")
        ax[0].grid(True)

        lines = [line1, line2, line3]
        labels = [line.get_label() for line in lines]
        ax[0].legend(lines, labels, loc='best')

        # === Bottom Plot: Sigma2 ===
        line4, = ax[1].plot(df["Temp"] / Tc, sigma2 / 2e8, label="RX", color='tab:blue')
        line5, = ax[1].plot(df["Temp"] / Tc, sigma2T, label="Trunin", color='tab:green')
        ax1 = ax[1].twinx()
        line6, = ax1.plot(tempS / Tc, s2S, '--', label="Mattis-Bardeen", color='tab:orange')

        ax[1].set_ylabel(r"$\sigma_2$")
        ax[1].set_xlabel("Reduced Temperature $T/T_c$")
        ax[1].grid(True)

        lines = [line4, line5, line6]
        labels = [line.get_label() for line in lines]
        ax[1].legend(lines, labels, loc='best')

        # === Title ===
        fig.suptitle(r"$\sigma_1$ and $\sigma_2$ vs Reduced Temperature")
        plt.tight_layout()
        self._handle_plot("plot_sigma_dual.png")
