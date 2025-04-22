# fg004_analysis.py
# FG004 Analysis Script — 11-3-2024

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from lmfit import Model, Parameters
import SCconductivity as sc
from analysis_utils import analysis_Helper  # <--- import helper

#? Change deltaLambda function in analysis_Helper
class MyAnalysisHelper(analysis_Helper):
    def deltaLambda(self, freq, temp, G=192):
        f0 = np.max(freq)
        return -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)

# Use your custom version
helper = MyAnalysisHelper(G=150)

# ---------------------- Constants & Setup ----------------------
G = 150
fileName = "data/FG004_throughTc.txt"
helper = MyAnalysisHelper(G)  # <--- create instance

# ---------------------- Load and Filter Data ----------------------
df = pd.read_csv(fileName, sep=r'\s+')
df.columns = ["Time", "Temp", "MKS1000", "LowerEdge1", "Bandwidth", "Freq_raw", "Q0", "LowerEdge2", "Loss", "Max_Freq"]

#! Correct the frequency using presusre information
df["Freq"] = df["Freq_raw"] - 750 * (df["MKS1000"] - 1000)

min_freq = df["Freq"].min()
df = df.query("Freq != @min_freq").reset_index(drop=True)

Tc_guess = 9.2
df_filtered = df.query("Temp >= @Tc_guess * 0.1 and Temp <= 9").reset_index(drop=True)

# --------------------------- Derived Quantities ---------------------------
RsData = helper.Rs(df_filtered["Q0"])
sigman = G / df_filtered["Q0"].iloc[0]
f0_ref = df_filtered["Freq"].iloc[0] + 11e3
XsData = helper.Xs(df_filtered["Freq"], f0_ref, sigman)
print(f"σₙ = {sigman:.4e}")

# ---------------------- Sigma Calculation ----------------------
sigma1, sigma2, sigma = helper.sigmaRX(RsData, XsData, df_filtered["Freq"].iloc[0])
sigma1T, sigma2T, sigmaT = helper.sigmaTrunin(RsData, XsData, sigman)

# ----------------------------- SC Conductivity -----------------------------
#! What are these numbers
freqS = df["Freq"].max()
Tc = 8.9
tempS = np.linspace(2, Tc - 1e-3, 1000)
Gamma = 0.06 * 0.05
sigman = 1 / (152e-9 * 1e-2)
G = 192

mySc = sc.SCconductivity(Tc, freqS, Gamma, tempS, sigman)
Q = mySc.Q()
deltaf = mySc.deltaf()
ZsS = mySc.Zs()
s1S, s2S = mySc.sigma()

# ---------------------------- deltaLambda ----------------------------
df_dl = df.query("Temp >= @Tc_guess * 0.1 and Temp <= @Tc - 0.3").reset_index(drop=True)
deltaL = helper.deltaLambda(df_dl["Freq"], df_dl["Temp"])

# ---------------------- Plotting ----------------------

plt.figure()
plt.plot(df["Temp"], df["Freq"], '.-', label="Corrected Freq")
plt.plot(df["Temp"], df["Freq_raw"], 'o', label="Raw Freq")
plt.xlabel('Temperature [K]')
plt.ylabel('Frequency [Hz]')
plt.xlim([4, 13])
plt.legend()

helper.plot_q0_vs_temp(df)

plt.figure()
plt.semilogy(df["Temp"], df["MKS1000"] - 1000, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [mbar]')
plt.xlim([4, 13])
plt.show()

####

plt.figure()
plt.plot(temp[idx1], deltaL * 1e10, '.-', label='Data')
plt.plot(temp[idx1], deltaLFit(temp[idx1], Tc - 0.2, 1000, 1, 1, 0), label='Initial Guess')
plt.xlabel("Temperature [K]")
plt.ylabel("$\Delta \lambda\ [\mathring{\mathrm{A}}]$")
plt.grid(True)

gmodel = Model(deltaLFit)
params = Parameters()
params.add('Tc', min=8.6, max=8.9, value=Tc)
params.add('lLondon', min=100, max=2000, value=500)
params.add('l', min=400, max=600, value=600)
params.add('eps', min=0.1, max=600, value=620, vary=False)
params.add('l0', min=100, max=600, value=610)

result = gmodel.fit(deltaL * 1e10, temp=temp[idx1], params=params)

print(result.fit_report())
plt.plot(temp[idx1], result.best_fit, label='Fit')
plt.legend()

# ---------------------------- Final Plots ----------------------------

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx], freq[idx] / 1e6, label='Data')
ax[0].set_xlim([7, 9])
ax1 = ax[0].twinx()
ax1.plot(tempS, deltaf / 1e6, '-.', color='orange', label='MB model')
ax1.set_xlim([7, 9])
ax[0].set_ylabel("Frequency [MHz]")
ax[1].semilogy(temp[idx], Q0[idx])
ax[1].semilogy(tempS, Q, '-')
ax[1].set_ylabel("Q$_0$")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].semilogy(temp[idx], RsData)
ax[0].semilogy(tempS, np.real(ZsS))
ax[0].set_ylabel("R$_s$ [$\Omega$]")
ax[1].plot(temp[idx], XsData / 1e-3)
ax1 = ax[1].twinx()
ax1.plot(tempS, np.imag(ZsS) / 1e-3, 'c')
ax[1].set_ylabel("X$_s$ [mΩ]")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx] / Tc, sigma1 / 2e8)
ax[0].plot(temp[idx] / Tc, sigma1T)
ax0 = ax[0].twinx()
ax0.plot(tempS / Tc, s1S)
ax[0].set_ylabel("σ₁")

ax[1].plot(temp[idx] / Tc, sigma2 / 2e8)
ax[1].plot(temp[idx] / Tc, sigma2T)
ax1 = ax[1].twinx()
ax1.plot(tempS / Tc, s2S)
ax[1].set_ylabel("σ₂")
ax[1].set_xlabel("T / Tc")

plt.tight_layout()
plt.show()
