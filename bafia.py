import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as const
from lmfit import Model, Parameters
import SCconductivity as sc
from analysis_utils import analysis_Helper  # <--- import helper

# ---------------------- Constants & Setup ----------------------
#! Here G=150, below is set to 192 ?!
G = 150
fileName = "data/20250128_FNAL_103.txt"
helper = analysis_Helper(G)  # <--- create instance

# ---------------------- Load and Filter Data ----------------------
df = pd.read_csv(fileName, sep=r'\s+')
df.columns = ["Time", "Temp", "MKS1000", "LowerEdge1", "Bandwidth", "Freq", "Q0", "LowerEdge2", "Loss", "Max_Freq"]

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

# ---------------------- SC Conductivity ----------------------
#! What are these numbers
freqS = 650e6
Tc = 8.7
tempS = np.linspace(2, Tc - 1e-3, 1000)
Gamma = 0.06
sigman_MB = 1 / (152e-9 * 1e-2)
G = 192

mySc = sc.SCconductivity(Tc, freqS, Gamma, tempS, sigman_MB)
Q = mySc.Q()
deltaf = mySc.deltaf()
ZsS = mySc.Zs()
s1S, s2S = mySc.sigma()

# ---------------------- deltaLambda ----------------------
df_dl = df.query("Temp >= @Tc_guess * 0.1 and Temp <= @Tc - 0.3").reset_index(drop=True)
deltaL = helper.deltaLambda(df_dl["Freq"], df_dl["Temp"])

# ---------------------- Fitting ----------------------
gmodel = Model(helper.deltaLFit)
params = Parameters()
params.add('Tc', min=8.6, max=8.9, value=Tc)
params.add('lLondon', min=100, max=2000, value=500)
params.add('l', min=400, max=600, value=600)
params.add('eps', min=0.1, max=600, value=620, vary=False)
params.add('l0', min=100, max=600, value=610)

result = gmodel.fit(deltaL * 1e10, temp=df_dl["Temp"], params=params)

# ---------------------- Plotting ----------------------
helper.plot_frequency_vs_temp(df)
helper.plot_q0_vs_temp(df)
helper.plot_dlambda_fit(df_dl["Temp"], deltaL, result)

print(result.fit_report())
plt.show()

fig, ax = plt.subplots(2, 1)
ax[0].plot(df_filtered["Temp"], df_filtered["Freq"] / 1e6)
ax[0].set_xlim([7, 9])
ax1 = ax[0].twinx()
ax1.plot(tempS, deltaf / 1e6, '-.')
ax1.set_xlim([7, 9])
ax[0].set_ylabel("Frequency [MHz]")

ax[1].semilogy(df_filtered["Temp"], df_filtered["Q0"])
ax[1].semilogy(tempS, Q, '-')
ax[1].set_ylabel("Q$_0$")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].semilogy(df_filtered["Temp"], RsData)
ax[0].semilogy(tempS, np.real(ZsS))
ax[0].set_ylabel("R$_s$ [$\Omega$]")

ax[1].plot(df_filtered["Temp"], XsData / 1e-3)
ax1 = ax[1].twinx()
ax1.plot(tempS, np.imag(ZsS) / 1e-3, 'c')
ax[1].set_ylabel("X$_s$ [m$\Omega$]")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].plot(df_filtered["Temp"] / Tc, sigma1 / 2e8, label="RX")
ax[0].plot(df_filtered["Temp"] / Tc, sigma1T, label="Trunin")
ax0 = ax[0].twinx()
ax0.plot(tempS / Tc, s1S, linestyle="--", label="Mattis-Bardeen")
ax[0].set_ylabel("$\sigma_1$")
ax[0].legend()

ax[1].plot(df_filtered["Temp"] / Tc, sigma2 / 2e8)
ax[1].plot(df_filtered["Temp"] / Tc, sigma2T)
ax1 = ax[1].twinx()
ax1.plot(tempS / Tc, s2S)
ax[1].set_ylabel("$\sigma_2$")
ax[1].set_xlabel("Temperature [K]")

plt.tight_layout()
plt.show()
