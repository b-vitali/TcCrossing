# fg004_analysis.py

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from lmfit import Model, Parameters
import SCconductivity as sc
from analysis_utils import analysis_Helper

# === Parameters ===
fileName = "FG004_throughTc.txt"
Tc = 8.9
Gamma = 0.06 * 0.05
sigman = 1 / (152e-9 * 1e-2)
G = 192

#! CHECK Slight change to the deltaLambda function
class analysis_Helper(analysis_Helper):
    def deltaLambda(self, freq, temp, G=192):
        f0 = np.max(freq)
        dL = -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)
        return dL

# === Instantiate Helper ===
bafia = analysis_Helper(G=G)

# === Load Data ===
data = np.loadtxt(fileName, skiprows=4)
temp = data[:, 1]
MKS1000 = data[:, 2]
freq1 = data[:, 5]
Q0 = data[:, 6]

# Pressure correction
freq = freq1 - 750 * (MKS1000 - 1000)

# Clean invalid data
valid_idx = np.where(freq != freq.min())[0]
temp = temp[valid_idx]
freq = freq[valid_idx]
freq1 = freq1[valid_idx]
MKS1000 = MKS1000[valid_idx]
Q0 = Q0[valid_idx]

# Filter for superconducting range
idx = np.where((temp <= 9) & (temp >= 0.1 * 9.2))
RsData = bafia.Rs(Q0[idx])
sigman_fit = G / Q0[idx[0][0]]
XsData = bafia.Xs(freq[idx], freq[np.min(idx)] + 11e3, sigman_fit)

sigma1, sigma2, sigma = bafia.sigmaRX(RsData, XsData, freq[np.min(idx)])
sigma1T, sigma2T, sigmaT = bafia.sigmaTrunin(RsData, XsData, sigman_fit)

# === Plot: Frequency vs Temperature ===
plt.figure()
plt.plot(temp, freq, '.-', label='Corrected')
plt.plot(temp, freq1, 'o', label='Raw')
plt.xlabel('Temperature [K]')
plt.ylabel('Frequency [Hz]')
plt.xlim([4, 13])
plt.legend()

# === Plot: Q0 vs Temperature ===
plt.figure()
plt.semilogy(temp, Q0, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('$Q_0$')
plt.xlim([4, 13])

# === Plot: Pressure vs Temperature ===
plt.figure()
plt.semilogy(temp, MKS1000 - 1000, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [mbar]')
plt.xlim([4, 13])

# === Superconductor Conductivity Model ===
freqS = np.max(freq)
tempS = np.linspace(2, Tc - 1e-3, 1000)
mySc = sc.SCconductivity(Tc, freqS, Gamma, tempS, sigman)
Q = mySc.Q()
deltaf = mySc.deltaf()
ZsS = mySc.Zs()
s1S, s2S = mySc.sigma()

# === Delta Lambda ===
idx1 = np.where((temp <= Tc - 0.3) & (temp >= 0.1 * 9.2))
deltaL = bafia.deltaLambda(freq[idx1], temp[idx1])

# Plot Delta Lambda
plt.figure()
plt.plot(temp[idx1], deltaL * 1e10, '.-', label='Data')
plt.xlabel("Temperature [K]")
plt.ylabel("$\Delta \lambda [\mathring{\mathrm{A}}]$")
plt.grid(True)

# Initial guess curve
plt.plot(temp[idx1], bafia.deltaLFit(temp[idx1], Tc - 0.2, 1000, 1, 1, 0), label='Initial Guess')

# === Fit Delta Lambda ===
gmodel = Model(bafia.deltaLFit)
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
plt.show()

# === Frequency and Q0 Comparison ===
fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx], freq[idx] / 1e6)
ax[0].set_xlim([7, 9])
ax1 = ax[0].twinx()
ax1.plot(tempS, deltaf / 1e6, '-.')
ax1.set_xlim([7, 9])
ax[0].set_ylabel("Frequency [MHz]")

ax[1].semilogy(temp[idx], Q0[idx])
ax[1].semilogy(tempS, Q, '-')
ax[1].set_ylabel(r"$Q_0$")
ax[1].set_xlabel("Temperature [K]")

# === Rs and Xs Plot ===
fig, ax = plt.subplots(2, 1)
ax[0].semilogy(temp[idx], RsData)
ax[0].semilogy(tempS, np.real(ZsS))
ax[0].set_ylabel(r"$R_s$ [$\Omega$]")

ax[1].plot(temp[idx], XsData / 1e-3)
ax2 = ax[1].twinx()
ax2.plot(tempS, np.imag(ZsS) / 1e-3, 'c')
ax[1].set_ylabel(r"$X_s$ [m$\Omega$]")
ax[1].set_xlabel("Temperature [K]")

# === Conductivity Plots ===
fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx] / Tc, sigma1 / 2e8, label='RX')
ax[0].plot(temp[idx] / Tc, sigma1T, label='Trunin')
ax02 = ax[0].twinx()
ax02.plot(tempS / Tc, s1S, linestyle='--', label='MB')
ax[0].set_ylabel(r"$\sigma_1$")
ax[0].legend()

ax[1].plot(temp[idx] / Tc, sigma2 / 2e8, label='RX')
ax[1].plot(temp[idx] / Tc, sigma2T, label='Trunin')
ax12 = ax[1].twinx()
ax12.plot(tempS / Tc, s2S, linestyle='--', label='MB')
ax[1].set_ylabel(r"$\sigma_2$")
ax[1].set_xlabel("Temperature [K]")
ax[1].legend()

plt.tight_layout()
plt.show()
