# fg004_analysis.py
# FG004 Analysis Script — 11-3-2024

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from lmfit import Model, Parameters
import SCconductivity as sc  # Ensure this module is available in the same directory

# ------------------------- Constants & Functions -------------------------

G = 150

def Rs(Q):
    return G / Q

def Xs(f, f0, X0):
    return -2 * G * (f - f0) / f0 + X0

def sigmaRX(Rs, Xs, freq0):
    omega = 2 * np.pi * freq0
    mu0 = const.mu_0
    sigma = omega * mu0 * (2 * Rs * Xs / (Rs**2 + Xs**2)**2 + (Xs**2 - Rs**2) / (Rs**2 + Xs**2)**2 * 1j)
    return np.real(sigma), np.imag(sigma), sigma

def sigmaTrunin(Rs, Xs, Rn):
    sigma1 = 4 * Rn**2 * Rs * Xs / (Rs**2 + Xs**2)**2
    sigma2 = 2 * Rn**2 * (Xs**2 - Rs**2) / (Rs**2 + Xs**2)**2
    return sigma1, sigma2, (sigma1 + 1j * sigma2)

def deltaLambda(freq, temp, G=192):
    f0 = np.max(freq)
    return -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)

def deltaLFit(temp, Tc, lLondon, l, eps, l0):
    return lLondon * np.sqrt(1 + eps / l) / np.sqrt(1 - (temp / Tc)**4) - l0

# ----------------------------- Load Data -----------------------------

fileName = "FG004_throughTc.txt"
data = np.loadtxt(fileName, skiprows=4)
temp = data[:, 1]
MKS1000 = data[:, 2]
freq1 = data[:, 5]
Q0 = data[:, 6]
freq = freq1 - 750 * (MKS1000 - 1000)

# Clean out invalid entries
mask = freq != freq.min()
temp, freq, freq1, MKS1000, Q0 = temp[mask], freq[mask], freq1[mask], MKS1000[mask], Q0[mask]

# --------------------------- Derived Quantities ---------------------------

idx = np.where((temp <= 9) & (temp >= 0.1 * 9.2))
RsData = Rs(Q0[idx])
sigman = G / Q0[idx[0][0]]
print(f"σₙ = {sigman:.4e}")

XsData = Xs(freq[idx], freq[np.min(idx)] + 11e3, sigman)
sigma1, sigma2, sigma = sigmaRX(RsData, XsData, freq[np.min(idx)])
sigma1T, sigma2T, sigmaT = sigmaTrunin(RsData, XsData, sigman)

# ------------------------------ Raw Plots ------------------------------

plt.figure()
plt.plot(temp, freq, '.-', label="Corrected Freq")
plt.plot(temp, freq1, 'o', label="Raw Freq")
plt.xlabel('Temperature [K]')
plt.ylabel('Frequency [Hz]')
plt.xlim([4, 13])
plt.legend()

plt.figure()
plt.semilogy(temp, Q0, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('$Q_0$')
plt.xlim([4, 13])

plt.figure()
plt.semilogy(temp, MKS1000 - 1000, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [mbar]')
plt.xlim([4, 13])

# ----------------------------- SC Conductivity -----------------------------

freqS = np.max(freq)
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

# ---------------------------- Δλ Fit Section ----------------------------

idx1 = np.where((temp <= Tc - 0.3) & (temp >= 0.1 * 9.2))
deltaL = deltaLambda(freq[idx1], temp[idx1])

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
