# bafia_analysis.py

import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameters
import SCconductivity as sc
from analysis_utils import analysis_Helper

# === Parameters ===
fileName = "20250128_FNAL_103.txt"
freqS = 650e6
Tc = 8.7
Gamma = 0.06
sigman = 1 / (152e-9 * 1e-2)
G = 192

# === Instantiate Helper ===
bafia = analysis_Helper(G=G)

# === Load Data ===
data = np.loadtxt(fileName, skiprows=4)
temp = data[:, 1]
freq = data[:, 5]
Q0 = data[:, 6]

# Clean invalid entries
temp = np.delete(temp, np.where(freq == freq.min()))
freq = np.delete(freq, np.where(freq == freq.min()))
Q0 = np.delete(Q0, np.where(freq == freq.min()))

idx = np.where((temp <= 9) & (temp >= 0.1 * 9.2))
RsData = bafia.Rs(Q0[idx])
sigman = G / Q0[idx[0][0]]
XsData = bafia.Xs(freq[idx], freq[np.min(idx)] + 11e3, sigman)

sigma1, sigma2, sigma = bafia.sigmaRX(RsData, XsData, freq[np.min(idx)])
sigma1T, sigma2T, sigmaT = bafia.sigmaTrunin(RsData, XsData, sigman)

# === Plotting frequency & Q0 ===
plt.figure()
plt.plot(temp, freq, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('Frequency [Hz]')
plt.xlim([4, 11])
plt.ylim([6.4982e8, 6.4984e8])

plt.figure()
plt.plot(temp, Q0, '.')
plt.xlabel('Temperature [K]')
plt.ylabel('$Q_0$')
plt.xlim([4, 11])

# === Superconductivity Model ===
tempS = np.linspace(2, Tc - 1e-3, 1000)
mySc = sc.SCconductivity(Tc, freqS, Gamma, tempS, sigman)
Q = mySc.Q()
deltaf = mySc.deltaf()
ZsS = mySc.Zs()
s1S, s2S = mySc.sigma()

# === Delta Lambda ===
idx1 = np.where((temp <= Tc - 0.3) & (temp >= 0.1 * 9.2))
deltaL = bafia.deltaLambda(freq[idx1], temp[idx1])

fig = plt.figure()
plt.plot(temp[idx1], deltaL * 1e10, label="Data")
plt.xlabel("Temperature [K]")
plt.ylabel("$\Delta \lambda [\mathring {\mathrm {A}}]$")
plt.grid(True)

plt.plot(temp[idx1], bafia.deltaLFit(temp[idx1], Tc - 0.2, 1000, 1, 1, 0), label="Initial Guess")

# === Fit delta lambda ===
gmodel = Model(bafia.deltaLFit)
params = Parameters()
params.add('Tc', min=8.6, max=8.9, value=Tc)
params.add('lLondon', min=100, max=2000, value=500)
params.add('l', min=400, max=600, value=600)
params.add('eps', min=0.1, max=600, value=620, vary=False)
params.add('l0', min=100, max=600, value=610)

result = gmodel.fit(deltaL * 1e10, temp=temp[idx1], params=params)

print(result.fit_report())
plt.plot(temp[idx1], result.best_fit, label="Fit")
plt.legend()
plt.show()

# === Additional Plots ===

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx], freq[idx] / 1e6)
ax[0].set_xlim([7, 9])
ax1 = ax[0].twinx()
ax1.plot(tempS, deltaf / 1e6, '-.')
ax1.set_xlim([7, 9])
ax[0].set_ylabel("Frequency [Hz]")

ax[1].semilogy(temp[idx], Q0[idx])
ax[1].semilogy(tempS, Q, '-')
ax[1].set_ylabel(r"$Q_0$")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].semilogy(temp[idx], RsData)
ax[0].semilogy(tempS, np.real(ZsS))
ax[0].set_ylabel(r"$R_s$ [$\Omega$]")

ax[1].plot(temp[idx], XsData / 1e-3)
ax2 = ax[1].twinx()
ax2.plot(tempS, np.imag(ZsS) / 1e-3, 'c')
ax[1].set_ylabel(r"$X_s$ [m$\Omega$]")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx] / Tc, sigma1 / 2e8, label='RX')
ax[0].plot(temp[idx] / Tc, sigma1T, label='Trunin')
ax02 = ax[0].twinx()
ax02.plot(tempS / Tc, s1S, label='Mattis-Bardeen', linestyle='--')
ax[0].set_ylabel(r"$\sigma_1$")
ax[0].legend()

ax[1].plot(temp[idx] / Tc, sigma2 / 2e8, label='RX')
ax[1].plot(temp[idx] / Tc, sigma2T, label='Trunin')
ax12 = ax[1].twinx()
ax12.plot(tempS / Tc, s2S, label='Mattis-Bardeen', linestyle='--')
ax[1].set_ylabel(r"$\sigma_2$")
ax[1].set_xlabel("Temperature [K]")
ax[1].legend()

plt.tight_layout()
plt.show()
