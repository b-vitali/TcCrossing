# bafia_analysis.py

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from lmfit import Model, Parameters
import SCconductivity as sc  # Ensure this module is available in the same directory

# ---------------------- Constants & Setup ----------------------
G = 150
fileName = "20250128_FNAL_103.txt"

# ---------------------- Functions ----------------------
def Rs(Q):
    return G/Q

def Xs(f, f0, X0):
    return -2*G*(f - f0)/f0 + X0

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
    f0 = freq[np.where(temp <= 5.0)[0][0]]
    dL = -G * (freq - f0) / (const.pi * const.mu_0 * f0**2)
    return dL

def deltaLFit(temp, Tc, lLondon, l, eps, l0):
    dl = lLondon * np.sqrt(1 + eps/l) * 1 / np.sqrt(1 - (temp/Tc)**4) - l0
    return dl

# ---------------------- Load Data ----------------------
data = np.loadtxt(fileName, skiprows=4)
temp = data[:, 1]
freq = data[:, 5]
Q0 = data[:, 6]

# Clean and filter data
temp = np.delete(temp, np.where(freq == freq.min()))
freq = np.delete(freq, np.where(freq == freq.min()))
Q0 = np.delete(Q0, np.where(freq == freq.min()))

idx = np.where((temp <= 9) & (temp >= 0.1 * 9.2))
RsData = Rs(Q0[idx])
sigman = G / Q0[idx[0][0]]
XsData = Xs(freq[idx], freq[np.min(idx)] + 11e3, sigman)

# ---------------------- Sigma Calculation ----------------------
sigma1, sigma2, sigma = sigmaRX(RsData, XsData, freq[np.min(idx)])
sigma1T, sigma2T, sigmaT = sigmaTrunin(RsData, XsData, sigman)

# ---------------------- SC Conductivity ----------------------
freqS = 650e6
Tc = 8.7
tempS = np.linspace(2, Tc - 1e-3, 1000)
Gamma = 0.06
sigman = 1 / (152e-9 * 1e-2)
G = 192
Rr = 4e-9

mySc = sc.SCconductivity(Tc, freqS, Gamma, tempS, sigman)
Q = mySc.Q()
deltaf = mySc.deltaf()
ZsS = mySc.Zs()
s1S, s2S = mySc.sigma()

# ---------------------- deltaLambda ----------------------
idx1 = np.where((temp <= Tc - 0.3) & (temp >= 0.1 * 9.2))
deltaL = deltaLambda(freq[idx1], temp[idx1])

# ---------------------- Fitting ----------------------
gmodel = Model(deltaLFit)
params = Parameters()
params.add('Tc', min=8.6, max=8.9, value=Tc)
params.add('lLondon', min=100, max=2000, value=500)
params.add('l', min=400, max=600, value=600)
params.add('eps', min=0.1, max=600, value=620, vary=False)
params.add('l0', min=100, max=600, value=610)

tEmp = temp[idx1]
result = gmodel.fit(deltaL * 1e10, temp=tEmp, params=params)

# ---------------------- Plotting ----------------------
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

fig = plt.figure()
plt.plot(temp[idx1], deltaL * 1e10)
plt.plot(temp[idx1], result.best_fit)
plt.xlabel("Temperature [K]")
plt.ylabel("$\Delta \lambda [\mathring {\mathrm {A}}]$")
plt.grid(True)
print(result.fit_report())

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx], freq[idx]/1e6)
ax[0].set_xlim([7, 9])
ax1 = ax[0].twinx()
ax1.plot(tempS, deltaf/1e6, '-.')
ax1.set_xlim([7, 9])
ax[0].set_ylabel("Frequency [Hz]")

ax[1].semilogy(temp[idx], Q0[idx])
ax[1].semilogy(tempS, Q, '-')
ax[1].set_ylabel("Q$_0$")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].semilogy(temp[idx], RsData)
ax[0].semilogy(tempS, np.real(ZsS))
ax[0].set_ylabel("R$_s$ [$\Omega$]")

ax[1].plot(temp[idx], XsData/1e-3)
ax1 = ax[1].twinx()
ax1.plot(tempS, np.imag(ZsS)/1e-3, 'c')
ax[1].set_ylabel("X$_s$ [m$\Omega$]")
ax[1].set_xlabel("Temperature [K]")

fig, ax = plt.subplots(2, 1)
ax[0].plot(temp[idx]/Tc, sigma1/2e8, label="RX")
ax[0].plot(temp[idx]/Tc, sigma1T, label="Trunin")
ax0 = ax[0].twinx()
ax0.plot(tempS/Tc, s1S, label="Mattis-Bardeen", linestyle="--")
ax[0].set_ylabel("$\sigma_1$")

ax[1].plot(temp[idx]/Tc, sigma2/2e8)
ax[1].plot(temp[idx]/Tc, sigma2T)
ax1 = ax[1].twinx()
ax1.plot(tempS/Tc, s2S)
ax[1].set_ylabel("$\sigma_2$")
ax[1].set_xlabel("Temperature [K]")

plt.tight_layout()
plt.show()
