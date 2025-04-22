# Superconductivity Analysis Scripts


This repository of python scripts for analyzing superconductivity data.
These scripts use experimental data to model and visualize the complex conductivity behavior.
The focus is on temperature and frequency-dependent properties.

> [!WARNING]
>
> This is a work in progress and stuff might be wrong. 
> I will try to place references where possible for the different approximation and simulations.

## Overview

There are two main scripts in this repository and two 'utility':

1. **BAFIA (arXiv:2103.10601)**  
2. **FG004 (Sertore 11-3-2024)**

The two main scripts perform basically the same analysis on different data:

- **Frequency-dependent Analysis:** It computes the frequency dependence of resistance (`Rs`), reactance (`Xs`), and complex conductivity (`sigmaRX`). In the second there is an additional correction to the frequency using the pressure information (`MKS1000`)
- **Modeling and Visualization:** The script utilizes the `SCconductivity` class for modeling the superconducting material's behavior and visualizes important properties such as frequency, temperature, and conductivity.
- **Curve Fitting:** The script fits a delta lambda model to the data and optimizes the parameters using `lmfit`.
  
3. **SCconductivity**
4. **analysis_utils**

The first is the script containing the tools to model and simulate the behaviour of the system while the second is a collection of tools for analyzing and representing the data.

## `analysis_Helper` Class

The `analysis_Helper` class provides tools for analyzing superconducting resonator data. 
It includes methods for calculating surface impedance, complex conductivity, and changes in penetration depth based on resonator measurements.

<details>
<summary>Details on this script</summary>

### Initialization
Simply import the `analysis_Helper.py` and create an istance of the class

Through this you will access the different methods.

```python
helper = analysis_Helper(G=150)
```

- `G`: Geometric factor of the resonator (default: 150).

### Methods

#### `Rs(Q)`

Calculates the surface resistance $R_s$ from the quality factor $Q$.

$$R_s = \frac{G}{Q}$$

#### `Xs(f, f0, X0)`

Calculates the surface reactance $X_s$ based on a frequency shift from a reference.

- `f`: Measured frequency  
- `f0`: Reference frequency  
- `X0`: Reference surface reactance

$$X_s = -2G \cdot \frac{f - f_0}{f_0} + X_0$$

#### `sigmaRX(Rs, Xs, freq0)`

Calculates the complex conductivity $\sigma = \sigma_1 + i\sigma_2$ using surface resistance and reactance.

- `Rs`: Surface resistance  
- `Xs`: Surface reactance  
- `freq0`: Reference frequency

$$\sigma = \omega \mu_0 \left( \frac{2 R_s X_s}{(R_s^2 + X_s^2)^2} + i \cdot \frac{X_s^2 - R_s^2}{(R_s^2 + X_s^2)^2} \right)$$

where:  
- $\omega = 2\pi f_0$  
- $\mu_0$: Permeability of free space

**Returns**:  
- `sigma1`: Real part  
- `sigma2`: Imaginary part  
- `sigma`: Complex conductivity

#### `sigmaTrunin(Rs, Xs, Rn)`

Estimates complex conductivity using the [Trunin approximation model](http://www.issp.ac.ru/lek/trunin/art60E.pdf).

- `Rs`: Surface resistance  
- `Xs`: Surface reactance  
- `Rn`: Normal-state resistance

$$\sigma_1 = \frac{4 R_n^2 R_s X_s}{(R_s^2 + X_s^2)^2}$$
$$\sigma_2 = \frac{2 R_n^2 (X_s^2 - R_s^2)}{(R_s^2 + X_s^2)^2}$$

**Returns**:  
- `sigma1`: Real part  
- `sigma2`: Imaginary part  
- `sigma`: Complex conductivity

#### `deltaLambda(freq, temp, G=192)`

Calculates the change in London penetration depth from frequency shift measurements.

- `freq`: Array of measured frequencies  
- `temp`: Corresponding temperatures  
- `G`: Geometry factor (default: 192)

$$\Delta \lambda(T) = -\frac{G (f - f_0)}{\pi \mu_0 f_0^2}$$

where $f_0$ is the frequency at the base temperature (e.g., $T \leq 5\,K$).

#### `deltaLFit(temp, Tc, lLondon, l, eps, l0)`

Fits the change in penetration depth $\Delta \lambda(T)$ using a standard superconducting model.

- `temp`: Temperature array  
- `Tc`: Critical temperature  
- `lLondon`: Zero-temperature London penetration depth  
- `l`: Thickness of the superconducting film  
- `eps`: Dielectric constant factor  
- `l0`: Reference penetration depth offset

$$\Delta \lambda(T) = \lambda_L \cdot \sqrt{1 + \frac{\varepsilon}{l}} \cdot \frac{1}{\sqrt{1 - \left(\frac{T}{T_c}\right)^4}} - \lambda_0$$

### Example Visualizations:
- Temperature vs Frequency plot.
- Temperature vs Pressure plot.
- Temperature vs Quality Factor (`Q0`).
- Temperature vs Resistance (`Rs`), Reactance (`Xs`), and Complex Conductivity (`sigma`).
- Curve fitting of the delta lambda model to extract superconductivity parameters.
</details>

## `SCconductivity` Class

A Python script for calculating and plotting various superconducting properties such as the real and imaginary parts of conductivity, quality factor, frequency shift, and the superconducting energy gap.


<details>
<summary>Details on this script</summary>

### Overview

The script calculates the following:

- **Real and imaginary parts of conductivity (σ₁ and σ₂)** as functions of temperature $T$ and frequency $\omega$.
- **Quality factor (Q)**, which quantifies the resonance sharpness of a system.
- **Frequency shift ($\Delta$f)** due to temperature variations.
- **Superconducting energy gap ($\Delta$)** as a function of temperature.

### Equations

#### Delta Zero ($\Delta_0$)

The energy gap at zero temperature $ \Delta_0 $ is given by:

$$\Delta_0 = 1.764 k_B T_c$$

where $ k_B $ is the Boltzmann constant and $ T_c $ is the critical temperature.

#### Fermi Distribution Function

The Fermi-Dirac distribution function $ f(E, T) $ is used to model the occupancy of energy states:

$$f(E, T) = \frac{1}{1 + e^{\frac{E}{k_B T}}}$$

#### Delta ($\Delta$) as a function of Temperature

The superconducting energy gap $ \Delta(T) $ depends on the temperature as:

$$\Delta(T) = \Delta_0 \sqrt{1 - \frac{T}{T_c}}$$

#### Determinant and Numerical Functions

The determinant function used in the calculation of conductivity is:

$$\text{det}(E, T) = (E + i \Gamma)^2 - \Delta(T)^2$$

where $ \Gamma $ is a phenomenological broadening factor.

The numerators and Green's functions used in the conductivity expressions are:

$$f_{\text{numer}}(E, T) = (1 - 2 f(E + \hbar \omega, T))$$

and

$$g(E, T) = \frac{(E + i \Gamma)((E + i \Gamma) + \hbar \omega) + \Delta(T)^2}{\sqrt{(E + i \Gamma + \hbar \omega)^2 - \Delta(T)^2}}$$

#### Conductivity Calculations

The real and imaginary parts of the conductivity, $ \sigma_1 $ and $ \sigma_2 $, are computed by integrating the respective functions over the energy range.

#### Quality Factor (Q)

The quality factor $ Q $ is related to the conductivities and is given by:

$$Q = \frac{G}{Z_s + R_r}$$

where $G$ is the geometric-specific constant, $ Z_s $ is the impedance, and $ R_r $ is the residual resistance.

#### Frequency Shift (Δf)

The frequency shift due to temperature changes is given by:

$$\Delta f = -\frac{(\text{Im}(Z_s) - \frac{G}{Q}) \cdot f}{2 G}$$

#### Superconducting Energy Gap in MeV

The energy gap $ \Delta(T) $ is converted into units of meV using:

$$\Delta(T) \ [\text{meV}] = \Delta(T) \ [\text{Joules}] \times 6.242 \times 10^{18}$$

### Code Usage

To use the code, simply initialize an instance of the `SCconductivity` class with the desired parameters, such as the critical temperature $ T_c $, frequency $f$, broadening parameter $\Gamma$, and temperature array. The following example demonstrates how to plot various properties:

```python
Tc = 9.2  # Critical temperature in K
freq = 650e6  # Frequency in Hz
temp = np.linspace(1.5, Tc - 1e-3, 1000)  # Temperature range from 1.5 K to Tc
Gamma = 0.0058 * 0  # Broadening factor (for zero broadening)
sigman = 1/(152e-9 * 1e-2)  # Conductivity parameter

sc = SCconductivity(Tc, freq, Gamma, temp, sigman)

# Plot sigma1 and sigma2 vs Temperature
s1, s2 = sc.sigma()
plt.plot(temp/Tc, np.real(s1), label='$\sigma_1$')
plt.plot(temp/Tc, np.real(s2), label='$\sigma_2$')
plt.xlabel('Temperature / Tc')
plt.ylabel('Conductivity')
plt.legend()
plt.show()
```
</details>