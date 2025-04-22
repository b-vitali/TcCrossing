## `SCconductivity` Class

A Python script for calculating and plotting various superconducting properties such as the real and imaginary parts of conductivity, quality factor, frequency shift, and the superconducting energy gap.

### Overview

The script calculates the following:

- **Real and imaginary parts of conductivity ($\sigma_1$ and $\sigma_2$)** as functions of temperature $T$ and frequency $\omega$.
- **Quality factor (Q)**, which quantifies the resonance sharpness of a system.
- **Frequency shift ($\Delta f$)** due to temperature variations.
- **Superconducting energy gap ($\Delta$)** as a function of temperature.

### Equations

#### Delta Zero ($\Delta_0$)

The energy gap at zero temperature $\Delta_0$ is given by:

$$\Delta_0 = 1.764 k_B T_c$$

where $k_B$ is the Boltzmann constant and $T_c$ is the critical temperature.

#### Fermi Distribution Function

The Fermi-Dirac distribution function $f(E, T)$ is used to model the occupancy of energy states:

$$f(E, T) = \frac{1}{1 + e^{\frac{E}{k_B T}}}$$

#### Delta ($\Delta$) as a function of Temperature

The superconducting energy gap $\Delta(T)$ depends on the temperature as:

$$\Delta(T) = \Delta_0 \sqrt{1 - \frac{T}{T_c}}$$

> [!CAUTION]
>
> This is an approx of the BCS equation, can I use it?
> 
> $$\Delta(T) \approx \Delta(0) \cdot \tanh\left[ 1.74 \sqrt{\frac{T_c}{T} - 1} \right]$$


#### Determinant and Numerical Functions

The determinant function used in the calculation of conductivity is:

$$\text{det}(E, T) = (E + i \Gamma)^2 - \Delta(T)^2$$

where $\Gamma$ is a phenomenological broadening factor.

The numerators and Green's functions used in the conductivity expressions are:

$$f_{\text{numer}}(E, T) = (1 - 2 f(E + \hbar \omega, T))$$

and

$$g(E, T) = \frac{(E + i \Gamma)((E + i \Gamma) + \hbar \omega) + \Delta(T)^2}{\sqrt{(E + i \Gamma + \hbar \omega)^2 - \Delta(T)^2}}$$

> [!TIP]
>
> Refer to: Gor’kov equations, Green’s function formalism for superconductors and Mattis-Bardeen

#### Conductivity Calculations

The real and imaginary parts of the conductivity, $\sigma_1$ and $\sigma_2$, are computed by integrating the respective functions over the energy range.

> [!TIP]
>
> Refer to: Mattis-Bardeen or Nam’s generalization of the BCS

#### Quality Factor (Q)

The quality factor $Q$ is related to the conductivities and is given by:

$$Q = \frac{G}{Z_s + R_r}$$

where $G$ is the geometric-specific constant, $Z_s$ is the impedance, and $R_r$ is the residual resistance.

#### Frequency Shift (Δf)

The frequency shift due to temperature changes is given by:

$$\Delta f = -\frac{(\text{Im}(Z_s) - \frac{G}{Q}) \cdot f}{2 G}$$

#### Superconducting Energy Gap in MeV

The energy gap $\Delta(T)$ is converted into units of meV using:

$$\Delta(T) \ [\text{meV}] = \Delta(T) \ [\text{Joules}] \times 6.242 \times 10^{18}$$

### Code Usage

To use the code, simply initialize an instance of the `SCconductivity` class with the desired parameters, such as the critical temperature $T_c$, frequency $f$, broadening parameter $\Gamma$, and temperature array. The following example demonstrates how to plot various properties:

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