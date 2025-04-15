# Superconductivity Analysis Scripts

This repository of python scripts for analyzing superconductivity data.
These scripts use experimental data to model and visualize the complex conductivity behavior.
The focus is on temperature and frequency-dependent properties.

## Overview

There are two main scripts in this repository:

1. **BAFIA (arXiv:2103.10601)**  
2. **FG004 (11-3-2024)**

### 1. **BAFIA (arXiv:2103.10601)**

This script analyzes data from the BAFIA article. 
The script includes:

- **Frequency-dependent Analysis:** It computes the frequency dependence of resistance (`Rs`), reactance (`Xs`), and complex conductivity (`sigmaRX`).
- **Modeling and Visualization:** The script utilizes the `SCconductivity` class for modeling the superconducting material's behavior and visualizes important properties such as frequency, temperature, and conductivity.
- **Curve Fitting:** The script fits a delta lambda model to the data and optimizes the parameters using `lmfit`.
  
Key calculations include:
- Resistance as a function of the quality factor (`Rs`).
- Complex conductivity calculation (`sigmaRX`).
- Temperature dependence of the superconducting material.
- Visualization of the temperature-frequency relationship, quality factor (`Q0`), and resistance/ reactance.

### 2. **FG004 (11-3-2024)**

The FG004 script is similar to the BAFIA script but uses different experimental data. 
This script also calculates and visualizes the superconducting material's properties, with an emphasis on pressure and temperature dependence. 

Key calculations include:
- The same functions for `Rs` and `Xs` as in BAFIA, but with additional handling of the pressure data (`MKS1000`).
- A similar model for complex conductivity (`sigmaRX`).
- Frequency adjustments for pressure-dependent data.
- Visualization of the pressure, temperature, frequency, quality factor (`Q0`), and superconductivity properties.

### Common Functions:
Both scripts use the following functions:
- `Rs(Q)`: Calculates the resistance from the quality factor (`Q`).
- `Xs(f, f0, X0)`: Calculates the reactance as a function of frequency, with tuning to a reference value.
- `sigmaRX(Rs, Xs, freq0)`: Computes the real and imaginary parts of the complex conductivity.
- `sigmaTrunin(Rs, Xs, Rn)`: Computes the conductivity at a transition point.

### Example Visualizations:
- Temperature vs Frequency plot.
- Temperature vs Quality Factor (`Q0`).
- Temperature vs Resistance (`Rs`), Reactance (`Xs`), and Complex Conductivity (`sigma`).
- Curve fitting of the delta lambda model to extract superconductivity parameters.