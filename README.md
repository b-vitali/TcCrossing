# Superconductivity Analysis Scripts


This repository of python scripts for analyzing superconductivity data.
These scripts use experimental data to model and visualize the complex conductivity behavior.
The focus is on temperature and frequency-dependent properties.

> [!WARNING]
>
> This is a work in progress and stuff might be wrong. 
> I will try to place references where possible for the different approximation and simulations.

## Overview

Aside from the data, markdown, .gitignore and legacy files, there are 4 scripts:

```
📁 TcCrossing
├── 📄 bafia.py
├── 📄 fg004.py
│
├── 📁 utils
│ ├── 🐍 analysis_Helper.py
│ └── 🐍 SCconductivity.py
```

1. **BAFIA (arXiv:2103.10601)**  
2. **FG004 (Sertore 11-3-2024)**

These are the two main scripts perform basically the same analysis on different data:

- **Frequency-dependent Analysis:** It computes the frequency dependence of resistance (`Rs`), reactance (`Xs`), and complex conductivity (`sigmaRX`). In the second there is an additional correction to the frequency using the pressure information (`MKS1000`)
- **Modeling and Visualization:** The script utilizes the `SCconductivity` class for modeling the superconducting material's behavior and visualizes important properties such as frequency, temperature, and conductivity.
- **Curve Fitting:** The script fits a delta lambda model to the data and optimizes the parameters using `lmfit`.
  
3. **SCconductivity**
4. **analysis_utils**

The first is the script containing the tools to model and simulate the behaviour of the system while the second is a collection of tools for analyzing and representing the data.

## `analysis_Helper` Class

The `analysis_Helper` class provides tools for analyzing superconducting resonator data. 
It includes methods for calculating surface impedance, complex conductivity, and changes in penetration depth based on resonator measurements.

📖 Read more: [Details on analysis_Helper](docs/analysis_Helper.md)


## `SCconductivity` Class

A Python script for calculating and plotting various superconducting properties such as the real and imaginary parts of conductivity, quality factor, frequency shift, and the superconducting energy gap.

📖 Read more: [Details on SCconductivity](docs/SCconductivity.md)