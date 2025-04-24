# `analysis_Helper` class

The `analysis_Helper` class provides tools for analyzing superconducting resonator data. 
It includes methods for calculating surface impedance, complex conductivity, and changes in penetration depth based on resonator measurements.

## Initialization
Simply import the `analysis_Helper.py` and create an istance of the class

Through this you will access the different methods


```
helper = analysis_Helper(G, save=save_plots)
```

- `save_plots` flags is passed to the `_handle_plot` function to take care of plots/save

- `G`: Geometric factor of the resonator (default: 150).

## Methods

### `_handle_plot(self, filename)`

This is a utility function to show or save generated plots.

```
def _handle_plot(self, filename):
    if self.save:
        full_path = os.path.join(self.save_dir, filename)
        plt.savefig(full_path)
        print(f"Saved plot to {full_path}")
        plt.close()
    else:
        plt.show()
```

### `Rs(Q)`

Calculates the surface resistance $R_s$ from the quality factor $Q$.

$$R_s = \frac{G}{Q}$$

### `Xs(f, f0, X0)`

Calculates the surface reactance $X_s$ based on a frequency shift from a reference.

- `f`: Measured frequency  
- `f0`: Reference frequency  
- `X0`: Reference surface reactance

$$X_s = -2G \cdot \frac{f - f_0}{f_0} + X_0$$

> [!TIP]
>
> Refer to: [arXiv:cond-mat/0110109](https://arxiv.org/abs/cond-mat/0110109)

### `sigmaRX(Rs, Xs, freq0)`

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

> [!TIP]
>
> Refer to: Electrodynamics of Solids: Optical Properties of Electrons in Matter. Cambridge University Press 2002

### `sigmaTrunin(Rs, Xs, Rn)`

Estimates complex conductivity using the Trunin approximation.

- `Rs`: Surface resistance  
- `Xs`: Surface reactance  
- `Rn`: Normal-state resistance

$$\sigma_1 = \frac{4 R_n^2 R_s X_s}{(R_s^2 + X_s^2)^2}$$
$$\sigma_2 = \frac{2 R_n^2 (X_s^2 - R_s^2)}{(R_s^2 + X_s^2)^2}$$

**Returns**:  
- `sigma1`: Real part  
- `sigma2`: Imaginary part  
- `sigma`: Complex conductivity

> [!TIP]
>
> Refer to: [Trunin approximation model](http://www.issp.ac.ru/lek/trunin/art60E.pdf).

### `deltaLambda(freq, temp, G=192)`

Calculates the change in London penetration depth from frequency shift measurements.

- `freq`: Array of measured frequencies  
- `temp`: Corresponding temperatures  
- `G`: Geometry factor (default: 192)

$$\Delta \lambda(T) = -\frac{G (f - f_0)}{\pi \mu_0 f_0^2}$$

where $f_0$ is the frequency at the base temperature (e.g., $T \leq 5\,K$).

> [!TIP]
>
> Refer to: [Brorson et al](https://arxiv.org/abs/cond-mat/9311027)

### `deltaLFit(temp, Tc, lLondon, l, eps, l0)`

Fits the change in penetration depth $\Delta \lambda(T)$ using a standard superconducting model.

- `temp`: Temperature array  
- `Tc`: Critical temperature  
- `lLondon`: Zero-temperature London penetration depth  
- `l`: Thickness of the superconducting film  
- `eps`: Dielectric constant factor  
- `l0`: Reference penetration depth offset

$$\Delta \lambda(T) = \lambda_L \cdot \sqrt{1 + \frac{\varepsilon}{l}} \cdot \frac{1}{\sqrt{1 - \left(\frac{T}{T_c}\right)^4}} - \lambda_0$$

> [!TIP]
>
> Refer to: Ciovati's [SUPERFIT](https://www.researchgate.net/publication/255216727_SUPERFIT_a_Computer_Code_to_Fit_Surface_Resistance_and_Penetration_Depth_of_a_Superconductor) 

## Example Visualizations:
The plot functions produce the following plots
- Temperature vs Frequency plot.
- Temperature vs Pressure plot.
- Temperature vs Quality Factor (`Q0`).
- Temperature vs Resistance (`Rs`), Reactance (`Xs`), and Complex Conductivity (`sigma`).
- Curve fitting of the delta lambda model to extract superconductivity parameters.