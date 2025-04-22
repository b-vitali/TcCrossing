import numpy as np
from scipy.integrate  import quad
import scipy.constants as const
import matplotlib.pyplot as plt

class SCconductivity:

    hbar = const.hbar
    k = const.Boltzmann

    def __init__(self, Tc, freq, Gamma, temp, sigman =  1/(152e-9*1e-2)):
        self.Tc = Tc
        self.Delta0 = self.Delta0F()
        self.sigman = sigman
        self.freq = freq
        self.omega = 2*const.pi*freq
        self.Gamma = Gamma*self.Delta0
        self.temp = temp


    def Delta0F(self):
        return 1.764*self.k*self.Tc

    def f(self, E, T):
        return 1/(1+np.exp(E/(self.k*T)))

    def Delta(self, T):
        mio = self.Delta0*1.74*np.sqrt(1-T/self.Tc)
        bafia = self.Delta0*np.sqrt(np.cos(1.571*(T/self.Tc)**2))
        return bafia

    def determ(self, E, T):
        return (E+1j*self.Gamma)**2-self.Delta(T)**2

    def numer(self, E, T):
        return (1-2*self.f(E+self.hbar*self.omega, T))

    def g(self, E, T):
        num = (E+1j*self.Gamma)*((E+1j*self.Gamma)+self.hbar*self.omega)+self.Delta(T)**2
        den = np.sqrt(((E+1j*self.Gamma)+self.hbar*self.omega)**2-self.Delta(T)**2)
        return num/den

    def sigma2n(self, T):
        
        def func(E, T):          
            func2 = (1-2*self.f(E+self.hbar*self.omega,T))*self.g(E, T)/np.sqrt(self.Delta(T)**2-(E+1j*self.Gamma)**2)
            
            return func2
        
        def func_real(E, T):
            return np.real(func(E, T))
        
        def func_imag(E, T):
            return np.imag(func(E, T))

        func_real_int, err_real_int = quad(func_real, self.Delta(T)-self.hbar*self.omega, self.Delta(T), (T))
        func_imag_int, err_imag_int = quad(func_imag, self.Delta(T)-self.hbar*self.omega, self.Delta(T), (T))
        
        sigma2val = 1/(self.hbar*self.omega)*(func_real_int+1j*func_imag_int)

        # print(sigma2val)

        return sigma2val

    def sigma1n(self, T):
        
        def func1(E, T):
            return ((self.f(E, T)-self.f(E+self.hbar*self.omega, T))*self.g(E, T))/np.sqrt(self.determ(E, T))

        def func2(E, T):
            return self.numer(E,  T)*self.g(E, T)/np.sqrt(self.determ(E, T))

        def func1_real(E, T):
            return np.real(func1(E, T))
        
        def func1_imag(E, T):
            return np.imag(func1(E, T))
        
        def func2_real(E, T):
            return np.real(func2(E, T))
        
        def func2_imag(E, T):
            return np.imag(func2(E, T))

        func1_real_int, err1_real_int = quad(func1_real, self.Delta(T), 100*self.Delta(T), (T))
        func1_imag_int, err1_imag_int = quad(func1_imag, self.Delta(T), 100*self.Delta(T), (T))
        func2_real_int, err2_real_int = quad(func2_real, self.Delta(T)-self.hbar*self.omega, -self.Delta(T), (T))
        func2_imag_int, err2_imag_int = quad(func2_imag, self.Delta(T)-self.hbar*self.omega, -self.Delta(T), (T))

        sigma1v = 2/(self.hbar*self.omega)*(func1_real_int+1j*func1_imag_int) #+1/(hbar*omega)*(func2_real_int+1j*func2_imag_int)

        return sigma1v

    def sigma(self):
        sigman=self.sigman
        s1=[self.sigma1n(t)*sigman for t in self.temp]
        s2=[self.sigma2n(t)*sigman for t in self.temp]
        # print(np.real(s1))
        # print(np.real(s2))

        return s1, s2

    def Zs(self):
        s1, s2 = self.sigma()
        Zs = np.sqrt(1j*self.omega*const.mu_0/(np.real(s1)-1j*np.real(s2)))
        return Zs

    def Q(self, G= 192, Rr=4e-9):
        Q = G/(np.real(self.Zs())+Rr)
        return Q
    
    def deltaf(self, G = 192, Q = 4e4):
        df = ((np.imag(self.Zs())-G/Q)*self.freq)/(-2*G)
        return df
    
def joules_to_mev(energy_joules):
    return energy_joules * 6.242e18  # Conversion factor from Joules to meV

if __name__ == '__main__':
    Tc = 9.2
    freq = 650e6
    temp = np.linspace(1.5, Tc-1e-3, 1000)
    Gamma = 0.0058*0
    sigman =  1/(152e-9*1e-2)
    G = 192
    Rr = 4e-9

    sc = SCconductivity(Tc, freq, Gamma, temp, sigman)

    # First Plot: σ₁ and σ₂ vs Temperature/Tc
    fig, ax = plt.subplots()
    s1, s2 = sc.sigma()
    p1, = ax.plot(temp/Tc, np.real(s1)/sigman, '-r', label="$\sigma_1$")
    ax1 = ax.twinx()
    p2, = ax1.plot(temp/Tc, np.real(s2)/sigman, '-b', label="$\sigma_2$")
    
    # Add labels, title, grid, and legend
    ax.set_xlabel('Temperature / Tc', fontsize=12)
    ax.set_ylabel("$\sigma_1$", color='tab:red', fontsize=12)
    ax1.set_ylabel("$\sigma_2$", color='tab:blue', fontsize=12)
    ax.set_title('R/Im Conductivities vs Temperature (normalized to Tc)', fontsize=14)
    ax.tick_params(axis='y', labelcolor='tab:red')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax.grid(True)
    
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax.transAxes)
    plt.tight_layout()  # Adjust layout to avoid overlap
    plt.show()

    # Second Plot: Q vs Temperature
    Q = sc.Q(G, Rr)
    plt.figure()
    plt.semilogy(temp, Q, '-g')
    plt.xlabel('Temperature [K]', fontsize=12)
    plt.ylabel('$Q_0$', fontsize=12)
    plt.title('Quality Factor (Q) vs Temperature', fontsize=14)
    plt.grid(True, which='both', linestyle='--', color='gray', alpha=0.5)
    plt.tight_layout()
    plt.show()

    # Third Plot: Δf vs Temperature
    deltaf = sc.deltaf()
    plt.figure()
    plt.plot(temp, deltaf, '-m')
    plt.xlabel('Temperature [K]', fontsize=12)
    plt.ylabel('$\Delta f$ [Hz]', fontsize=12)
    plt.title('Frequency Shift (Δf) vs Temperature', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Fourth Plot: Δ vs Temperature
    plt.figure()
    plt.plot(temp, joules_to_mev(sc.Delta(temp)), '-c')
    plt.xlabel('Temperature [K]', fontsize=12)
    plt.ylabel('$\Delta$ [meV]', fontsize=12)
    plt.title('Superconducting Energy Gap (Δ) vs Temperature', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.show()