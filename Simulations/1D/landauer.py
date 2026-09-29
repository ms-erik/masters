import numpy as np
from constants import E_F, V_bias, Temperature_K, EV_TO_JOULE, q, h, K_BT

def fermi_dirac(E_meV, mu_meV, T_kelvin):
    num = (E_meV - mu_meV)* 1e-3 * EV_TO_JOULE
    den = (K_BT * T_kelvin)
    f_e = 1.0 / (1.0 + np.exp(num / den))
    return f_e

def calculate_current(energies_meV, transmission_array):
    mu_l = E_F + V_bias/2
    mu_r = E_F - V_bias/2
    
    janela_fermi = fermi_dirac(energies_meV, mu_l, Temperature_K)-fermi_dirac(energies_meV, mu_r, Temperature_K)

    integrando = transmission_array * janela_fermi
    energies_J = energies_meV * 1e-3 * EV_TO_JOULE

    integral = np.trapezoid(y = integrando, x=energies_J)

    prefator = (2.0 * q)/ h

    I_amperes = prefator * integral
    return I_amperes

