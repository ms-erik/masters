from constants import *

def k_nm(E_meV, V_eV, m_eff_ratio): 

    E_J = E_meV * q
    V_J = V_eV * q
    m_eff = m_eff_ratio * m0

    aux = 2 * m_eff * (E_J - V_J + 0j)
    k = np.sqrt(aux) / hbar
    k_nm = k * m_to_nm 
    return k_nm

