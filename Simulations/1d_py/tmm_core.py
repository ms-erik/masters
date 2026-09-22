import numpy as np
from constants import HBAR, M0, EV_TO_JOULE, NM_TO_METER

def calculate_k(E_eV, V_eV, m_eff_ratio):
    """
    Calculate the wave vector k in a given layer.
    """
    E_J = E_eV * EV_TO_JOULE
    V_J = V_eV * EV_TO_JOULE
    m_star = m_eff_ratio * M0
    
    val = 2 * m_star * (E_J - V_J + 0j)
    k = np.sqrt(val) / HBAR
    return k

def propagation_matrix(k, d_nm):
    """
    Propagation matrix for a layer of thickness d.
    """
    d_m = d_nm * NM_TO_METER
    phi = k * d_m
    
    P = np.array([
        [np.exp(1j * phi), 0],
        [0, np.exp(-1j * phi)]
    ], dtype=complex)
    
    return P

def interface_matrix(k1, k2, m1_ratio, m2_ratio):
    """
    Interface matrix applying BenDaniel-Duke boundary conditions.
    """
    beta1 = k1 / m1_ratio
    beta2 = k2 / m2_ratio
    
    if abs(beta2) < 1e-30:
        beta2 = 1e-30 + 1e-30j

    M11 = 0.5 * (1 + beta1/beta2)
    M12 = 0.5 * (1 - beta1/beta2)
    M21 = 0.5 * (1 - beta1/beta2)
    M22 = 0.5 * (1 + beta1/beta2)
    
    I = np.array([
        [M11, M12],
        [M21, M22]
    ], dtype=complex)
    
    return I

def calculate_transfer_matrix(E_eV, layers):
    """
    layers: list of dicts [{'V': eV, 'm': m_ratio, 'd': nm}, ...]
    The first and last layers are considered semi-infinite contacts (d is ignored).
    Returns the total transfer matrix and wavevectors of contacts.
    """
    k_L = calculate_k(E_eV, layers[0]['V'], layers[0]['m'])
    
    M_tot = np.eye(2, dtype=complex)
    
    k_prev = k_L
    m_prev = layers[0]['m']
    
    for i in range(1, len(layers) - 1):
        layer = layers[i]
        k_curr = calculate_k(E_eV, layer['V'], layer['m'])
        
        I_mat = interface_matrix(k_prev, k_curr, m_prev, layer['m'])
        P_mat = propagation_matrix(k_curr, layer['d'])
        
        M_tot = M_tot @ I_mat @ P_mat
        
        k_prev = k_curr
        m_prev = layer['m']
    
    k_R = calculate_k(E_eV, layers[-1]['V'], layers[-1]['m'])
    I_final = interface_matrix(k_prev, k_R, m_prev, layers[-1]['m'])
    
    M_tot = M_tot @ I_final
    
    return M_tot, k_L, k_R

def calculate_transmission(E_eV, layers):
    """
    Calculates Transmission T and Reflection R probabilities.
    """
    M_tot, k_L, k_R = calculate_transfer_matrix(E_eV, layers)
    
    M11 = M_tot[0, 0]
    M21 = M_tot[1, 0]
    
    t = 1.0 / M11
    r = M21 / M11
    
    m_L = layers[0]['m']
    m_R = layers[-1]['m']
    
    v_L = np.real(k_L / m_L)
    v_R = np.real(k_R / m_R)
    
    if v_L == 0:
        return 0.0, 1.0
        
    T = (v_R / v_L) * (np.abs(t)**2)
    R = np.abs(r)**2
    
    return T, R
