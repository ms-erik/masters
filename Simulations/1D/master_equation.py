import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

def add_rate(K, to_idx, from_idx, rate):
    """Adiciona a taxa de conversão"""
    K[to_idx  , from_idx] += rate
    K[from_idx, from_idx] -= rate

def build_rate_matrix():
    """Transições - Símbolo / Taxa (10-6)Hz"""
    
    k_ex     = 60e6       # Excitação do laser
    k_rad    = 77e6       # Decaimento radiativo
    k_isc0   = 4e6        # ISC do E0
    k_isc1   = 45e6       # ISC do E1
    k_s_to_0 = 3.5e6      # return to g0
    k_s_to_1 = 0.7e6      # return to g1
    k_ion_0  = 1e6        # turn into NV0
    k_ion_1  = 5e6
    k_ion_s  = 1e6
    k_rec_0  = 1e6        # recapture electron
    k_rec_1  = 0.2e6

    # initialize a 6x6 matrix
    k = np.zeros((6,6))

    # add values to matrix
    # excitations
    add_rate(k, 2, 0, k_ex) # E0 ok
    add_rate(k, 3, 1, k_ex) # E1 ok

    # G0
    add_rate(k, 0, 2, k_rad) # decay
    add_rate(k, 0, 4, k_s_to_0) # return from singleto
    add_rate(k, 0, 5, k_rec_0) # recapture electron

    # G1
    add_rate(k, 1, 3, k_rad)  #decay
    add_rate(k, 1, 4, k_s_to_1) # return from singleto
    add_rate(k, 1, 5, k_rec_1) # recapture electron

    # Singleto
    add_rate(k, 4, 2, k_isc0) # isc 0
    add_rate(k, 4, 3, k_isc1) # isc 1

    # lost electron
    add_rate(k, 5, 2, k_ion_0) # e0 -> c0
    add_rate(k, 5, 3, k_ion_1) # e1 -> c0
    add_rate(k, 5, 4, k_ion_s) # s -> c0

    return k

def solve_time_dynamics(time, p0_initial):
    k = build_rate_matrix()
    p_nv_minus_history = []

    for t in time:
        p_t = expm(k*t) @ p0_initial

        prob_minus = np.sum(p_t[0:5])
        p_nv_minus_history.append(prob_minus)

    return np.array(p_nv_minus_history)


t = np.linspace(0, 3e-6, 100)
p0 = [1, 0, 0, 0, 0, 0]
p0_linha = [0, 1, 0, 0, 0, 0]

s0 = solve_time_dynamics(t, p0)
s0_linha = solve_time_dynamics(t, p0_linha)

plt.figure(figsize=(8,5))
plt.plot(t, s0, color='red', label = "ms=0")
plt.plot(t, s0_linha, color='blue', label='ms=1')
plt.title("charge probability")
plt.xlabel("time")
plt.ylabel("nv minus")
plt.show()

