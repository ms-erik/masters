import numpy as np
import matplotlib.pyplot as plt
from tmm_core import calculate_transmission
from device_builder import build_arrw_structure, build_single_barrier

def test_single_barrier():
    print("Testing Single Barrier...")
    layers = build_single_barrier()
    
    energies = np.linspace(0.01, 0.5, 500)
    T_list = []
    
    for E in energies:
        T, R = calculate_transmission(E, layers)
        T_list.append(T)
        assert np.isclose(T + R, 1.0, atol=1e-6), f"Unitarity failed at E={E}: T+R={T+R}"
        
    plt.figure(figsize=(8, 4))
    plt.plot(energies, T_list, label='Transmission')
    plt.axvline(x=0.3, color='r', linestyle='--', label='Barrier Height (0.3 eV)')
    plt.title('Single Barrier Transmission')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Transmission Probability T(E)')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.savefig('single_barrier.png')
    print("Saved 'single_barrier.png'")

def test_arrw():
    print("Testing ARR≡W Structure...")
    layers = build_arrw_structure(num_reflectors_per_side=3)
    
    energies = np.linspace(0.201, 0.5, 2000)
    T_list = []
    
    for E in energies:
        T, R = calculate_transmission(E, layers)
        T_list.append(T)
        
    plt.figure(figsize=(8, 4))
    plt.plot(energies, T_list)
    plt.title('ARR≡W Transmission (Mock parameters)')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Transmission Probability T(E)')
    plt.grid(True)
    plt.savefig('arrw_transmission.png')
    print("Saved 'arrw_transmission.png'")

if __name__ == "__main__":
    test_single_barrier()
    test_arrw()
