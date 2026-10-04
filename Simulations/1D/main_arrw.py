import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from device_builder import DeviceParams
from tm_transport import identity_scatter, interface_scatter, propagation_scatter, star_product
from materials import get_material_properties
from landauer import calculate_current, fermi_dirac
from constants import E_F, V_bias, Temperature_K

def calculate_transmission(V_base, m_eff, transmission):
    for E in energies_meV:

        S_total = identity_scatter()

        k_prev = params.calculate_k(E, V_base[0], m_eff[0])
        m_prev = m_eff[0]

        for i in range(1, len(z)):

            k_curr = params.calculate_k(E, V_base[i], m_eff[i])
            m_curr = m_eff[i]

            # interface
            S_interface = interface_scatter(k_prev, k_curr, m_prev, m_curr)
            S_total = star_product(S_total, S_interface)

            S_prop = propagation_scatter(k_curr, dz)
            S_total = star_product(S_total, S_prop)

            k_prev = k_curr
            m_prev = m_curr

        # Transmission probability
        k_l = params.calculate_k(E, V_base[0], m_eff[0])
        k_r = params.calculate_k(E, V_base[-1], m_eff[-1])

        v_l = np.real(k_l / m_eff[0])
        v_r = np.real(k_r / m_eff[-1])

        #  T = (v_r / v_l) * n |tlr|²
        if v_l == 0:
            T = 0.0
        else:
            T = (v_r / v_l) * (np.abs(S_total.tLR)**2)

        transmission.append(T)

# constructor instance
params = DeviceParams(targetEnergy_meV=150.0) 
params.compute_thickness()

# layers' thickness
top_thickness = 3.0
buffer_thickness = 5.5
film_thickness = params.d_high_nm
inner_buffer_thickness = 6.0
central_guide_thickness = params.d_guide_nm
bottom_thickness = 8.0

# define list to keep thickness and potentials
thickness  = []
potentials = []
masses = []

# top access
thickness.append(top_thickness)
potentials.append(params.v_high)
masses.append(params.m_high)

for _ in range(5):
    thickness.append(buffer_thickness)
    potentials.append(params.v_buffer)
    masses.append(params.m_buffer)

    thickness.append(film_thickness)
    potentials.append(params.v_high)
    masses.append(params.m_high)


thickness.append(inner_buffer_thickness)
potentials.append(params.v_buffer)
masses.append(params.m_buffer)

thickness.append(central_guide_thickness)
potentials.append(params.v_guide)
masses.append(params.m_guide)

thickness.append(inner_buffer_thickness)
potentials.append(params.v_buffer)
masses.append(params.m_buffer)

for _ in range(5):

    thickness.append(film_thickness)
    potentials.append(params.v_high)
    masses.append(params.m_high)

    thickness.append(buffer_thickness)
    potentials.append(params.v_buffer)
    masses.append(params.m_buffer)

thickness.append(bottom_thickness)
potentials.append(params.v_high)
masses.append(params.m_high)

# do the cumulative sum of the thickness
boundaries = np.cumsum(thickness)

# it gives the location of each layer, for example:
# t1 = 4nm, t2 = 7nm, t3 = 2nm, t4 = 5nm
#boundaries = [ 4, 6, 14, 16, 21]
# so we can know the exact place of the interfaces

# create grid
dz = 0.1 # nm
z_max = boundaries[-1]
z = np.arange(0, z_max, dz)
V_base = np.zeros_like(z)
V_base = V_base.astype(complex)
m_eff = np.zeros_like(z)
delta_V_NV = 0.09 * np.exp(-z / 12.0)

# mapping
# basically where each point z should be placed in the boundaries vector
layers_indices = np.searchsorted(boundaries, z)

eletron_energy_meV = 150



for i in range(len(z)):
    idx = layers_indices[i]
    V_base[i] = potentials[idx]
    m_eff[i] = masses[idx]

## add noise to the system
white_noise = np.random.randn(len(z))
smooth_noise = ndimage.gaussian_filter(white_noise, sigma=5.0)
V_disorder = 1.5 * (smooth_noise / np.std(smooth_noise))

V_base = V_base + V_disorder

# add inelastic scattering rate
Gamma_meV = 0.1
V_base = V_base -1j * (Gamma_meV / 2.0)

V_total_A = V_base + delta_V_NV
V_total_B = V_base - delta_V_NV


# energy sweep
#energies_meV = np.linspace(110, 190, 400)
energies_meV = np.linspace(170.0, 173.5, 400)
transmission_a = []
transmission_b = []

calculate_transmission(V_total_A, m_eff, transmission_a)
calculate_transmission(V_total_B, m_eff, transmission_b)

mu_l = E_F + V_bias/2
mu_r = E_F - V_bias/2

janela_fermi = fermi_dirac(energies_meV, mu_l, Temperature_K)-fermi_dirac(energies_meV, mu_r, Temperature_K)

integrando_a = transmission_a * janela_fermi
integrando_b = transmission_b * janela_fermi

I_A = calculate_current(energies_meV, transmission_a)
I_B = calculate_current(energies_meV, transmission_b)

# alternating the gate voltage
V_gates = np.linspace(-3, 3, 60)
currents_A = []
currents_B = []

for Vg in V_gates:
    print(f"Calculando Gate: {Vg:.2f} meV...")

    # move the potential based on the gate voltage
    v_base_gate = V_base - Vg

    # Add nv center perturbation
    V_total_A = v_base_gate + delta_V_NV
    V_total_B = v_base_gate - delta_V_NV

    trans_a = []
    trans_b = []
    
    calculate_transmission(V_total_A, m_eff, trans_a)
    calculate_transmission(V_total_B, m_eff, trans_b)

    # calculate laundauer current
    I_A = calculate_current(energies_meV, np.array(trans_a))
    I_B = calculate_current(energies_meV, np.array(trans_b))

    currents_A.append(I_A)
    currents_B.append(I_B)

# Current graph
plt.figure(figsize=(8,5))
plt.plot(V_gates, currents_A, color='blue', label='Spin A (NV- extra)')
plt.plot(V_gates, currents_B, color='red', label='Spin B (NV- menor)')
plt.title("Sensor de Spin: corrente vs gate voltage")
plt.xlabel("V_gate (meV)")
plt.ylabel("Corrente (Amperes)")
plt.legend()
plt.show()


######################################################################
#plt.figure(figsize=(10,4))
#plt.plot(energies_meV, integrando_a, color='blue', label='Spin A (NV- extra)')
#plt.plot(energies_meV, integrando_a, color='red', label='Spin B (NV- menor)')
#plt.legend()
#plt.title("Base Potential Profile")
#plt.xlabel(" current ")
#plt.ylabel("Energy (meV)")
#plt.show()
