from dataclasses import dataclass
from constants import hbar, m0, EV_TO_JOULE, nm_to_m, m_to_nm
import numpy as np
from materials import get_material_properties

# Constantes físicas iniciais
targetEnergy_meV = 150

# Potencial e massas do InGaAs
# (High-k) filmes refletores

#Potencial e massa do InGaAlAs
# low-k - buffer central

# k = sqrt(2m* (E-V))/ hbar

class DeviceParams:
    def __init__(self, targetEnergy_meV=150.0):
        self.targetEnergy_meV = targetEnergy_meV

        # aluminium fraction (y) per region 
        self.y_high = 0.0
        self.y_buffer = 0.08
        self.y_guide = 0.10

        self.m_high, self.v_high     = get_material_properties(self.y_high) 
        self.m_buffer, self.v_buffer = get_material_properties(self.y_buffer)
        self.m_guide, self.v_guide   = get_material_properties(self.y_guide)

        self.d_high_nm = 0.0
        self.d_guide_nm = 0.0

    def calculate_k(self, E_meV, V_meV, m_eff):
        # convert to J
        e_j = E_meV * EV_TO_JOULE * 1e-3
        v_j = V_meV * EV_TO_JOULE * 1e-3
        m_star = m_eff * m0

        k = np.sqrt(2 * m_star * complex(e_j - v_j)) / hbar
        return k

    def compute_thickness(self):
        k_high = self.calculate_k(self.targetEnergy_meV, self.v_high, self.m_high)
        k_guide = self.calculate_k(self.targetEnergy_meV, self.v_guide, self.m_guide)

        d_high = np.real((0.5 * np.pi) / k_high)
        d_guide = np.real((2 * np.pi) / k_guide)

        self.d_high_nm = d_high * m_to_nm
        self.d_guide_nm = d_guide * m_to_nm
