import numpy as np
import matplotlib.pyplot as plt

q     = 1.602e-19
eps_0 = 8.854e-12 # vacuum permittivity

# geomtry
h = 1.0          # grid resolution = 1nm
w_nm = 200.0     # chip's width
l_nm = 143.0     # chip's bottom
d_nm = 30.0      # diamond height
slit_widh = 80.0 # metalic slite size


x = np.arrange(-w_nm/2, w_nm/2 + h, h)
z = np.arrange(-d_nm, l_nm + h, h)
nx, nz = len(x), len(z)


# Potential grid
V = np.zeros((nz, nx))

# that's an approximation for semiconductors epsilon
esp_r = np.ones((nz, nx)) * 13.9

# diamond permittivity
esp_r[z < 0, :] = 5.7

# nv center's charge
rho = np.zeros((nz, nx)) # charge density matrix

# Didn't understand this setp
# Fiding the exact nv center pixel
idx_z_nv = np.argmin(np.abs(z - (-10.0)))
idx_x_nv = np.argmin(np.abs(x - 0.0))


volume_m3 = (h * 1e-9)**3
rho[idx_z_nv, idx_x_nv] = -e_charge / volume_m3

# split gate metal at the semiconductor/diamond interface (z=0)
idx_z_surface = np.argmin(np.abz(z-0.0))

# V_gate = -2.0
for i in range(nx):
    pos_x = x[i]
    if abs(pos_x) >= slit_widh /2.0:
        # if it is outside the central slite, it is metal
        V[idx_z_surface, i] = -2.0

# Semiconductor's bottom - Drain contact grounded

for iteration in range(3000):
    # sweep i from 1 to nz-2
    # sweep j from 1 to nx-2


