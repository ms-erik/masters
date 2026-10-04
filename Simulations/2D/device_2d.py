import kwant                                                                                  
import matplotlib.pyplot as plt                                                               
import numpy as np                                                                            
from constants import HBAR2_M0
from materials import get_material
                                                                                                
# lattice parameters
a = 1.0 # set grid to 1 nm
W_nm = 50.0 # channel width
L_nm = 143.0 # channel length

# band energy based on the effective mass
# t = haber^2 / (2 m* a²)

def make_system():                                                                            
    # squared lattice
    lat = kwant.lattice.square(a)                                                             
    sys = kwant.Builder()                                                                     

    def onsite(site, V_gate, nv_state):
        x, z = site.pos
        V_z, m_z = get_material(z)

        # lateral confinment 
        k_confine = 0.5
        V_lateral = k_confine * (x**2)

        # local cinetic energy
        t = HBAR2_M0 / (2.0 * m_z * (a**2))
        # 2d  -> sum 47 to the base potential
        return V_z + 4.0 * t
    
    def hopping_z(site1, site2, V_gate, nv_state):
        # vertical hopping
        z1 = site1.pos[1]
        z2 = site2.pos[1]

        _, m_z1 = get_material(z1)
        _, m_z2 = get_material(z2)

        # in the tight band formalism, the hopping energy 
        # between adjacent sites is inversely proportional 
        # to the effective mass. In this case ( by conservation )
        # it is the spatial avarage between them
        m_avg = (m_z1 + m_z2) / 2.0
        # now the band energy based on the effective mass
        t = HBAR2_M0 / (2.0 * m_avg * (a**2))
        return -t

    def hopping_x(site1, site2, V_gate, nv_state):
        # horizontal jump - does not change with height z
        z = site1.pos[1]
        _, m = get_material(z)

        # the energy jump uses the current mass, no stress
        t = HBAR2_M0 / (2.0 * m * (a**2))
        return -t
                                                                                                    
                                                                                                
    # define geometry
    def shape(pos):                                                                           
        x, z = pos                                                                            
        return (-W_nm/2 <= x <= W_nm/2) and (0 <= z <= L_nm)                                  
                                                                                                
    # population the grid
    sys[lat.shape(shape, (0, 0))] = onsite
    sys[kwant.builder.HoppingKind((0, 1), lat, lat)] = hopping_z
    sys[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopping_x

    # Building the top lead - Source
    # this -a says which side the lead points to infinity
    sym_top = kwant.TranslationalSymmetry((0, -a))
    lead_top = kwant.Builder(sym_top)
                                                                                                
    def lead_shape(pos):                                                                      
        x, z = pos                                                                            
        return (-W_nm/2 <= x <= W_nm/2)                                                       
                                                                                                
    lead_top[lat.shape(lead_shape, (0, 0))] = onsite
    lead_top[kwant.builder.HoppingKind((0,1), lat, lat)] = hopping_z
    lead_top[kwant.builder.HoppingKind((1,0), lat, lat)] = hopping_x
                                                                                                
    # constucts the inferior lead
    sym_bot = kwant.TranslationalSymmetry((0, a))                                             
    lead_bot = kwant.Builder(sym_bot)                                                         
    lead_bot[lat.shape(lead_shape, (0, 0))] = onsite
    lead_bot[kwant.builder.HoppingKind((0,1), lat, lat)] = hopping_z
    lead_bot[kwant.builder.HoppingKind((1,0), lat,lat)] = hopping_x

    # attach the terminals to the central device
    sys.attach_lead(lead_top)                                                                 
    sys.attach_lead(lead_bot)                                                                 
                                                                                                
    return sys                                                                                
                                                                                                
dispositivo = make_system()
dispositivo_final = dispositivo.finalized()
                                                                                                
kwant.plot(dispositivo)
plt.show()
