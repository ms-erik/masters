import numpy as np
from constants import m0
from device_builder import DeviceParams

params = DeviceParams(targetEnergy_meV=150.0)                                                 
params.compute_thickness()                                                                    
                                                                                                
thickness = [                                                                                 
    3.0,                                                                                      
    5.5, params.d_high_nm, 5.5, params.d_high_nm, 5.5, params.d_high_nm, 5.5, params.d_high_nm,
5.5, params.d_high_nm,                                                                          
    6.0,                                                                                      
    params.d_guide_nm,                                                                        
    6.0,                                                                                      
    params.d_high_nm, 5.5, params.d_high_nm, 5.5, params.d_high_nm, 5.5, params.d_high_nm, 5.5,
params.d_high_nm, 5.5,                                                                          
    8.0                                                                                       
]                                                                                             
boundaries = np.cumsum(thickness)                                                             
                                                                                                
potentials = [                                                                                
    params.v_high,                                                                            
    params.v_buffer, params.v_high, params.v_buffer, params.v_high, params.v_buffer, params.  
v_high, params.v_buffer, params.v_high, params.v_buffer, params.v_high,                         
    params.v_buffer, params.v_guide, params.v_buffer,                                         
    params.v_high, params.v_buffer, params.v_high, params.v_buffer, params.v_high, params.    
v_buffer, params.v_high, params.v_buffer, params.v_high, params.v_buffer,                       
    params.v_high                                                                             
]                                                                                             
                                                                                                
masses = [                                                                                    
    params.m_high,                                                                            
    params.m_buffer, params.m_high, params.m_buffer, params.m_high, params.m_buffer, params.  
m_high, params.m_buffer, params.m_high, params.m_buffer, params.m_high,                         
    params.m_buffer, params.m_guide, params.m_buffer,                                         
    params.m_high, params.m_buffer, params.m_high, params.m_buffer, params.m_high, params.    
m_buffer, params.m_high, params.m_buffer, params.m_high, params.m_buffer,                       
    params.m_high
]
    	

def get_material(z):
    if z < 0:
        return potentials[0], masses[0]

    idx = np.searchsorted(boundaries, z)

    if idx >= len(potentials):
        idx = -1

    return potentials[idx], masses[idx]
