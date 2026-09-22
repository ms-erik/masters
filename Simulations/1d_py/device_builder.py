import numpy as np

def build_arrw_structure(num_reflectors_per_side=3):
    """
    Builds a mock ARR≡W structure for initial testing.
    Returns a list of layer dictionaries.
    """
    contact = {'V': 0.0, 'm': 0.05, 'd': 0}
    
    layer_high_k = {'V': 0.0, 'm': 0.05, 'd': 5.0} # 5 nm
    layer_low_k  = {'V': 0.2, 'm': 0.07, 'd': 3.0} # 3 nm
    
    central_guide = {'V': 0.2, 'm': 0.07, 'd': 15.0} # 15 nm
    
    layers = []
    
    layers.append(contact)
    for _ in range(num_reflectors_per_side):
        layers.append(layer_high_k)
        layers.append(layer_low_k)
        
    layers.append(central_guide)
    
    for _ in range(num_reflectors_per_side):
        layers.append(layer_low_k)
        layers.append(layer_high_k)
        
    layers.append(contact)
    
    return layers

def build_single_barrier():
    """
    Builds a single barrier for testing evanescent decay.
    """
    contact = {'V': 0.0, 'm': 0.05, 'd': 0}
    barrier = {'V': 0.3, 'm': 0.07, 'd': 5.0}
    
    return [contact, barrier, contact]
