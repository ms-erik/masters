import kwant
import numpy as np
import matplotlib.pyplot as plt
from device_kwant_2d import make_system

def main():
    print("Construindo o sistema Kwant 2D...")
    # Largura de 40nm, Profundidade de 80nm
    syst = make_system(W=40, L=80, a=1.0)
    syst_final = syst.finalized()
    
    # Gerar e salvar a visualização da malha (muito útil para debugar a geometria!)
    print("Gerando gráfico da rede tight-binding...")
    fig = kwant.plot(syst, show=False)
    fig.savefig('kwant_lattice.png', dpi=300, bbox_inches='tight')
    print("Gráfico salvo como 'kwant_lattice.png'.")
    
    # ---------------------------------------------------------
    # Testando o Transporte
    # ---------------------------------------------------------
    # Energia do elétron (em meV)
    energy_meV = 150.0  
    
    print(f"Calculando a S-Matrix para E = {energy_meV} meV...")
    
    # Cenário 1: Centro NV no estado negativo (NV-) e Gate 0V
    params_nv_minus = dict(V_gate=0.0, nv_state='NV-')
    smatrix_minus = kwant.smatrix(syst_final, energy=energy_meV, params=params_nv_minus)
    T_minus = smatrix_minus.transmission(1, 0) # Lead 1 (Drain) <- Lead 0 (Source)
    
    # Cenário 2: Centro NV neutro (NV0) e Gate 0V
    params_nv_zero = dict(V_gate=0.0, nv_state='NV0')
    smatrix_zero = kwant.smatrix(syst_final, energy=energy_meV, params=params_nv_zero)
    T_zero = smatrix_zero.transmission(1, 0)
    
    print("\n=== Resultados da Transmissão (T) ===")
    print(f"Estado NV- (Com Carga) : {T_minus:.5f}")
    print(f"Estado NV0 (Sem Carga) : {T_zero:.5f}")
    
    if T_minus != T_zero:
        print("\nSucesso! O estado do spin modular a transmissão quântica 2D!")
    else:
        print("\nAviso: Não houve modulação. Verifique o perfil de potencial.")

if __name__ == '__main__':
    main()
