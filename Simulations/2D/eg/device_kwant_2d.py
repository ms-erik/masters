import kwant
import numpy as np

def make_system(W=40, L=80, a=1.0):
    """
    Constrói a malha de tight-binding 2D do filtro ARR≡W.
    W: Largura total (ao longo de x) em nm.
    L: Comprimento/Profundidade (ao longo de z) em nm.
    a: constante de rede em nm.
    """
    # 1. Definindo a malha quadrada (1 orbital por sítio)
    lat = kwant.lattice.square(a, norbs=1)
    syst = kwant.Builder()
    
    # Energia de hopping (t) baseada na massa efetiva (simplificado InGaAs m* ~ 0.041)
    # t = hbar^2 / (2 * m* * a^2). Em meV para a=1nm, t vale aprox 930 meV.
    t = 930.0 
    
    # Formato geométrico (Região de Espalhamento / Canal)
    def scattering_region(pos):
        x, z = pos
        return (-W/2 <= x <= W/2) and (0 <= z <= L)
        
    # Perfil de Potencial 2D Analítico (Opção 1)
    def potential_2d(x, z, V_gate, nv_state):
        # O gate atua fortemente nas bordas (simulando um split-gate com janela central)
        # Se |x| > W/4, o elétron sente a voltagem do Gate. Na janela central (|x| < W/4), é zero.
        v_g = V_gate if abs(x) > W/4 else 0.0
        
        # O centro NV está no centro da janela na superfície (x=0, z=0)
        q_nv = 1.0 if nv_state == 'NV-' else 0.0
        # Decaimento exponencial 2D (distância radial r)
        r = np.sqrt(x**2 + z**2)
        v_nv = q_nv * 0.09 * np.exp(-r / 12.0)
        
        return v_g + v_nv
        
    # Função "On-site" (Energia que o elétron sente em cada ponto exato da malha)
    def onsite(site, V_gate, nv_state):
        x, z = site.pos
        
        # V_banda = o degrau estrutural da heteroestrutura (os filmes InGaAlAs).
        # Neste primeiro teste, faremos zero para validar o transporte 2D.
        V_banda = 0.0
        V_eletr = potential_2d(x, z, V_gate, nv_state)
        
        # Fundo da banda de condução p/ malha quadrada: 4*t + V(x,z)
        return 4 * t + V_banda + V_eletr
        
    # Preenchendo a região de espalhamento com os sítios e conectando (hopping -t)
    syst[lat.shape(scattering_region, (0, L/2))] = onsite
    syst[lat.neighbors()] = -t
    
    # 2. Construindo os Contatos (Leads Semi-infinitos)
    
    # Lead Topo (Source) - Vem de z = -infinito em direção a z=0
    sym_top = kwant.TranslationalSymmetry((0, -a))
    lead_top = kwant.Builder(sym_top)
    
    def lead_shape_top(pos):
        x, z = pos
        return (-W/2 <= x <= W/2)
        
    # O contato de entrada é metálico/altamente dopado, V = 0
    lead_top[lat.shape(lead_shape_top, (0, 0))] = 4 * t
    lead_top[lat.neighbors()] = -t
    
    # Lead Base (Drain) - Sai de z = L e vai para z = infinito
    sym_bottom = kwant.TranslationalSymmetry((0, a))
    lead_bottom = kwant.Builder(sym_bottom)
    
    def lead_shape_bottom(pos):
        x, z = pos
        return (-W/2 <= x <= W/2)
        
    lead_bottom[lat.shape(lead_shape_bottom, (0, 0))] = 4 * t
    lead_bottom[lat.neighbors()] = -t
    
    # Anexando (Lead 0 é Topo/Source, Lead 1 é Base/Drain)
    syst.attach_lead(lead_top)
    syst.attach_lead(lead_bottom)
    
    return syst
