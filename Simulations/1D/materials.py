from constants import m0

def get_material_properties(y):

    m_star = m0* 0.08125 * y + 0.041 * m0
    m_eff = m_star / m0

    V = 1062.5 * y

    return m_eff, V




