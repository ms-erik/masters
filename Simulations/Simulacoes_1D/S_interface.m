function S = S_interface(k_L, m_L, k_R, m_R)
% S_INTERFACE: Calcula a Matriz de Espalhamento (S-Matrix) para uma
% interface simples entre uma camada Esquerda (L) e uma Direita (R).
%
% Entradas:
%   k_L, k_R : número de onda nas camadas Esquerda e Direita.
%   m_L, m_R : massas efetivas nas camadas Esquerda e Direita.
%
% Saída:
%   S : Matriz 2x2 de espalhamento [r_L, t_RL; t_LR, r_R]

    % A conservação de corrente usa as "velocidades" v = k / m*
    v_L = k_L / m_L;
    v_R = k_R / m_R;
    
    % Coeficientes de Reflexão (r) e Transmissão (t)
    r = (v_L - v_R) / (v_L + v_R);
    t_LR = 2 * v_L / (v_L + v_R);
    t_RL = 2 * v_R / (v_L + v_R);
    
    % Construção da S-matrix
    S = zeros(2,2);
    S(1,1) = r;          % Reflexão de volta para L (S11)
    S(2,2) = -r;         % Reflexão de volta para R (S22)
    S(2,1) = t_LR;       % Transmissão de L para R (S21)
    S(1,2) = t_RL;       % Transmissão de R para L (S12)
end
