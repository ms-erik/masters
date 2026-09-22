function S = S_propagation(k, d)
% S_PROPAGATION: Calcula a Matriz de Espalhamento para a propagação
% (ou decaimento evanescente) de um elétron através de uma camada.
%
% Entradas:
%   k : número de onda (real para propagante, imaginário para evanescente)
%   d : espessura da camada
%
% Saída:
%   S : Matriz 2x2 de espalhamento para a propagação.

    % O fator de fase que acumula na viagem pela camada
    % Se k for imaginário, exp(1i * k * d) se torna exp(-|k|*d), um decaimento real.
    phase_factor = exp(1i * k * d);
    
    % Construção da S-matrix
    % Como o meio é homogêneo, não há reflexão (S11 = 0, S22 = 0)
    S = zeros(2,2);
    S(1,1) = 0;
    S(2,2) = 0;
    S(2,1) = phase_factor; % Transmissão no sentido positivo
    S(1,2) = phase_factor; % Transmissão no sentido negativo
end
