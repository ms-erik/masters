function S_total = S_cascade(S_A, S_B)
% S_CASCADE: Combina (cascateia) duas Matrizes de Espalhamento S_A e S_B
% usando o Produto Estrela de Redheffer. É numericamente estável.
%
% S_A conecta a região 1 à região 2.
% S_B conecta a região 2 à região 3.
% S_total conectará a região 1 diretamente à região 3.

    % Extraindo as componentes das matrizes de entrada
    r11_A = S_A(1,1); t22_A = S_A(1,2);
    t11_A = S_A(2,1); r22_A = S_A(2,2);
    
    r11_B = S_B(1,1); t22_B = S_B(1,2);
    t11_B = S_B(2,1); r22_B = S_B(2,2);
    
    % O denominador captura as múltiplas reflexões entre as matrizes A e B
    % (Fabry-Perot interno)
    den = 1 - r22_A * r11_B;
    
    % Construindo a matriz resultante
    S_total = zeros(2,2);
    
    % Reflexão Total Esquerda -> Esquerda
    S_total(1,1) = r11_A + (t22_A * r11_B * t11_A) / den;
    
    % Reflexão Total Direita -> Direita
    S_total(2,2) = r22_B + (t11_B * r22_A * t22_B) / den;
    
    % Transmissão Total Esquerda -> Direita (Crucial para a corrente)
    S_total(2,1) = (t11_B * t11_A) / den;
    
    % Transmissão Total Direita -> Esquerda
    S_total(1,2) = (t22_A * t22_B) / den;
end
