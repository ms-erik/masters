function current = landauer_current(energy_meV, Tnone, TA, TB, bias_V, params)
%LANDAUER_CURRENT Current from the Landauer formula.
%
% I = (2e/h) integral T(E) [fL(E)-fR(E)] dE.

c = params.constants;
E_J = energy_meV(:).' * 1e-3 * c.q_C;
EF_J = params.current.fermiEnergy_meV * 1e-3 * c.q_C;
kBT = c.kB_JK * params.current.temperature_K;

muL = EF_J + 0.5 * c.q_C * bias_V;
muR = EF_J - 0.5 * c.q_C * bias_V;

fL = 1 ./ (1 + exp((E_J - muL) ./ kBT));
fR = 1 ./ (1 + exp((E_J - muR) ./ kBT));
window = fL - fR;

h_Js = 2 * pi * c.hbar_Js;
prefactor = 2 * c.q_C / h_Js;

current.none_A = prefactor * trapz(E_J, Tnone(:).' .* window);
current.spinA_A = prefactor * trapz(E_J, TA(:).' .* window);
current.spinB_A = prefactor * trapz(E_J, TB(:).' .* window);
current.deltaI_AminusB_A = current.spinA_A - current.spinB_A;
current.deltaI_AminusNone_A = current.spinA_A - current.none_A;
current.bias_V = bias_V;
end

