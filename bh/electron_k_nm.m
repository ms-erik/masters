function k_nm = electron_k_nm(E_meV, V_meV, mEff, constants, params)
%ELECTRON_K_NM Complex electron wavevector in nm^-1.
%
% V_meV may be complex. A small negative imaginary part can be used as a
% phenomenological absorbing/decoherence potential.

dE_J = (E_meV - V_meV) * 1e-3 * constants.q_C;
m_kg = mEff * constants.m0_kg;
k_m = sqrt(2 * m_kg .* dE_J) / constants.hbar_Js;
k_nm = k_m * 1e-9;

small = abs(k_nm) < params.numeric.minAbsK_nm;
k_nm(small) = params.numeric.minAbsK_nm + 0i;
end
