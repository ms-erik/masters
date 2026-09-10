function [K, labels] = nv_charge_master_matrix(nv)
%NV_CHARGE_MASTER_MATRIX Rate matrix for a compact NV-/NV0 model.
%
% State vector convention:
%   p = [G0; G1; E0; E1; S; C0]
%
% G0 : NV- triplet ground population in m_s = 0
% G1 : NV- triplet ground population in m_s = +-1 lumped manifold
% E0 : NV- excited triplet population reached from G0
% E1 : NV- excited triplet population reached from G1
% S  : effective singlet/intersystem-crossing manifold
% C0 : neutral charge state NV0
%
% The matrix K is column-stochastic in the master-equation sense:
%   dp/dt = K p,
% and each column sums to zero.

r = nv.rates;
labels = ["G0","G1","E0","E1","S","NV0"];
K = zeros(6, 6);

% Optical excitation within NV-.
K = add_rate(K, 3, 1, r.k_ex_532_s);      % G0 -> E0
K = add_rate(K, 4, 2, r.k_ex_532_s);      % G1 -> E1

% Radiative decay.
K = add_rate(K, 1, 3, r.k_rad_s);         % E0 -> G0
K = add_rate(K, 2, 4, r.k_rad_s);         % E1 -> G1

% Spin-dependent intersystem crossing.
K = add_rate(K, 5, 3, r.k_isc_0_s);       % E0 -> S
K = add_rate(K, 5, 4, r.k_isc_1_s);       % E1 -> S

% Singlet decay polarizes into m_s=0 preferentially.
K = add_rate(K, 1, 5, r.k_s_to_0_s);      % S -> G0
K = add_rate(K, 2, 5, r.k_s_to_1_s);      % S -> G1

% Photoionization to NV0.
K = add_rate(K, 6, 3, r.k_ion_0_s);       % E0 -> NV0
K = add_rate(K, 6, 4, r.k_ion_1_s);       % E1 -> NV0
K = add_rate(K, 6, 5, r.k_ion_singlet_s); % S  -> NV0

% Recombination from NV0 into NV-.
K = add_rate(K, 1, 6, r.k_rec_to_0_s);    % NV0 -> G0
K = add_rate(K, 2, 6, r.k_rec_to_1_s);    % NV0 -> G1

% Optional microwave-driven mixing during the readout window.
if isfield(r, "k_mw_mix_s") && r.k_mw_mix_s > 0
    K = add_rate(K, 2, 1, r.k_mw_mix_s);  % G0 <-> G1
    K = add_rate(K, 1, 2, r.k_mw_mix_s);
end
end

function K = add_rate(K, to, from, rate)
if rate == 0
    return;
end
K(to, from) = K(to, from) + rate;
K(from, from) = K(from, from) - rate;
end
