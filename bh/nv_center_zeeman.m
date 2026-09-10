function nv = nv_center_zeeman(params)
%NV_CENTER_ZEEMAN Axial Zeeman levels of the NV-center spin triplet.
%
% This module describes the NV defect spin, not the propagating electron in
% the semiconductor ARROW. It is meant to track which NV transition is being
% prepared/read by ODMR while the ARROW model computes the resulting change
% in transmission.

if ~isfield(params, "nvZeeman") || ~params.nvZeeman.enable
    nv.enabled = false;
    return;
end

B = params.nvZeeman.B_T(:).';
axisNV = params.nvZeeman.axis(:).';
axisNorm = norm(axisNV);
if axisNorm == 0
    axisNV = [0 0 1];
else
    axisNV = axisNV ./ axisNorm;
end

Bparallel_T = dot(B, axisNV);
Bmag_T = norm(B);
Bperp_T = sqrt(max(Bmag_T.^2 - Bparallel_T.^2, 0));
D_GHz = params.nvZeeman.D_GHz;
gamma_GHz_T = params.nvZeeman.gamma_GHz_T;

ms = [-1 0 1];
freqFromMs0_GHz = D_GHz .* (ms.^2) + gamma_GHz_T .* Bparallel_T .* ms;
freqFromMs0_GHz = freqFromMs0_GHz - freqFromMs0_GHz(ms == 0);

GHz_to_meV = 4.135667696e-3;
energyFromMs0_meV = freqFromMs0_GHz .* GHz_to_meV;

nv.enabled = true;
nv.B_T = B;
nv.Bmag_T = Bmag_T;
nv.axis = axisNV;
nv.Bparallel_T = Bparallel_T;
nv.Bperp_T = Bperp_T;
nv.D_GHz = D_GHz;
nv.gamma_GHz_T = gamma_GHz_T;
nv.ms = ms;
nv.frequencyFromMs0_GHz = freqFromMs0_GHz;
nv.energyFromMs0_meV = energyFromMs0_meV;
nv.f_0_to_minus_GHz = D_GHz - gamma_GHz_T .* Bparallel_T;
nv.f_0_to_plus_GHz = D_GHz + gamma_GHz_T .* Bparallel_T;
nv.E_0_to_minus_meV = nv.f_0_to_minus_GHz .* GHz_to_meV;
nv.E_0_to_plus_meV = nv.f_0_to_plus_GHz .* GHz_to_meV;
nv.stateA_ms = params.nvZeeman.stateA_ms;
nv.stateB_ms = params.nvZeeman.stateB_ms;
nv.stateA_frequency_GHz = frequency_for_ms(nv, nv.stateA_ms);
nv.stateB_frequency_GHz = frequency_for_ms(nv, nv.stateB_ms);
nv.deltaFrequency_AminusB_GHz = nv.stateA_frequency_GHz - nv.stateB_frequency_GHz;
end

function f = frequency_for_ms(nv, stateMs)
[found, idx] = ismember(stateMs, nv.ms);
if found
    f = nv.frequencyFromMs0_GHz(idx);
else
    f = NaN;
end
end
