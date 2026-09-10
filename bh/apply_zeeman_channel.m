function channelDevice = apply_zeeman_channel(device, params, spinProjection)
%APPLY_ZEEMAN_CHANNEL Add diagonal electron Zeeman shift to one spin channel.
%
% spinProjection = +1 for the sigma_B = +1 channel and -1 for sigma_B = -1.
% The shift is DeltaE_Z(z) = spinProjection * 0.5 * g(z) * mu_B * |B|.
% This is the exact reduction of a diagonal 2x2 Zeeman Hamiltonian when
% spin-flip terms are neglected.

channelDevice = device;

if ~isfield(params, "zeeman") || ~params.zeeman.enable
    return;
end

Bmag_T = norm(params.zeeman.B_T);
if Bmag_T == 0
    return;
end

muB_meV_T = 5.7883818060e-2;
shiftGrid_meV = spinProjection .* 0.5 .* device.gFactor .* muB_meV_T .* Bmag_T;

channelDevice.V_meV = device.V_meV + shiftGrid_meV;
channelDevice.Vbase_meV = device.Vbase_meV + shiftGrid_meV;

layerShift_meV = spinProjection .* 0.5 .* device.layers.gFactor .* muB_meV_T .* Bmag_T;
channelDevice.layers.V_meV = device.layers.V_meV + layerShift_meV;
channelDevice.layers.Vleft_meV = device.layers.Vleft_meV + ...
    spinProjection .* 0.5 .* device.layers.gLeft .* muB_meV_T .* Bmag_T;
channelDevice.layers.Vright_meV = device.layers.Vright_meV + ...
    spinProjection .* 0.5 .* device.layers.gRight .* muB_meV_T .* Bmag_T;

if isfield(device, "baseLayers")
    channelDevice.baseLayers.V_meV = device.baseLayers.V_meV + layerShift_meV;
    channelDevice.baseLayers.Vleft_meV = channelDevice.layers.Vleft_meV;
    channelDevice.baseLayers.Vright_meV = channelDevice.layers.Vright_meV;
end

channelDevice.zeemanShift_meV = shiftGrid_meV;
channelDevice.zeemanSpinProjection = spinProjection;
end
