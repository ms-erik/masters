function spin = make_nv_layer_perturbation(z_nm, spinState, params)
%MAKE_NV_LAYER_PERTURBATION Effective perturbation from a 2D NV sheet.
%
% The NV layer is on the lower face of a removable diamond slab. Since the
% diamond is not the transport channel, the perturbation is represented as a
% two-stage effective potential:
%   1) geometric/dielectric coupling from the NV sheet to the interface;
%   2) screening only inside the semiconductor stack.
%
% Model:
%   deltaV_NV(z,ms) = C_interface * deltaV_ms *
%       (d0/(d_NV + d0))^n * exp(-z/lambda_scr), z >= 0.
%
% This is still phenomenological, but it separates the diamond/interface
% coupling from semiconductor screening.

spinState = lower(string(spinState));

switch spinState
    case "a"
        amplitude_meV = params.spin.deltaV_A_meV;
    case "b"
        amplitude_meV = params.spin.deltaV_B_meV;
    otherwise
        amplitude_meV = 0;
end

d = params.spin.distanceFromStack_nm;
lambda = params.spin.screeningLength_nm;
d0 = params.spin.distanceRegularization_nm;
nPower = params.spin.distancePower;
sheetSigma = max(params.spin.sheetThickness_nm / 2.355, 0.2);

% z=0 is the semiconductor surface. The NV sheet is above it at z=-d.
zNV = -d;
penetration = exp(-max(z_nm, 0) ./ lambda) .* (z_nm >= 0);

% The NV sheet itself is not part of the electronic current path. Coupling
% to the semiconductor is split into an interface transmission coefficient
% and a distance law in the diamond/interface region.
interfaceWeight = local_interface_coupling(params);
distanceWeight = (d0 ./ (max(d, 0) + d0)).^nPower;
stackTail = interfaceWeight .* distanceWeight .* penetration;
profile = stackTail;
profile = params.spin.lateralFillFactor .* profile;

spin.deltaV_meV = amplitude_meV .* profile;
spin.phase_rad = params.spin.phaseAmplitude_rad .* profile .* sign(amplitude_meV);
spin.nvSheetPosition_nm = zNV;
spin.nvSheetSigma_nm = sheetSigma;
spin.interfaceCoupling = interfaceWeight;
spin.distanceWeight = distanceWeight;
spin.screeningLength_nm = lambda;
spin.distancePower = nPower;
end

function cInterface = local_interface_coupling(params)
mode = "manual";
if isfield(params.spin, "interfaceCouplingMode")
    mode = lower(string(params.spin.interfaceCouplingMode));
end

switch mode
    case "permittivity"
        epsDiamond = params.spin.epsilonDiamond;
        epsSemi = params.spin.epsilonSemiconductor;
        cInterface = 2 .* epsDiamond ./ (epsDiamond + epsSemi);
    otherwise
        cInterface = params.spin.interfaceCoupling;
end

cInterface = max(cInterface, 0);
end
