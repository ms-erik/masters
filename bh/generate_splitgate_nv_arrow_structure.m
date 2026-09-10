function device = generate_splitgate_nv_arrow_structure(params, spinState)
%GENERATE_SPLITGATE_NV_ARROW_STRUCTURE Build the corrected vertical device.
%
% The 1D profile is the centerline below the metal-free optical/NV window.
% The continuous bottom contact is represented by the bottom reservoir. The
% split top gate/contact contributes only an optional smooth electrostatic
% shift, since no top metal crosses the active window.

if nargin < 2
    spinState = "none";
end
spinState = string(spinState);

z = params.vertical.zMin_nm:params.vertical.dz_nm:params.vertical.zMax_nm;
Vstack = splitgate_piecewise_potential(z, params.vertical.boundaries_nm, params.vertical.potential_meV);
meff = splitgate_piecewise_potential(z, params.vertical.boundaries_nm, params.vertical.massEff);
gFactor = splitgate_piecewise_potential(z, params.vertical.boundaries_nm, params.vertical.gFactor);
Vgate = splitgate_effective_gate_shift(z, params);
spin = make_nv_layer_perturbation(z, spinState, params);
Vdisorder = make_disorder_potential(z, params);
Vdecoherence = make_decoherence_potential(z, params);

Vbase = Vstack + Vgate + Vdisorder + Vdecoherence;
Vtotal = Vbase + spin.deltaV_meV;

device.x_nm = z(:).';
device.z_nm = z(:).';
device.Vbase_meV = Vbase(:).';
device.Vstack_meV = Vstack(:).';
device.Vgate_meV = Vgate(:).';
device.Vdisorder_meV = Vdisorder(:).';
device.Vdecoherence_meV = Vdecoherence(:).';
device.V_meV = Vtotal(:).';
device.deltaVspin_meV = spin.deltaV_meV(:).';
device.phaseSpin_rad = spin.phase_rad(:).';
device.spinState = spinState;
device.mEff = meff(:).';
device.gFactor = gFactor(:).';
device.leftPotential_meV = params.leads.leftPotential_meV;
device.rightPotential_meV = params.leads.rightPotential_meV;
device.dx_nm = params.vertical.dz_nm;
device.coordinateName = "z";
device.layout = params.layout;
device.layers = splitgate_grid_to_layers(device.z_nm, device.V_meV, device.mEff, device.gFactor, ...
    device.leftPotential_meV, device.rightPotential_meV);
device.baseLayers = splitgate_grid_to_layers(device.z_nm, device.Vbase_meV, device.mEff, device.gFactor, ...
    device.leftPotential_meV, device.rightPotential_meV);
device.layerBoundaries_nm = params.vertical.boundaries_nm;
device.layerNames = params.vertical.layerNames;
end

function V = splitgate_piecewise_potential(z, boundaries, potentials)
V = potentials(1) * ones(size(z));
for n = 1:numel(boundaries)
    V(z >= boundaries(n)) = potentials(n + 1);
end
V(z < 0) = potentials(1);
end

function Vgate = splitgate_effective_gate_shift(z, params)
Vgate = zeros(size(z));
if ~params.gate.enableSplitGateShift
    return;
end
if isfield(params.gate, "model") && lower(string(params.gate.model)) == "poisson2d"
    Vgate = poisson2d_split_gate_centerline(z, params);
else
    if params.gate.topGateShift_meV == 0
        return;
    end
    tail = exp(-max(z, 0) ./ params.gate.screeningLength_nm) .* (z >= 0);
    Vgate = params.gate.topGateShift_meV .* params.gate.windowAttenuation .* tail;
end
end

function Vdis = make_disorder_potential(z, params)
Vdis = zeros(size(z));
if ~params.disorder.enable || params.disorder.rms_meV == 0
    return;
end
rng(params.disorder.seed);
raw = randn(size(z));
sigmaPts = max(params.disorder.correlationLength_nm / params.vertical.dz_nm, 1);
halfWidth = max(ceil(5 * sigmaPts), 1);
grid = -halfWidth:halfWidth;
kernel = exp(-0.5 * (grid ./ sigmaPts).^2);
kernel = kernel ./ sum(kernel);
smooth = conv(raw, kernel, "same");
smooth = smooth ./ std(smooth(z >= 0));
Vdis(z >= 0) = params.disorder.rms_meV .* smooth(z >= 0);
Vdis(z < 0) = 0;
end

function Vabs = make_decoherence_potential(z, params)
Vabs = zeros(size(z));
if ~params.decoherence.enable || params.decoherence.gamma_meV == 0
    return;
end
mask = z >= 0;
if lower(string(params.decoherence.region)) == "guide"
    names = lower(string(params.vertical.layerNames));
    guideIdx = find(contains(names, "guide"), 1);
    if ~isempty(guideIdx)
        z0 = params.vertical.boundaries_nm(guideIdx);
        z1 = params.vertical.boundaries_nm(guideIdx + 1);
        mask = z >= z0 & z <= z1;
    end
end
Vabs(mask) = -1i * params.decoherence.gamma_meV / 2;
end

function Vgate = poisson2d_split_gate_centerline(z, params)
%POISSON2D_SPLIT_GATE_CENTERLINE Simple Laplace solver for the split pads.
% This is an approximate electrostatic pre-solver. It maps the electrostatic
% potential on the centerline to an electron band-edge shift in meV.
x_um = -params.gate.poisson.xHalfWidth_um:params.gate.poisson.dx_um:params.gate.poisson.xHalfWidth_um;
z_nm_grid = 0:params.gate.poisson.dz_nm:params.vertical.stackThickness_nm;
nx = numel(x_um);
nz = numel(z_nm_grid);
phi = zeros(nz, nx);

halfWindow = params.layout.windowWidth_um / 2;
topMask = abs(x_um) >= halfWindow;
phi(1, topMask) = params.layout.topPadVoltage_V;
phi(1, ~topMask) = 0;
phi(end, :) = params.layout.bottomContactVoltage_V;
phi(:, 1) = params.layout.topPadVoltage_V;
phi(:, end) = params.layout.topPadVoltage_V;

fixed = false(nz, nx);
fixed(1, topMask) = true;
fixed(end, :) = true;
fixed(:, 1) = true;
fixed(:, end) = true;

omega = params.gate.poisson.relaxation;
tol = params.gate.poisson.tolerance_V;
dx2 = (params.gate.poisson.dx_um * 1000)^2;
dz2 = params.gate.poisson.dz_nm^2;
den = 2 / dx2 + 2 / dz2;

for iter = 1:params.gate.poisson.maxIter
    maxDelta = 0;
    % Metal-free top window: approximate open/Neumann boundary.
    phi(1, ~topMask) = phi(2, ~topMask);
    for iz = 2:nz-1
        for ix = 2:nx-1
            if fixed(iz, ix)
                continue;
            end
            candidate = ((phi(iz, ix-1) + phi(iz, ix+1)) / dx2 + ...
                (phi(iz-1, ix) + phi(iz+1, ix)) / dz2) / den;
            newValue = (1 - omega) * phi(iz, ix) + omega * candidate;
            maxDelta = max(maxDelta, abs(newValue - phi(iz, ix)));
            phi(iz, ix) = newValue;
        end
    end
    if maxDelta < tol
        break;
    end
end

[~, ix0] = min(abs(x_um));
phiCenter = interp1(z_nm_grid, phi(:, ix0), max(z, 0), "linear", "extrap");
Vgate = -1000 .* params.gate.leverArm .* phiCenter .* (z >= 0);
end

function layers = splitgate_grid_to_layers(z, V, mEffProfile, gFactorProfile, Vtop, Vbottom)
dz = diff(z);
zc = 0.5 * (z(1:end-1) + z(2:end));
Vc = 0.5 * (V(1:end-1) + V(2:end));
mc = 0.5 * (mEffProfile(1:end-1) + mEffProfile(2:end));
gc = 0.5 * (gFactorProfile(1:end-1) + gFactorProfile(2:end));

layers.xCenter_nm = zc;
layers.width_nm = dz;
layers.V_meV = Vc;
layers.mEff = mc;
layers.gFactor = gc;
layers.Vleft_meV = Vtop;
layers.Vright_meV = Vbottom;
layers.mLeft = mEffProfile(1);
layers.mRight = mEffProfile(end);
layers.gLeft = gFactorProfile(1);
layers.gRight = gFactorProfile(end);
end
