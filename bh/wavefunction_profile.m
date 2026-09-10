function wf = wavefunction_profile(device, energy_meV, params)
%WAVEFUNCTION_PROFILE Reconstruct incident + reflected + transmitted wave.

layers = device.layers;
c = params.constants;
x = device.x_nm;

transport = tm_transport(device, energy_meV, params);
r = transport.S.rL;

kLeft = electron_k_nm(energy_meV, layers.Vleft_meV, layers.mLeft, c, params);
A = 1;
B = r;
xBoundary = x(1);

psi = zeros(size(x));
kPrev = kLeft;
mPrev = layers.mLeft;

for n = 1:numel(layers.width_nm)
    kLayer = electron_k_nm(energy_meV, layers.V_meV(n), layers.mEff(n), c, params);
    [C, D] = match_interface(A, B, kPrev, kLayer, mPrev, layers.mEff(n));
    xStart = x(n);
    xEnd = x(n + 1);
    mask = x >= xStart & x <= xEnd;
    xi = x(mask) - xStart;
    psi(mask) = C .* exp(1i * kLayer .* xi) + D .* exp(-1i * kLayer .* xi);
    A = C .* exp(1i * kLayer .* layers.width_nm(n));
    B = D .* exp(-1i * kLayer .* layers.width_nm(n));
    kPrev = kLayer;
    mPrev = layers.mEff(n);
    xBoundary = xEnd;
end

kRight = electron_k_nm(energy_meV, layers.Vright_meV, layers.mRight, c, params);
[C, D] = match_interface(A, B, kPrev, kRight, mPrev, layers.mRight);
mask = x >= xBoundary;
psi(mask) = C .* exp(1i * kRight .* (x(mask) - xBoundary)) + ...
    D .* exp(-1i * kRight .* (x(mask) - xBoundary));

wf.energy_meV = energy_meV;
wf.x_nm = x;
wf.psi = psi;
wf.probabilityDensity = abs(psi).^2;
wf.transmission = transport.T;
wf.reflection = transport.R;
end

function [C, D] = match_interface(A, B, k1, k2, m1, m2)
q1 = k1 / m1;
q2 = k2 / m2;
sumPsi = A + B;
fluxPsi = q1 * (A - B) / q2;
C = 0.5 * (sumPsi + fluxPsi);
D = 0.5 * (sumPsi - fluxPsi);
end

