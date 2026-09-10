function out = tm_transport(device, energy_meV, params)
%TM_TRANSPORT Coherent 1D transport through a multilayer structure.
%
% The framework is based on transfer-matrix physics. For numerical stability,
% layers are cascaded as equivalent two-port scattering matrices, avoiding
% overflow from large evanescent factors. The equivalent transfer matrix is
% returned for diagnostics.

layers = device.layers;
c = params.constants;

kLeft = electron_k_nm(energy_meV, layers.Vleft_meV, layers.mLeft, c, params);
kRight = electron_k_nm(energy_meV, layers.Vright_meV, layers.mRight, c, params);

S = struct("rL", 0, "tLR", 1, "tRL", 1, "rR", 0);
kPrev = kLeft;
mPrev = layers.mLeft;

for n = 1:numel(layers.width_nm)
    kLayer = electron_k_nm(energy_meV, layers.V_meV(n), layers.mEff(n), c, params);
    Sint = interface_scatter(kPrev, kLayer, mPrev, layers.mEff(n), params);
    Sprop = propagation_scatter(kLayer, layers.width_nm(n));
    S = star_product(S, Sint);
    S = star_product(S, Sprop);
    kPrev = kLayer;
    mPrev = layers.mEff(n);
end

Sout = interface_scatter(kPrev, kRight, mPrev, layers.mRight, params);
S = star_product(S, Sout);

vRatio = real((kRight / layers.mRight) / (kLeft / layers.mLeft));
if abs(real(kLeft)) < params.numeric.minAbsK_nm || vRatio < 0
    vRatio = 0;
end

T = vRatio * abs(S.tLR)^2;
R = abs(S.rL)^2;

if abs(S.tRL) < eps
    M = [Inf Inf; Inf Inf];
else
    M = [S.tLR - S.rL * S.rR / S.tRL, S.rR / S.tRL; ...
        -S.rL / S.tRL, 1 / S.tRL];
end

out.T = min(max(real(T), 0), params.numeric.maxTransmission);
out.R = min(max(real(R), 0), params.numeric.maxTransmission);
out.S = S;
out.M = M;
out.detM = det(M);
out.kLeft_nm = kLeft;
out.kRight_nm = kRight;
end

function S = interface_scatter(k1, k2, m1, m2, params)
% BenDaniel-Duke effective-mass boundary conditions.
q1 = k1 / m1;
q2 = k2 / m2;
den = q1 + q2 + params.numeric.interfaceRegularization;
S.rL = (q1 - q2) / den;
S.tLR = 2 * q1 / den;
S.rR = (q2 - q1) / den;
S.tRL = 2 * q2 / den;
end

function S = propagation_scatter(k, width_nm)
p = exp(1i * k * width_nm);
S.rL = 0;
S.rR = 0;
S.tLR = p;
S.tRL = p;
end

function S = star_product(A, B)
den1 = 1 - B.rL * A.rR;
den2 = 1 - A.rR * B.rL;
S.rL = A.rL + A.tRL * B.rL * A.tLR / den1;
S.tLR = B.tLR * A.tLR / den1;
S.rR = B.rR + B.tLR * A.rR * B.tRL / den2;
S.tRL = A.tRL * B.tRL / den2;
end

