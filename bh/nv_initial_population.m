function p0 = nv_initial_population(ms, nv)
%NV_INITIAL_POPULATION Initial population after optical/MW spin preparation.
%
% The charge is assumed to start mostly in NV-.  The spin preparation is
% encoded by assigning the requested m_s state to G0 or G1 with finite
% preparation purity.  This is intentionally explicit so the user can later
% replace it by measured initialization fidelities.

purity = nv.initialSpinPurity;
chargeMinus = nv.initialChargeNVminus;

p0 = zeros(6, 1);
if ms == 0
    p0(1) = chargeMinus * purity;
    p0(2) = chargeMinus * (1 - purity);
else
    p0(1) = chargeMinus * (1 - purity);
    p0(2) = chargeMinus * purity;
end

% Residual neutral population.
p0(6) = 1 - chargeMinus;
p0 = p0 ./ sum(p0);
end
