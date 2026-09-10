%% Parameter sweep for the NV master-equation coupling
%
% This script explores how DeltaV0 depends on physical/calibration
% quantities without rerunning the full transport calculation at every
% point.  Use it to decide plausible ranges for the expensive ARROW runs.

clear; clc; close all;
packageDir = fileparts(mfilename('fullpath'));
addpath(packageDir, "-begin");
cd(packageDir);

base = multipair_electronic_arrow_params("baseline");

density_cm2 = logspace(9, 12, 80);
pulse_ns = linspace(20, 1000, 90);
deltaV_AB_meV = zeros(numel(pulse_ns), numel(density_cm2));
deltaP_AB = zeros(size(deltaV_AB_meV));

for it = 1:numel(pulse_ns)
    for in = 1:numel(density_cm2)
        params = base;
        params.spin.surfaceDensity_cm2 = density_cm2(in);
        params.spin.masterEquation.readoutPulse_s = pulse_ns(it) * 1e-9;
        params.spin.masterEquation.verbose = false;
        [~, master] = update_spin_deltaV_from_master_equation(params);
        deltaV_AB_meV(it, in) = master.deltaV_AminusB_meV;
        deltaP_AB(it, in) = master.deltaP_AminusB;
    end
end

figure("Color", "w", "Name", "Master-equation DeltaV0 sweep");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
imagesc(log10(density_cm2), pulse_ns, deltaV_AB_meV);
axis xy; box on; grid on;
xlabel("log_{10}(n_{NV} / cm^{-2})");
ylabel("readout pulse (ns)");
title("\DeltaV_{0,A-B} from master equation (meV)");
cb = colorbar;
cb.Label.String = "\DeltaV_{A-B} (meV)";
colormap(turbo);

nexttile;
imagesc(log10(density_cm2), pulse_ns, deltaP_AB);
axis xy; box on; grid on;
xlabel("log_{10}(n_{NV} / cm^{-2})");
ylabel("readout pulse (ns)");
title("\DeltaP_{NV^-}(A-B)");
cb = colorbar;
cb.Label.String = "\DeltaP_{NV^-}";
colormap(turbo);

save("Resultados_MasterEquation_CouplingSweep.mat", ...
    "base", "density_cm2", "pulse_ns", "deltaV_AB_meV", "deltaP_AB");
