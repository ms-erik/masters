function plot_splitgate_nv_arrow_results(devices, results, wf, params)
%PLOT_SPLITGATE_NV_ARROW_RESULTS Plots for split-top/continuous-bottom device.

labelA = "NV state A";
labelB = "NV state B";
if isfield(params, "spin")
    if isfield(params.spin, "stateA_label")
        labelA = params.spin.stateA_label;
    end
    if isfield(params.spin, "stateB_label")
        labelB = params.spin.stateB_label;
    end
end

figure("Color", "w", "Name", "Split-top / continuous-bottom NV-ARROW geometry");
tiledlayout(4, 1, "TileSpacing", "compact");

nexttile;
plot(devices.none.z_nm, devices.none.Vstack_meV, "Color", [0.55 0.55 0.55], "LineWidth", 1.1); hold on;
plot(devices.none.z_nm, devices.none.Vbase_meV, "k-", "LineWidth", 1.4);
plot(devices.spinA.z_nm, devices.spinA.V_meV, "r-", "LineWidth", 1.0);
plot(devices.spinB.z_nm, devices.spinB.V_meV, "b-", "LineWidth", 1.0);
xline(0, "Color", [0.2 0.6 0.2], "LineStyle", "--", "LineWidth", 1.0);
ylabel("V(z) (meV)");
title("Centerline potential below metal-free NV window");
legend("ARROW stack", "no NV + gate", labelA, labelB, "surface", "Location", "best");
grid on; box on;

nexttile;
plot(devices.none.z_nm, devices.none.Vgate_meV, "Color", [0.65 0.35 0.05], "LineWidth", 1.2); hold on;
plot(devices.spinA.z_nm, devices.spinA.deltaVspin_meV, "r-", "LineWidth", 1.3);
plot(devices.spinB.z_nm, devices.spinB.deltaVspin_meV, "b-", "LineWidth", 1.3);
xline(-params.spin.distanceFromStack_nm, "m--", "LineWidth", 1.0);
xline(0, "Color", [0.2 0.6 0.2], "LineStyle", "--", "LineWidth", 1.0);
ylabel("shift (meV)");
title("Effective split-gate and 2D-NV perturbations");
legend("split top gate effective shift", labelA, labelB, "2D NV layer", "surface", "Location", "best");
grid on; box on;

nexttile;
draw_splitgate_layer_bands(params);
xlabel("z (nm)");
ylabel("layer");
title("Deposited horizontal ARROW stack");
grid on; box on;

nexttile;
draw_contact_cartoon(params);
axis off;
title("Contact concept: split top gate/contact, continuous bottom contact");

figure("Color", "w", "Name", "Split-gate NV-ARROW transmission");
tiledlayout(3, 1, "TileSpacing", "compact");

nexttile;
plot(results.energy_meV, results.none.T, "k-", "LineWidth", 1.1); hold on;
plot(results.energy_meV, results.spinA.T, "r-", "LineWidth", 1.1);
plot(results.energy_meV, results.spinB.T, "b-", "LineWidth", 1.1);
ylabel("T(E)");
title("Vertical coherent tunneling through split-top/continuous-bottom device");
legend("no NV perturbation", labelA, labelB, "Location", "best");
grid on; box on;
ylim([0, min(1.1, max([results.none.T results.spinA.T results.spinB.T]) * 1.05 + 0.02)]);

nexttile;
plot(results.energy_meV, results.deltaT_AminusB, "m-", "LineWidth", 1.2);
ylabel("\DeltaT_{A-B}");
title("Spin-sensitive differential transmission");
grid on; box on;

nexttile;
plot(results.energy_meV, results.spinA.R, "r-", "LineWidth", 1.0); hold on;
plot(results.energy_meV, results.spinB.R, "b-", "LineWidth", 1.0);
xlabel("Energy (meV)");
ylabel("R(E)");
title("Reflection");
legend(labelA, labelB, "Location", "best");
grid on; box on;

figure("Color", "w", "Name", "Split-gate wavefunction at resonance");
maxProb = max([wf.none.probabilityDensity wf.spinA.probabilityDensity wf.spinB.probabilityDensity]);
if maxProb == 0
    maxProb = 1;
end
scale = max(devices.none.Vbase_meV) / maxProb;
plot(devices.none.z_nm, devices.none.Vbase_meV, "Color", [0.4 0.4 0.4], "LineWidth", 1.2); hold on;
plot(wf.none.x_nm, wf.none.probabilityDensity * scale, "k-", "LineWidth", 1.0);
plot(wf.spinA.x_nm, wf.spinA.probabilityDensity * scale, "r-", "LineWidth", 1.0);
plot(wf.spinB.x_nm, wf.spinB.probabilityDensity * scale, "b-", "LineWidth", 1.0);
xline(0, "Color", [0.2 0.6 0.2], "LineStyle", "--", "LineWidth", 1.0);
xlabel("z (nm)");
ylabel("V (meV), scaled |\psi|^2");
title(sprintf("Centerline probability density near E = %.6g meV", wf.spinA.energy_meV));
legend("V(z)", "no NV", labelA, labelB, "surface", "Location", "best");
grid on; box on;

figure("Color", "w", "Name", "Split-gate current response");
bar(categorical(["none", labelA, labelB, "\DeltaI A-B"]), ...
    [results.current.none_A, results.current.spinA_A, ...
    results.current.spinB_A, results.current.deltaI_AminusB_A]);
ylabel("Current (A)");
title(sprintf("Landauer current at V_b = %.3g V, E_F = %.4g meV", ...
    results.current.bias_V, params.current.fermiEnergy_meV));
grid on; box on;

if ~isempty(results.resonances.energy_meV)
    figure("Color", "w", "Name", "Split-gate resonance map");
    scatter(results.resonances.energy_meV, results.resonances.Q, ...
        42, results.resonances.prominence, "filled");
    xlabel("Resonance energy (meV)");
    ylabel("Estimated Q = E/FWHM");
    title("Detected resonances for " + labelA);
    cb = colorbar;
    cb.Label.String = "Prominence";
    grid on; box on;
end

if isfield(results, "nvZeeman") && results.nvZeeman.enabled
    plot_nv_zeeman_levels(results.nvZeeman, labelA, labelB);
end

if isfield(results, "zeeman") && results.zeeman.enabled
    figure("Color", "w", "Name", "Zeeman spin-resolved transmission");
    tiledlayout(2, 1, "TileSpacing", "compact");
    nexttile;
    plot(results.energy_meV, results.spinA.zeeman.up.T, "r-", "LineWidth", 1.0); hold on;
    plot(results.energy_meV, results.spinA.zeeman.down.T, "r--", "LineWidth", 1.0);
    plot(results.energy_meV, results.spinB.zeeman.up.T, "b-", "LineWidth", 1.0);
    plot(results.energy_meV, results.spinB.zeeman.down.T, "b--", "LineWidth", 1.0);
    ylabel("T(E)");
    title(sprintf("ARROW electron Zeeman channels, |B| = %.3g T", results.zeeman.Bmag_T));
    legend(labelA + " electron up", labelA + " electron down", ...
        labelB + " electron up", labelB + " electron down", "Location", "best");
    grid on; box on;

    nexttile;
    plot(results.energy_meV, results.zeeman.deltaT_up_AminusB, "Color", [0.7 0 0.7], "LineWidth", 1.0); hold on;
    plot(results.energy_meV, results.zeeman.deltaT_down_AminusB, "Color", [0.2 0.2 0.8], "LineWidth", 1.0);
    xlabel("Energy (meV)");
    ylabel("\DeltaT_{A-B}");
    title("ARROW electron spin-channel differential transmission");
    legend("electron up channel", "electron down channel", "Location", "best");
    grid on; box on;
end
end

function plot_nv_zeeman_levels(nv, labelA, labelB)
figure("Color", "w", "Name", "NV-center Zeeman levels");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
bar(categorical(["m_s=-1", "m_s=0", "m_s=+1"]), nv.frequencyFromMs0_GHz);
ylabel("frequency from m_s=0 (GHz)");
title(sprintf("NV-center axial Zeeman levels, B_parallel = %.4g T", nv.Bparallel_T));
grid on; box on;

nexttile;
bar(categorical(["0 -> -1", "0 -> +1"]), [nv.f_0_to_minus_GHz nv.f_0_to_plus_GHz]);
ylabel("ODMR frequency (GHz)");
title(sprintf("%s: m_s=%d, %s: m_s=%d", labelA, nv.stateA_ms, labelB, nv.stateB_ms));
grid on; box on;
end

function draw_splitgate_layer_bands(params)
boundaries = params.vertical.boundaries_nm;
names = params.vertical.layerNames;
for n = 1:numel(names)
    z0 = boundaries(n);
    z1 = boundaries(n + 1);
    patch([z0 z1 z1 z0], [0 0 1 1] + n, layer_color(n), ...
        "EdgeColor", [0.25 0.25 0.25], "FaceAlpha", 0.7);
    text(0.5 * (z0 + z1), n + 0.5, names(n), ...
        "HorizontalAlignment", "center", "FontSize", 8);
    hold on;
end
ylim([0.8 numel(names) + 1.3]);
xlim([params.vertical.zMin_nm params.vertical.zMax_nm]);
set(gca, "YTick", []);
xline(0, "Color", [0.2 0.6 0.2], "LineStyle", "--", "LineWidth", 1.0);
end

function draw_contact_cartoon(params)
rectangle("Position", [0.05 0.70 0.30 0.12], "FaceColor", [0.95 0.70 0.18], "EdgeColor", "k");
rectangle("Position", [0.65 0.70 0.30 0.12], "FaceColor", [0.95 0.70 0.18], "EdgeColor", "k");
rectangle("Position", [0.36 0.63 0.28 0.22], "FaceColor", [0.85 1.00 1.00], "EdgeColor", [0.2 0.6 0.6]);
rectangle("Position", [0.12 0.35 0.76 0.22], "FaceColor", [0.40 0.75 0.85], "EdgeColor", "k");
rectangle("Position", [0.05 0.17 0.90 0.10], "FaceColor", [0.95 0.70 0.18], "EdgeColor", "k");
text(0.20, 0.86, "split top gate/contact", "HorizontalAlignment", "center", "FontSize", 8);
text(0.80, 0.86, "same potential", "HorizontalAlignment", "center", "FontSize", 8);
text(0.50, 0.88, "metal-free NV window", "HorizontalAlignment", "center", "FontSize", 8);
text(0.50, 0.47, "electronic ARROW stack", "HorizontalAlignment", "center", "FontSize", 8);
text(0.50, 0.10, "continuous bottom contact", "HorizontalAlignment", "center", "FontSize", 8);
text(0.50, 0.02, sprintf("window width = %.3g um", params.layout.windowWidth_um), ...
    "HorizontalAlignment", "center", "FontSize", 8);
xlim([0 1]);
ylim([0 1]);
end

function c = layer_color(n)
palette = [
    0.80 0.86 0.90
    0.95 0.55 0.38
    0.48 0.68 0.86
    0.90 0.42 0.34
    0.20 0.75 0.82
    0.90 0.42 0.34
    0.48 0.68 0.86
    0.95 0.55 0.38
    0.80 0.86 0.90];
c = palette(1 + mod(n - 1, size(palette, 1)), :);
end
