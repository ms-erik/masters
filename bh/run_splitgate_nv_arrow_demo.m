%% Split-top-gate / continuous-bottom-contact NV-ARROW demo
% Geometry represented:
%   532 nm pump from above
%   removable diamond with a 2D NV layer on its lower face
%   metal-free central window below the NV layer
%   split top gate/contact pads at the same potential, outside the window
%   horizontal-layer electronic ARROW stack
%   continuous bottom contact
%
% The 1D transport coordinate z represents the vertical column under the
% metal-free window. The split top gate is included as an effective smooth
% electrostatic shift, not as a metal layer in the current path.

if ~exist("params", "var")
    clearvars; clc; close all;
    packageDir = fileparts(mfilename('fullpath'));
    addpath(packageDir, "-begin");
    cd(packageDir);
    params = multipair_electronic_arrow_params("baseline");
else
    clc; close all;
    packageDir = fileparts(mfilename('fullpath'));
    addpath(packageDir, "-begin");
    cd(packageDir);
end

devices.none = generate_splitgate_nv_arrow_structure(params, "none");
devices.spinA = generate_splitgate_nv_arrow_structure(params, "A");
devices.spinB = generate_splitgate_nv_arrow_structure(params, "B");

energy = linspace(params.energy.min_meV, params.energy.max_meV, params.energy.nCoarse);
coarse = solve_transport_spectrum(devices, energy, params);
resGuess = detect_resonances(coarse.energy_meV, coarse.spinA.T, params.analysis);
energy = refine_energy_grid(energy, resGuess, params.energy);

results = solve_transport_spectrum(devices, energy, params);
results.resonances = detect_resonances(results.energy_meV, results.spinA.T, params.analysis);
results.antiresonances = detect_resonances(results.energy_meV, 1 - results.spinA.T, params.analysis);

if ~isempty(results.resonances.energy_meV)
    [~, bestIdx] = max(results.resonances.Q);
    probeEnergy = results.resonances.energy_meV(bestIdx);
else
    [~, bestIdx] = max(results.spinA.T);
    probeEnergy = results.energy_meV(bestIdx);
end

wfDevices = devices;
if isfield(params, "zeeman") && params.zeeman.enable && norm(params.zeeman.B_T) > 0
    if lower(string(params.zeeman.injection)) == "down"
        zeemanProjection = -1;
    else
        zeemanProjection = +1;
    end
    wfDevices.none = apply_zeeman_channel(devices.none, params, zeemanProjection);
    wfDevices.spinA = apply_zeeman_channel(devices.spinA, params, zeemanProjection);
    wfDevices.spinB = apply_zeeman_channel(devices.spinB, params, zeemanProjection);
end

wf.none = wavefunction_profile(wfDevices.none, probeEnergy, params);
wf.spinA = wavefunction_profile(wfDevices.spinA, probeEnergy, params);
wf.spinB = wavefunction_profile(wfDevices.spinB, probeEnergy, params);

physicalFermiEnergy = params.current.fermiEnergy_meV;
autoAlignToProbe = isfield(params.current, "autoAlignToProbe") && ...
    params.current.autoAlignToProbe;
if autoAlignToProbe
    params.current.fermiEnergy_meV = probeEnergy;
end
results.current = landauer_current(results.energy_meV, ...
    results.none.T, results.spinA.T, results.spinB.T, ...
    params.current.bias_V, params);
paramsAtProbe = params;
paramsAtProbe.current.fermiEnergy_meV = probeEnergy;
results.currentAtProbe = landauer_current(results.energy_meV, ...
    results.none.T, results.spinA.T, results.spinB.T, ...
    params.current.bias_V, paramsAtProbe);
results.current.physicalFermiEnergy_meV = physicalFermiEnergy;
results.current.probeEnergy_meV = probeEnergy;
results.nvZeeman = nv_center_zeeman(params);
results.metrics = extract_spin_sensitivity_metrics(results, params);

plot_splitgate_nv_arrow_results(devices, results, wf, params);

if params.sweep.enableDistanceSweep
    sweep = parameter_sweep_splitgate_nv_distance(params, params.sweep.nvDistanceValues_nm);
    figure("Color", "w", "Name", "Split-gate NV distance sweep");
    imagesc(sweep.energy_meV, sweep.distance_nm, sweep.deltaT);
    axis xy; box on;
    xlabel("Energy (meV)");
    ylabel("NV layer to stack distance (nm)");
    title("\DeltaT(E,d_{NV}) under metal-free window");
    cb = colorbar;
    cb.Label.String = "\DeltaT";
    colormap(turbo);
end

disp("Split-top-gate / continuous-bottom-contact NV-ARROW demo complete.");
disp("Detected spin-A resonances:");
disp(struct2table(results.resonances));
disp("Spin-sensitivity metrics:");
disp(results.metrics);
if ~isempty(results.metrics.resonanceTable)
    disp("Matched spin-dependent resonance shifts:");
    disp(results.metrics.resonanceTable);
end
if isfield(results, "nvZeeman") && results.nvZeeman.enabled
    disp("NV-center Zeeman levels and ODMR frequencies:");
    disp(struct2table(rmfield(results.nvZeeman, ...
        {'axis', 'ms', 'frequencyFromMs0_GHz', 'energyFromMs0_meV'})));
end
