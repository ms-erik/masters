function sweep = parameter_sweep_splitgate_nv_distance(params, distanceValues_nm)
%PARAMETER_SWEEP_SPLITGATE_NV_DISTANCE Sweep NV-stack distance.

energy = linspace(params.energy.min_meV, params.energy.max_meV, min(params.energy.nCoarse, 700));
deltaT = zeros(numel(distanceValues_nm), numel(energy));

p = params;
p.sweep.enableDistanceSweep = false;

for n = 1:numel(distanceValues_nm)
    p.spin.distanceFromStack_nm = distanceValues_nm(n);
    devices.none = generate_splitgate_nv_arrow_structure(p, "none");
    devices.spinA = generate_splitgate_nv_arrow_structure(p, "A");
    devices.spinB = generate_splitgate_nv_arrow_structure(p, "B");
    result = solve_transport_spectrum(devices, energy, p);
    deltaT(n, :) = result.deltaT_AminusB;
end

sweep.energy_meV = energy;
sweep.distance_nm = distanceValues_nm;
sweep.deltaT = deltaT;
end

