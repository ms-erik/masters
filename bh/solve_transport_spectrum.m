function results = solve_transport_spectrum(devices, energy_meV, params)
%SOLVE_TRANSPORT_SPECTRUM Compute T(E), R(E), and DeltaT.

energy_meV = energy_meV(:).';
names = string(fieldnames(devices));
results.energy_meV = energy_meV;

for s = 1:numel(names)
    name = char(names(s));
    if zeeman_enabled(params)
        upDevice = apply_zeeman_channel(devices.(name), params, +1);
        downDevice = apply_zeeman_channel(devices.(name), params, -1);
        up = solve_one_device(upDevice, energy_meV, params);
        down = solve_one_device(downDevice, energy_meV, params);
        switch lower(string(params.zeeman.injection))
            case "up"
                combined = up;
            case "down"
                combined = down;
            otherwise
                combined.T = 0.5 * (up.T + down.T);
                combined.R = 0.5 * (up.R + down.R);
                combined.detM = 0.5 * (up.detM + down.detM);
        end
        results.(name).T = real(combined.T);
        results.(name).R = real(combined.R);
        results.(name).detM = combined.detM;
        results.(name).zeeman.up = up;
        results.(name).zeeman.down = down;
    else
        solved = solve_one_device(devices.(name), energy_meV, params);
        results.(name).T = real(solved.T);
        results.(name).R = real(solved.R);
        results.(name).detM = solved.detM;
    end
end

results.deltaT_AminusB = results.spinA.T - results.spinB.T;
results.deltaT_AminusNone = results.spinA.T - results.none.T;
results.deltaT_BminusNone = results.spinB.T - results.none.T;

if zeeman_enabled(params)
    results.zeeman.enabled = true;
    results.zeeman.B_T = params.zeeman.B_T;
    results.zeeman.Bmag_T = norm(params.zeeman.B_T);
    results.zeeman.injection = params.zeeman.injection;
    results.zeeman.deltaT_up_AminusB = results.spinA.zeeman.up.T - results.spinB.zeeman.up.T;
    results.zeeman.deltaT_down_AminusB = results.spinA.zeeman.down.T - results.spinB.zeeman.down.T;
else
    results.zeeman.enabled = false;
end
end

function solved = solve_one_device(device, energy_meV, params)
T = zeros(size(energy_meV));
R = zeros(size(energy_meV));
detM = zeros(size(energy_meV));
for n = 1:numel(energy_meV)
    out = tm_transport(device, energy_meV(n), params);
    T(n) = out.T;
    R(n) = out.R;
    detM(n) = out.detM;
end
solved.T = real(T);
solved.R = real(R);
solved.detM = detM;
end

function tf = zeeman_enabled(params)
tf = isfield(params, "zeeman") && params.zeeman.enable && ...
    isfield(params.zeeman, "computeSpinResolved") && params.zeeman.computeSpinResolved && ...
    norm(params.zeeman.B_T) > 0;
end
