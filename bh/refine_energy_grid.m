function energy = refine_energy_grid(baseEnergy, resonances, energyParams)
%REFINE_ENERGY_GRID Add dense windows around candidate resonances.

energy = baseEnergy(:).';
if isempty(resonances.energy_meV)
    return;
end

for n = 1:numel(resonances.energy_meV)
    e0 = resonances.energy_meV(n);
    local = linspace(e0 - energyParams.refineHalfWindow_meV, ...
        e0 + energyParams.refineHalfWindow_meV, energyParams.refinePointsPerWindow);
    local = local(local >= energyParams.min_meV & local <= energyParams.max_meV);
    energy = [energy local]; %#ok<AGROW>
end

energy = unique(sort(energy));
end

