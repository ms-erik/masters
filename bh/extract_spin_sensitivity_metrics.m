function metrics = extract_spin_sensitivity_metrics(results, params)
%EXTRACT_SPIN_SENSITIVITY_METRICS Quantify whether the device senses spin.
%
% Returns resonance energies for no-spin, spin A, spin B, their splitting,
% max(abs(DeltaT)), and current contrast.

resNone = detect_resonances(results.energy_meV, results.none.T, params.analysis);
resA = detect_resonances(results.energy_meV, results.spinA.T, params.analysis);
resB = detect_resonances(results.energy_meV, results.spinB.T, params.analysis);

metrics.maxAbsDeltaT = max(abs(results.deltaT_AminusB));
[metrics.maxAbsDeltaTValue, idx] = max(abs(results.deltaT_AminusB));
metrics.energyAtMaxAbsDeltaT_meV = results.energy_meV(idx);

metrics.nResNone = numel(resNone.energy_meV);
metrics.nResA = numel(resA.energy_meV);
metrics.nResB = numel(resB.energy_meV);

metrics.resonanceTable = match_spin_resonances(resNone, resA, resB);

if isfield(results, "current")
    metrics.currentNone_A = results.current.none_A;
    metrics.currentA_A = results.current.spinA_A;
    metrics.currentB_A = results.current.spinB_A;
    metrics.deltaI_AminusB_A = results.current.deltaI_AminusB_A;
    denom = max(abs([results.current.spinA_A, results.current.spinB_A, eps]));
    metrics.relativeCurrentContrast = metrics.deltaI_AminusB_A / denom;
else
    metrics.currentNone_A = NaN;
    metrics.currentA_A = NaN;
    metrics.currentB_A = NaN;
    metrics.deltaI_AminusB_A = NaN;
    metrics.relativeCurrentContrast = NaN;
end

if isfield(results, "zeeman") && results.zeeman.enabled
    metrics.zeemanEnabled = true;
    metrics.zeemanBmag_T = results.zeeman.Bmag_T;
    metrics.maxAbsDeltaT_up = max(abs(results.zeeman.deltaT_up_AminusB));
    metrics.maxAbsDeltaT_down = max(abs(results.zeeman.deltaT_down_AminusB));
else
    metrics.zeemanEnabled = false;
    metrics.zeemanBmag_T = 0;
    metrics.maxAbsDeltaT_up = NaN;
    metrics.maxAbsDeltaT_down = NaN;
end

if isfield(results, "nvZeeman") && results.nvZeeman.enabled
    metrics.nvZeemanEnabled = true;
    metrics.nvBparallel_T = results.nvZeeman.Bparallel_T;
    metrics.nvODMR_0_to_minus_GHz = results.nvZeeman.f_0_to_minus_GHz;
    metrics.nvODMR_0_to_plus_GHz = results.nvZeeman.f_0_to_plus_GHz;
    metrics.nvStateA_ms = results.nvZeeman.stateA_ms;
    metrics.nvStateB_ms = results.nvZeeman.stateB_ms;
    metrics.nvDeltaFrequency_AminusB_GHz = results.nvZeeman.deltaFrequency_AminusB_GHz;
else
    metrics.nvZeemanEnabled = false;
    metrics.nvBparallel_T = NaN;
    metrics.nvODMR_0_to_minus_GHz = NaN;
    metrics.nvODMR_0_to_plus_GHz = NaN;
    metrics.nvStateA_ms = NaN;
    metrics.nvStateB_ms = NaN;
    metrics.nvDeltaFrequency_AminusB_GHz = NaN;
end

metrics.passesQualitativeSpinSensing = metrics.maxAbsDeltaT > 0.01 && ...
    ~isempty(metrics.resonanceTable);
end

function tableOut = match_spin_resonances(resNone, resA, resB)
if isempty(resA.energy_meV) || isempty(resB.energy_meV)
    tableOut = table();
    return;
end

rows = [];
for ia = 1:numel(resA.energy_meV)
    eA = resA.energy_meV(ia);
    [dB, ib] = min(abs(resB.energy_meV - eA));
    if isempty(dB)
        continue;
    end
    if isempty(resNone.energy_meV)
        eNone = NaN;
        dNone = NaN;
        qNone = NaN;
    else
        [dNone, in] = min(abs(resNone.energy_meV - eA));
        eNone = resNone.energy_meV(in);
        qNone = resNone.Q(in);
    end
    rows = [rows; eA, resB.energy_meV(ib), eNone, ...
        eA - resB.energy_meV(ib), dB, dNone, ...
        resA.Q(ia), resB.Q(ib), qNone]; %#ok<AGROW>
end

if isempty(rows)
    tableOut = table();
else
    tableOut = array2table(rows, "VariableNames", ...
        ["E_A_meV", "E_B_meV", "E_none_meV", "DeltaE_AminusB_meV", ...
        "NearestABDistance_meV", "NearestNoneDistance_meV", ...
        "Q_A", "Q_B", "Q_none"]);
end
end
