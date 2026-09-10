function resonances = detect_resonances(energy_meV, signal, analysis)
%DETECT_RESONANCES Locate resonances and estimate FWHM/Q.

energy_meV = energy_meV(:);
signal = signal(:);

if numel(energy_meV) < 5
    resonances = empty_resonance();
    return;
end

try
    [pks, locs, widths, prominences] = findpeaks(signal, energy_meV, ...
        "MinPeakProminence", analysis.minPeakProminence, ...
        "MinPeakDistance", analysis.minPeakDistance_meV, ...
        "WidthReference", "halfheight");
catch
    [pks, idx] = local_findpeaks(signal, analysis.minPeakProminence);
    locs = energy_meV(idx);
    widths = local_widths(energy_meV, signal, idx);
    prominences = pks - min(signal);
end

if isempty(locs)
    resonances = empty_resonance();
    return;
end

widths = max(widths(:), analysis.fwhmFloor_meV);
resonances.energy_meV = locs(:);
resonances.peakT = pks(:);
resonances.fwhm_meV = widths(:);
resonances.Q = locs(:) ./ widths(:);
resonances.prominence = prominences(:);
end

function r = empty_resonance()
r.energy_meV = [];
r.peakT = [];
r.fwhm_meV = [];
r.Q = [];
r.prominence = [];
end

function [pks, idx] = local_findpeaks(y, minProm)
idx = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
pks = y(idx);
keep = pks - min(y) >= minProm;
idx = idx(keep);
pks = pks(keep);
end

function widths = local_widths(x, y, idx)
widths = zeros(size(idx));
for n = 1:numel(idx)
    half = 0.5 * y(idx(n));
    left = find(y(1:idx(n)) <= half, 1, "last");
    rightRel = find(y(idx(n):end) <= half, 1, "first");
    if isempty(left) || isempty(rightRel)
        widths(n) = median(diff(x));
    else
        right = idx(n) + rightRel - 1;
        widths(n) = x(right) - x(left);
    end
end
end

