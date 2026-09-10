function diagnostics = strict_arrow_diagnostics(params, makePlots)
%STRICT_ARROW_DIAGNOSTICS Test whether leakage follows ARROW antiresonance.
%
% The decisive test is a sweep of the thin high-k film thickness. At fixed
% energy, leakage transmission should show minima close to:
%   k_high d = (m + 1/2) pi.

if nargin < 2
    makePlots = true;
end

E0 = params.arrow.targetEnergy_meV;
c = params.constants;
db = params.material.database;

kHigh = local_k(E0, db.highK.conductionEdge_meV, db.highK.electronMass_m0, c);
kBuffer = local_k(E0, db.buffer.conductionEdge_meV, db.buffer.electronMass_m0, c);
kGuide = local_k(E0, db.guide.conductionEdge_meV, db.guide.electronMass_m0, c);

diagnostics.targetEnergy_meV = E0;
diagnostics.kHigh_nm = kHigh;
diagnostics.kBuffer_nm = kBuffer;
diagnostics.kGuide_nm = kGuide;
diagnostics.kReservoir_nm = kHigh;
diagnostics.kOrderingPass = kGuide < kBuffer && kBuffer < kHigh;
diagnostics.allLayersPropagating = E0 > max(params.vertical.potential_meV);
diagnostics.highKPhase_rad = kHigh*params.arrow.highKFilmThickness_nm;
diagnostics.guidePhase_rad = kGuide*params.arrow.guideThickness_nm;
diagnostics.highKPhaseOverPi = diagnostics.highKPhase_rad/pi;
diagnostics.guidePhaseOverPi = diagnostics.guidePhase_rad/pi;

dVals = params.arrowSweep.highKThickness_nm;
energies = params.arrowSweep.energy_meV;
Tmap = zeros(numel(dVals),numel(energies));
Ttarget = zeros(size(dVals));
guideFraction = zeros(size(dVals));
mirrorReflectance = zeros(size(dVals));

for id = 1:numel(dVals)
    p = local_set_high_k_thickness(params,dVals(id));
    device = generate_splitgate_nv_arrow_structure(p,"none");
    for ie = 1:numel(energies)
        out = tm_transport(device,energies(ie),p);
        Tmap(id,ie) = out.T;
    end
    out = tm_transport(device,E0,p);
    Ttarget(id) = out.T;
    wf = wavefunction_profile(device,E0,p);
    guideIdx = p.arrow.guideLayerIndex;
    z0 = p.vertical.boundaries_nm(guideIdx);
    z1 = p.vertical.boundaries_nm(guideIdx + 1);
    stackMask = wf.x_nm >= 0 & wf.x_nm <= p.vertical.stackThickness_nm;
    guideMask = wf.x_nm >= z0 & wf.x_nm <= z1;
    guideFraction(id) = trapz(wf.x_nm(guideMask),wf.probabilityDensity(guideMask)) / ...
        max(trapz(wf.x_nm(stackMask),wf.probabilityDensity(stackMask)),eps);
    mirrorReflectance(id) = local_single_film_reflectance( ...
        E0,dVals(id),db.buffer,db.highK,params);
end

diagnostics.highKThickness_nm = dVals;
diagnostics.energy_meV = energies;
diagnostics.Tmap = Tmap;
diagnostics.Ttarget = Ttarget;
diagnostics.guideFraction = guideFraction;
diagnostics.mirrorReflectance = mirrorReflectance;

maxOrder = floor(kHigh*max(dVals)/pi - 0.5);
orders = 0:max(maxOrder,0);
predicted = (orders+0.5)*pi/kHigh;
predicted = predicted(predicted >= min(dVals) & predicted <= max(dVals));
diagnostics.predictedAntiresonantThickness_nm = predicted;

[~,minIdx] = local_minima(Ttarget);
diagnostics.observedLeakageMinima_nm = dVals(minIdx);
[~,mirrorMaxIdx] = local_maxima(mirrorReflectance);
diagnostics.observedMirrorMaxima_nm = dVals(mirrorMaxIdx);
if isempty(predicted) || isempty(mirrorMaxIdx)
    diagnostics.meanMirrorPeakError_nm = NaN;
    diagnostics.mirrorAntiresonancePass = false;
else
    errors = arrayfun(@(x) min(abs(dVals(mirrorMaxIdx)-x)),predicted);
    diagnostics.meanMirrorPeakError_nm = mean(errors);
    diagnostics.mirrorAntiresonancePass = diagnostics.meanMirrorPeakError_nm < ...
        median(diff(dVals))*2.5;
end
diagnostics.wholeDeviceLeakageModulation = ...
    (max(Ttarget)-min(Ttarget))/max(max(Ttarget),eps);
diagnostics.strictArrowPass = diagnostics.kOrderingPass && ...
    diagnostics.allLayersPropagating && diagnostics.mirrorAntiresonancePass && ...
    diagnostics.wholeDeviceLeakageModulation > 0.1;

if makePlots
    figure("Color","w","Name","Strict electronic ARROW diagnostics");
    tiledlayout(3,1,"TileSpacing","compact");
    nexttile;
    imagesc(energies,dVals,Tmap); axis xy; colorbar;
    xlabel("Energy (meV)"); ylabel("high-k film thickness (nm)");
    title("Leakage map T(E,d_{AR})");
    nexttile;
    yyaxis left;
    semilogy(dVals,max(Ttarget,1e-14),"k-","LineWidth",1.3); hold on;
    ylabel("whole-device leakage T");
    yyaxis right;
    plot(dVals,mirrorReflectance,"b-","LineWidth",1.2);
    ylabel("isolated-film reflectance");
    for x = predicted
        xline(x,"r--");
    end
    xlabel("high-k film thickness (nm)");
    title(sprintf("Mirror antiresonance and complete-device leakage at E = %.1f meV",E0));
    grid on;
    nexttile;
    plot(dVals,guideFraction,"b-","LineWidth",1.3); hold on;
    for x = predicted
        xline(x,"r--");
    end
    xlabel("high-k film thickness (nm)"); ylabel("guide probability fraction");
    title("Central-guide localization"); grid on;
end
end

function p = local_set_high_k_thickness(p,d_nm)
p.vertical.thickness_nm(p.arrow.highKLayerIndices) = d_nm;
p.vertical.boundaries_nm = [0 cumsum(p.vertical.thickness_nm)];
p.vertical.stackThickness_nm = sum(p.vertical.thickness_nm);
p.vertical.zMax_nm = p.vertical.stackThickness_nm + 8;
end

function k_nm = local_k(E,V,m,c)
k_nm = sqrt(2*m*c.m0_kg*(E-V)*1e-3*c.q_C)/c.hbar_Js*1e-9;
end

function [values,idx] = local_minima(y)
idx = find(y(2:end-1) < y(1:end-2) & y(2:end-1) <= y(3:end))+1;
values = y(idx);
end

function [values,idx] = local_maxima(y)
idx = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end))+1;
values = y(idx);
end

function R = local_single_film_reflectance(E,d,buffer,highK,params)
layers.width_nm = d;
layers.V_meV = highK.conductionEdge_meV;
layers.mEff = highK.electronMass_m0;
layers.gFactor = highK.gFactor;
layers.Vleft_meV = buffer.conductionEdge_meV;
layers.Vright_meV = buffer.conductionEdge_meV;
layers.mLeft = buffer.electronMass_m0;
layers.mRight = buffer.electronMass_m0;
layers.gLeft = buffer.gFactor;
layers.gRight = buffer.gFactor;
device.layers = layers;
out = tm_transport(device,E,params);
R = out.R;
end
