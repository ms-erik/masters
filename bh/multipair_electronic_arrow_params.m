function params = multipair_electronic_arrow_params(preset)
%MULTIPAIR_ELECTRONIC_ARROW_PARAMS Five antiresonant films per side.
%
% This is a strict propagating-layer ARROW/Bragg hybrid. The high-k InGaAs
% films are not potential barriers. They are phase reflectors embedded in
% lower-k InGaAlAs buffers. All active layers remain propagating at the
% nominal operating energy.

if nargin < 1
    preset = "baseline";
end
preset = lower(string(preset));
db = strict_arrow_material_database();

params.constants.hbar_Js = 1.054571817e-34;
params.constants.q_C = 1.602176634e-19;
params.constants.m0_kg = 9.1093837015e-31;
params.constants.kB_JK = 1.380649e-23;

params.material.name = ...
    "five-pair strict electronic ARROW: InGaAs/InGaAlAs on InP";
params.material.substrate = db.substrate.name;
params.material.temperature_K = 4.2;
params.material.massModel = "layer-resolved constant Gamma mass";
params.material.bandModel = "propagating quaternary layers";
params.material.database = db;
params.material.designStatus = ...
    "MBE-oriented phase design; quaternary band line-up requires calibration";

params.arrow.targetEnergy_meV = 150;
params.arrow.numberOfReflectorFilmsPerSide = 5;
params.arrow.phaseOrderHighK = 0;
params.arrow.phaseOrderGuide = 2;
params.arrow.definition = ...
    "low-k guide with five thin high-k antiresonant films per side";

kHigh = local_k_nm(params.arrow.targetEnergy_meV, ...
    db.highK.conductionEdge_meV, db.highK.electronMass_m0, params.constants);
kGuide = local_k_nm(params.arrow.targetEnergy_meV, ...
    db.guide.conductionEdge_meV, db.guide.electronMass_m0, params.constants);
dHigh_nm = (params.arrow.phaseOrderHighK + 0.5)*pi/kHigh;
dGuide_nm = params.arrow.phaseOrderGuide*pi/kGuide;

if preset == "fast"
    buffer_nm = 4.5;
    innerBuffer_nm = 5.0;
else
    buffer_nm = 5.5;
    innerBuffer_nm = 6.0;
end

nPairs = params.arrow.numberOfReflectorFilmsPerSide;
names = "n+ InGaAs top access";
types = "high";
thickness = 3;

for n = 1:nPairs
    names(end+1) = "InGaAlAs buffer L" + n; %#ok<AGROW>
    types(end+1) = "buffer"; %#ok<AGROW>
    thickness(end+1) = buffer_nm; %#ok<AGROW>
    names(end+1) = "thin InGaAs AR film L" + n; %#ok<AGROW>
    types(end+1) = "high"; %#ok<AGROW>
    thickness(end+1) = dHigh_nm; %#ok<AGROW>
end

names(end+1) = "InGaAlAs inner buffer L";
types(end+1) = "buffer";
thickness(end+1) = innerBuffer_nm;
names(end+1) = "InGaAlAs low-k central guide";
types(end+1) = "guide";
thickness(end+1) = dGuide_nm;
guideIndex = numel(names);
names(end+1) = "InGaAlAs inner buffer R";
types(end+1) = "buffer";
thickness(end+1) = innerBuffer_nm;

for n = nPairs:-1:1
    names(end+1) = "thin InGaAs AR film R" + n; %#ok<AGROW>
    types(end+1) = "high"; %#ok<AGROW>
    thickness(end+1) = dHigh_nm; %#ok<AGROW>
    names(end+1) = "InGaAlAs buffer R" + n; %#ok<AGROW>
    types(end+1) = "buffer"; %#ok<AGROW>
    thickness(end+1) = buffer_nm; %#ok<AGROW>
end

names(end+1) = "n+ InGaAs bottom access";
types(end+1) = "high";
thickness(end+1) = 8;

params.vertical.layerNames = names;
params.vertical.layerTypes = types;
params.vertical.thickness_nm = thickness;
params.arrow.guideLayerIndex = guideIndex;
params.arrow.highKLayerIndices = find(types == "high" & ...
    ~contains(names,"access"));

[Vlayers, mlayers, glayers, materialNames] = local_profiles(types, db);
params.vertical.potential_meV = ...
    [db.highK.conductionEdge_meV Vlayers db.highK.conductionEdge_meV];
params.vertical.massEff = ...
    [db.highK.electronMass_m0 mlayers db.highK.electronMass_m0];
params.vertical.gFactor = [db.highK.gFactor glayers db.highK.gFactor];
params.vertical.materialNames = ["InGaAs lead" materialNames "InGaAs lead"];
params.vertical.boundaries_nm = [0 cumsum(thickness)];
params.vertical.stackThickness_nm = sum(thickness);
params.vertical.zMin_nm = -8;
params.vertical.zMax_nm = params.vertical.stackThickness_nm + 8;
params.vertical.dz_nm = 0.10;

params.arrow.highKFilmThickness_nm = dHigh_nm;
params.arrow.guideThickness_nm = dGuide_nm;
params.arrow.bufferPotential_meV = db.buffer.conductionEdge_meV;
params.arrow.guidePotential_meV = db.guide.conductionEdge_meV;
params.arrow.highKPotential_meV = db.highK.conductionEdge_meV;
params.arrow.allLayersPropagatingAtTarget = ...
    params.arrow.targetEnergy_meV > max(params.vertical.potential_meV);

params.leads.leftPotential_meV = db.highK.conductionEdge_meV;
params.leads.rightPotential_meV = db.highK.conductionEdge_meV;

params.epitaxy.substrate = db.substrate.name;
params.epitaxy.growthMethod = "MBE";
params.epitaxy.contactDopingSpecies = "Si";
params.epitaxy.topAccessDoping_cm3 = 2.2e18;
params.epitaxy.bottomAccessDoping_cm3 = 2.2e18;
params.epitaxy.activeRegionBackground_cm3 = 1e15;
params.epitaxy.quaternaryModel = db.interpolation;
params.epitaxy.bufferAlFraction = db.buffer.alFraction;
params.epitaxy.guideAlFraction = db.guide.alFraction;
params.epitaxy.bandEdgeUncertainty_meV = 20;

params.diamond.thickness_um = 50;
params.diamond.nvGeometry = "2D shallow-NV layer on removable diamond face";
params.diamond.pumpWavelength_nm = 532;
params.diamond.isMetallized = false;
params.diamond.epsilon_r = 5.7;

params.layout.topContactType = "split lateral top contact outside window";
params.layout.bottomContactType = "continuous bottom contact";
params.layout.metalFreeWindow = true;
params.layout.windowWidth_um = 50;
params.layout.topPadSeparation_um = 50;
params.layout.topPadVoltage_V = 0;
params.layout.bottomContactVoltage_V = 0;
params.layout.electrostaticModel = "flat-band multipair ARROW baseline";

params.gate.enableSplitGateShift = true;
params.gate.topGateShift_meV = 0;
params.gate.screeningLength_nm = 20;
params.gate.windowAttenuation = 0.25;
params.gate.model = "analytic";
params.gate.poisson.xHalfWidth_um = 80;
params.gate.poisson.dx_um = 2;
params.gate.poisson.dz_nm = 1;
params.gate.poisson.maxIter = 4000;
params.gate.poisson.tolerance_V = 1e-6;
params.gate.poisson.relaxation = 1.75;
params.gate.leverArm = 1;

params.spin.deltaV_A_meV = 0.9;
params.spin.deltaV_B_meV = -0.9;
params.spin.distanceFromStack_nm = 2;
params.spin.sheetThickness_nm = 2;
params.spin.screeningLength_nm = 12;
params.spin.distancePower = 1;
params.spin.distanceRegularization_nm = 1;
params.spin.interfaceCouplingMode = "permittivity";
params.spin.interfaceCoupling = 0.58;
params.spin.epsilonDiamond = 5.7;
params.spin.epsilonSemiconductor = db.highK.staticPermittivity;
params.spin.lateralFillFactor = 1;
params.spin.phaseAmplitude_rad = 0;
params.spin.calibrationMode = "two-stage interface coupling";
params.spin.surfaceDensity_cm2 = 1e11;
params.spin.effectiveDipole_D = 1;
params.spin.stateA_label = "NV m_s = 0";
params.spin.stateB_label = "NV m_s = -1";
params.spin.stateA_ms = 0;
params.spin.stateB_ms = -1;

% New spin-charge master-equation block.
% The legacy fields deltaV_A_meV and deltaV_B_meV above are still kept as a
% manual fallback.  In the master-equation workflow they are overwritten by
% update_spin_deltaV_from_master_equation.m according to
%   m_s -> Delta P_NVminus(m_s) -> Delta V0(m_s).
params.spin.masterEquation.enable = true;
params.spin.masterEquation.mode = "finite_pulse";       % "finite_pulse" or "steady_state"
params.spin.masterEquation.readoutPulse_s = 300e-9;
params.spin.masterEquation.referenceMode = "midpoint";  % keeps A/B symmetric around no-NV baseline
params.spin.masterEquation.energySign = +1;             % +1 for electron potential-energy convention
params.spin.masterEquation.capacitanceMode = "parallel_plate";
params.spin.masterEquation.interfaceDielectricThickness_nm = 5.0;
params.spin.masterEquation.interfaceEpsilon_r = 8.0;
params.spin.masterEquation.capacitance_F_m2 = NaN;       % used when capacitanceMode="manual"
params.spin.masterEquation.initialSpinPurity = 0.95;
params.spin.masterEquation.initialChargeNVminus = 0.98;
params.spin.masterEquation.verbose = true;

% Effective rate constants for a compact NV-/NV0 spin-charge model.
% These are deliberately grouped here because they should be calibrated
% against measured charge-readout/PDMR data, not treated as immutable
% material constants.
params.spin.masterEquation.rates.k_ex_532_s = 60e6;      % optical excitation, NV-
params.spin.masterEquation.rates.k_rad_s = 77e6;         % triplet radiative decay
params.spin.masterEquation.rates.k_isc_0_s = 4e6;        % E0 -> singlet
params.spin.masterEquation.rates.k_isc_1_s = 45e6;       % E+-1 -> singlet
params.spin.masterEquation.rates.k_s_to_0_s = 3.5e6;     % singlet -> G0
params.spin.masterEquation.rates.k_s_to_1_s = 0.7e6;     % singlet -> G+-1
params.spin.masterEquation.rates.k_ion_0_s = 1.0e6;      % E0 -> NV0
params.spin.masterEquation.rates.k_ion_1_s = 5.0e6;      % E+-1 -> NV0
params.spin.masterEquation.rates.k_ion_singlet_s = 1.0e6;% singlet -> NV0
params.spin.masterEquation.rates.k_rec_to_0_s = 1.0e6;   % NV0 -> G0
params.spin.masterEquation.rates.k_rec_to_1_s = 0.2e6;   % NV0 -> G+-1
params.spin.masterEquation.rates.k_mw_mix_s = 0;         % optional microwave mixing during readout

params.nvZeeman.enable = true;
params.nvZeeman.B_T = [0 0 0.05];
params.nvZeeman.axis = [0 0 1];
params.nvZeeman.D_GHz = 2.870;
params.nvZeeman.gamma_GHz_T = 28.02495164;
params.nvZeeman.includeTransverseMixing = false;
params.nvZeeman.stateA_ms = 0;
params.nvZeeman.stateB_ms = -1;

params.zeeman.enable = false;
params.zeeman.B_T = [0 0 0];
params.zeeman.injection = "unpolarized";
params.zeeman.computeSpinResolved = true;
params.disorder.enable = false;
params.disorder.rms_meV = 0.15;
params.disorder.correlationLength_nm = 2;
params.disorder.seed = 7;
params.decoherence.enable = false;
params.decoherence.gamma_meV = 0.02;
params.decoherence.region = "stack";

params.energy.min_meV = 110;
params.energy.max_meV = 190;
params.energy.nCoarse = 1300;
params.energy.refinePointsPerWindow = 900;
params.energy.refineHalfWindow_meV = 0.5;
if preset == "fast"
    params.energy.nCoarse = 380;
    params.energy.refinePointsPerWindow = 260;
    params.vertical.dz_nm = 0.16;
end

params.analysis.minPeakProminence = 0.003;
params.analysis.minPeakDistance_meV = 0.08;
params.analysis.fwhmFloor_meV = 1e-6;

params.current.contactFermiEstimate_meV = local_fermi_meV( ...
    params.epitaxy.topAccessDoping_cm3,db.highK.electronMass_m0,params.constants);
params.current.fermiEnergy_meV = 171.6;                                     % params.current.contactFermiEstimate_meV;
params.current.autoAlignToProbe = false;
params.current.bias_V = 1e-3;
params.current.temperature_K = 4.2;

params.numeric.minAbsK_nm = 1e-12;
params.numeric.interfaceRegularization = 1e-14;
params.numeric.maxTransmission = 1e6;
params.sweep.enableDistanceSweep = false;
params.sweep.nvDistanceValues_nm = linspace(1,20,20);
params.arrowSweep.enable = true;
params.arrowSweep.highKThickness_nm = linspace(2.0,7.0,31);
params.arrowSweep.energy_meV = linspace(125,175,160);
end

function [V,m,g,names] = local_profiles(types,db)
V = zeros(size(types));
m = zeros(size(types));
g = zeros(size(types));
names = strings(size(types));
for i = 1:numel(types)
    switch types(i)
        case "buffer"
            a = db.buffer;
        case "guide"
            a = db.guide;
        otherwise
            a = db.highK;
    end
    V(i) = a.conductionEdge_meV;
    m(i) = a.electronMass_m0;
    g(i) = a.gFactor;
    names(i) = a.name;
end
end

function k_nm = local_k_nm(E,V,m,c)
k_nm = sqrt(2*m*c.m0_kg*(E-V)*1e-3*c.q_C)/c.hbar_Js*1e-9;
end

function EF = local_fermi_meV(n_cm3,m,c)
n = n_cm3*1e6;
EF = c.hbar_Js^2*(3*pi^2*n)^(2/3)/(2*m*c.m0_kg)/c.q_C*1e3;
end
