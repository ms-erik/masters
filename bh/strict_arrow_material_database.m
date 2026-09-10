function db = strict_arrow_material_database()
%STRICT_ARROW_MATERIAL_DATABASE Lattice-matched InGaAlAs design model.
%
% All alloys keep the indium fraction close to 0.53 and exchange Ga for Al:
%   In0.53 Al_y Ga_(0.47-y) As.
% This keeps the first design approximately lattice matched to InP while
% providing a continuously tunable conduction-band edge.
%
% The quaternary interpolation below is a transparent design model, not a
% substitute for HRXRD/composition measurements or a calibrated k.p model.

db.substrate.name = "semi-insulating InP (001)";

db.InGaAs.name = "In0.53Ga0.47As";
db.InGaAs.alFraction = 0;
db.InGaAs.conductionEdge_meV = 0;
db.InGaAs.electronMass_m0 = 0.041;
db.InGaAs.staticPermittivity = 13.9;
db.InGaAs.gFactor = -4.3;

db.InAlAs.name = "In0.52Al0.48As";
db.InAlAs.alFraction = 0.48;
db.InAlAs.conductionEdge_meV = 510;
db.InAlAs.electronMass_m0 = 0.080;
db.InAlAs.staticPermittivity = 12.42;
db.InAlAs.gFactor = 0.5;

db.buffer = interpolate_alloy(0.08, db);
db.buffer.name = "In0.53Al0.08Ga0.39As";
db.buffer.role = "low-k propagating buffer";

db.guide = interpolate_alloy(0.10, db);
db.guide.name = "In0.53Al0.10Ga0.37As";
db.guide.role = "lower-k central electronic guide";

db.highK = db.InGaAs;
db.highK.role = "thin high-k antiresonant reflector";

db.interpolation = ...
    "linear interpolation from InGaAs to InAlAs at fixed In fraction";
db.references = [
    "Peng et al., J. Appl. Phys. 60, 1709 (1986), DOI 10.1063/1.337262"
    "Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001), DOI 10.1063/1.1368156"
];
end

function alloy = interpolate_alloy(alFraction, db)
fraction = alFraction / db.InAlAs.alFraction;
alloy.alFraction = alFraction;
alloy.conductionEdge_meV = fraction * db.InAlAs.conductionEdge_meV;
alloy.electronMass_m0 = db.InGaAs.electronMass_m0 + fraction * ...
    (db.InAlAs.electronMass_m0 - db.InGaAs.electronMass_m0);
alloy.staticPermittivity = db.InGaAs.staticPermittivity + fraction * ...
    (db.InAlAs.staticPermittivity - db.InGaAs.staticPermittivity);
alloy.gFactor = db.InGaAs.gFactor + fraction * ...
    (db.InAlAs.gFactor - db.InGaAs.gFactor);
end
