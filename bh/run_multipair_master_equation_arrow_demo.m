%% Five-pair strict electronic ARROW with NV spin-charge master equation
%
% This is the new workflow.  It does not impose DeltaV_A and DeltaV_B as
% independent phenomenological parameters.  Instead it computes them from:
%
%   master equation -> Delta P_NVminus(m_s) -> Delta V0(m_s)
%
% and then reuses the existing 1D transport / Landauer pipeline.

clear; clc; close all;
packageDir = fileparts(mfilename('full11path'));
addpath(packageDir, "-begin");
cd(packageDir);

params = multipair_electronic_arrow_params("baseline");
[params, nvMaster] = update_spin_deltaV_from_master_equation(params);

fprintf("\nFive-pair strict electronic ARROW + master equation\n");
fprintf("Reflector films per side: %d\n", params.arrow.numberOfReflectorFilmsPerSide);
fprintf("Target energy: %.3f meV; estimated EF: %.3f meV\n", ...
    params.arrow.targetEnergy_meV, params.current.fermiEnergy_meV);
fprintf("Master-equation DeltaV_A = %.6g meV\n", params.spin.deltaV_A_meV);
fprintf("Master-equation DeltaV_B = %.6g meV\n", params.spin.deltaV_B_meV);
fprintf("DeltaV_AminusB = %.6g meV\n\n", nvMaster.deltaV_AminusB_meV);

run_splitgate_nv_arrow_demo;
diagnostics = strict_arrow_diagnostics(params, params.arrowSweep.enable);
plot_nv_master_equation_summary(nvMaster);

disp("Multipair strict-ARROW diagnostics:");
disp(rmfield(diagnostics, {'Tmap','Ttarget','guideFraction','mirrorReflectance'}));

save("Resultados_Multipair_MasterEquation_ARROW.mat", ...
    "params", "nvMaster", "devices", "energy", "coarse", "resGuess", "results", ...
    "probeEnergy", "wf", "wfDevices", "bestIdx", "diagnostics");
