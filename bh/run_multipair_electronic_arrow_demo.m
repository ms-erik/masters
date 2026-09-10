%% Five-pair strict electronic ARROW
clear; clc; close all;
packageDir = fileparts(mfilename('fullpath'));
addpath(packageDir,"-begin");
cd(packageDir);

params = multipair_electronic_arrow_params("baseline");

fprintf("\nFive-pair strict electronic ARROW\n");
fprintf("Reflector films per side: %d\n",params.arrow.numberOfReflectorFilmsPerSide);
fprintf("Target energy: %.3f meV; estimated EF: %.3f meV\n", ...
    params.arrow.targetEnergy_meV,params.current.fermiEnergy_meV);
fprintf("High-k film: %.4f nm; guide: %.4f nm\n", ...
    params.arrow.highKFilmThickness_nm,params.arrow.guideThickness_nm);
fprintf("Buffer/guide Al fractions: %.3f / %.3f\n\n", ...
    params.epitaxy.bufferAlFraction,params.epitaxy.guideAlFraction);

run_splitgate_nv_arrow_demo;
diagnostics = strict_arrow_diagnostics(params,params.arrowSweep.enable);

disp("Multipair strict-ARROW diagnostics:");
disp(rmfield(diagnostics,{'Tmap','Ttarget','guideFraction','mirrorReflectance'}));

save("Resultados_Multipair_Strict_Electronic_ARROW.mat", ...
    "params","devices","energy","coarse","resGuess","results", ...
    "probeEnergy","wf","wfDevices","bestIdx","diagnostics");
