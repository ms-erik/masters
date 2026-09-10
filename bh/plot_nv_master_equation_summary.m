function plot_nv_master_equation_summary(master)
%PLOT_NV_MASTER_EQUATION_SUMMARY Diagnostic plots for the master equation.

labels = categorical(master.stateA.labels);
labels = reordercats(labels, cellstr(master.stateA.labels));

figure("Color", "w", "Name", "NV spin-charge master equation");
tiledlayout(2, 1, "TileSpacing", "compact");

nexttile;
bar(labels, [master.stateA.p master.stateB.p]);
ylabel("population");
title("Final populations after readout pulse");
legend("state A", "state B", "Location", "best");
grid on; box on;

nexttile;
bar(categorical(["P_{NV-} A", "P_{NV-} B", "\DeltaV_A", "\DeltaV_B"]), ...
    [master.P_NVminus_A, master.P_NVminus_B, master.deltaV_A_meV, master.deltaV_B_meV]);
ylabel("population or meV");
title(sprintf("\\DeltaP_{A-B}=%.3g,  \\DeltaV_{A-B}=%.3g meV", ...
    master.deltaP_AminusB, master.deltaV_AminusB_meV));
grid on; box on;
end
