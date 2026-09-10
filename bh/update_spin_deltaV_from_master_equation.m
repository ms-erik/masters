function [params, master] = update_spin_deltaV_from_master_equation(params)
%UPDATE_SPIN_DELTAV_FROM_MASTER_EQUATION Compute DeltaV0 from NV dynamics.
%
% This function replaces the phenomenological amplitudes
% params.spin.deltaV_A_meV and params.spin.deltaV_B_meV by values obtained
% from a spin-charge master equation:
%
%   m_s -> Delta P_NVminus(m_s) -> Delta sigma -> Delta V0.
%
% For a 2D NV sheet:
%   Delta V0 = s * e^2 n_NV DeltaP_NVminus / C_A,
%
% where Delta V0 is the electron potential-energy amplitude in meV, n_NV is
% the areal NV density, C_A is the effective interface capacitance per area,
% and s is a sign convention.  The default reference is the midpoint between
% states A and B, so the "none" device remains the unperturbed baseline.

if ~isfield(params.spin, "masterEquation") || ~params.spin.masterEquation.enable
    master.enabled = false;
    return;
end

nv = params.spin.masterEquation;
c = params.constants;

stateA_ms = params.spin.stateA_ms;
stateB_ms = params.spin.stateB_ms;

master.stateA = solve_nv_charge_master(stateA_ms, nv);
master.stateB = solve_nv_charge_master(stateB_ms, nv);
master.enabled = true;

PA = master.stateA.P_NVminus;
PB = master.stateB.P_NVminus;

switch lower(string(nv.referenceMode))
    case "stateb"
        Pref = PB;
    case "statea"
        Pref = PA;
    case "zero"
        Pref = 0;
    otherwise
        Pref = 0.5 * (PA + PB);
end

n_m2 = params.spin.surfaceDensity_cm2 * 1e4;
C_A = local_interface_capacitance_F_m2(nv, c);
energyScale_J = c.q_C^2 * n_m2 / C_A;
energyScale_meV = energyScale_J / c.q_C * 1e3;

signConvention = nv.energySign;
deltaV_A_meV = signConvention * energyScale_meV * (PA - Pref);
deltaV_B_meV = signConvention * energyScale_meV * (PB - Pref);

params.spin.deltaV_A_meV = deltaV_A_meV;
params.spin.deltaV_B_meV = deltaV_B_meV;
params.spin.calibrationMode = "spin-charge master equation";

master.P_NVminus_A = PA;
master.P_NVminus_B = PB;
master.deltaP_AminusB = PA - PB;
master.referenceP_NVminus = Pref;
master.surfaceDensity_cm2 = params.spin.surfaceDensity_cm2;
master.capacitance_F_m2 = C_A;
master.energyScale_meV_per_unit_DeltaP = energyScale_meV;
master.deltaV_A_meV = deltaV_A_meV;
master.deltaV_B_meV = deltaV_B_meV;
master.deltaV_AminusB_meV = deltaV_A_meV - deltaV_B_meV;

params.spin.masterEquationResult = master;

if isfield(nv, "verbose") && nv.verbose
    fprintf("\nNV spin-charge master equation\n");
    fprintf("  mode: %s, readout pulse: %.3g ns\n", string(nv.mode), nv.readoutPulse_s * 1e9);
    fprintf("  P_NVminus(A, m_s=%g) = %.6f\n", stateA_ms, PA);
    fprintf("  P_NVminus(B, m_s=%g) = %.6f\n", stateB_ms, PB);
    fprintf("  DeltaP_AminusB       = %.6g\n", master.deltaP_AminusB);
    fprintf("  C_A                  = %.6g F/m^2\n", C_A);
    fprintf("  scale                = %.6g meV per unit DeltaP\n", energyScale_meV);
    fprintf("  deltaV_A / deltaV_B  = %.6g / %.6g meV\n\n", deltaV_A_meV, deltaV_B_meV);
end
end

function C_A = local_interface_capacitance_F_m2(nv, c)
mode = "parallel_plate";
if isfield(nv, "capacitanceMode")
    mode = lower(string(nv.capacitanceMode));
end

switch mode
    case "manual"
        C_A = nv.capacitance_F_m2;
    otherwise
        eps0 = 8.8541878128e-12;
        d_m = nv.interfaceDielectricThickness_nm * 1e-9;
        C_A = eps0 * nv.interfaceEpsilon_r / d_m;
end

if ~isfinite(C_A) || C_A <= 0
    error("Invalid interface capacitance C_A. Check params.spin.masterEquation.");
end
end
