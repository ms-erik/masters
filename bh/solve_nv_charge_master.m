function out = solve_nv_charge_master(ms, nv)
%SOLVE_NV_CHARGE_MASTER Solve the NV spin-charge master equation.
%
% The default is a finite optical readout pulse, because spin-to-charge
% contrast generally depends on the prepared initial spin state.  A pure
% steady-state solution is also provided for controls and calibration.

[K, labels] = nv_charge_master_matrix(nv);
p0 = nv_initial_population(ms, nv);

mode = lower(string(nv.mode));
switch mode
    case "steady_state"
        A = [K; ones(1, numel(labels))];
        b = [zeros(numel(labels), 1); 1];
        p = A \ b;
    otherwise
        t = nv.readoutPulse_s;
        p = expm(K * t) * p0;
end

p = real(p);
p(p < 0 & p > -1e-12) = 0;
p = p ./ sum(p);

out.ms = ms;
out.labels = labels;
out.K = K;
out.p0 = p0;
out.p = p;
out.P_NVminus = sum(p(1:5));
out.P_NV0 = p(6);
out.P_G0 = p(1);
out.P_G1 = p(2);
out.P_excited = p(3) + p(4);
out.P_singlet = p(5);
out.mode = mode;
out.readoutPulse_s = nv.readoutPulse_s;
end
