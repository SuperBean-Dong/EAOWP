function report = p44_verify_storage_dispatch(Pwind,epsilon,sp,Emax,Eini,energy,charge,discharge,delivery,curtail,requireExclusive)
%VERIFY_STORAGE_LP_DISPATCH Check every original physical constraint numerically.
% Slightly above the default solver feasibility tolerance to accommodate
% roundoff/row scaling. This is a residual acceptance test, not a model change.
Pwind = double(Pwind(:));
energy = double(energy(:));
charge = double(charge(:));
discharge = double(discharge(:));
delivery = double(delivery(:));
curtail = double(curtail(:));
T = numel(Pwind);
report = struct('feasible',false,'max_violation',Inf,'tolerance',NaN, ...
    'dimension_check',false,'finite_check',false);
report.dimension_check = isscalar(Emax)&&isscalar(Eini)&& ...
    numel(energy)==T+1&&numel(charge)==T&&numel(discharge)==T&& ...
    numel(delivery)==T&&numel(curtail)==T&&T>0;
if ~report.dimension_check; return; end
report.finite_check = all(isfinite([Pwind;epsilon;Emax;Eini;energy;charge;discharge;delivery;curtail]));
if ~report.finite_check; return; end
scale = max([1;abs(Pwind);abs(epsilon);abs(Emax);abs(Eini);abs(energy)]);
report.tolerance = 1e-5+1e-8*scale;
M = max(max(Pwind),epsilon);
report.capacity_nonnegativity = max(0,-Emax);
report.flow_nonnegativity = max([0;-charge;-discharge;-delivery;-curtail]);
report.charge_power_capacity = max([0;charge-Emax/sp.duration_h]);
report.discharge_power_capacity = max([0;discharge-Emax/sp.duration_h]);
report.charge_M_bound = max([0;charge-M]);
report.discharge_M_bound = max([0;discharge-M]);
report.export_bound = max([0;delivery-epsilon]);
report.initial_soc_lower = max(0,sp.soc_min*Emax-Eini);
report.initial_soc_upper = max(0,Eini-sp.soc_max*Emax);
report.initial_state_equality = abs(energy(1)-Eini);
report.power_balance = max([0;abs(Pwind+discharge-charge-delivery-curtail)]);
report.energy_recurrence = max([0;abs(energy(2:end)- ...
    (energy(1:end-1)*(1-sp.self_discharge_per_h)+sp.eta_ch*charge-discharge/sp.eta_dis))]);
report.soc_lower = max([0;sp.soc_min*Emax-energy]);
report.soc_upper = max([0;energy-sp.soc_max*Emax]);
report.terminal_energy = max(0,Eini-energy(end));
if requireExclusive
    report.simultaneous_charge_discharge = max([0;min(charge,discharge)]);
else
    report.simultaneous_charge_discharge = 0;
end
report.max_violation = max([report.capacity_nonnegativity,report.flow_nonnegativity, ...
    report.charge_power_capacity,report.discharge_power_capacity,report.charge_M_bound, ...
    report.discharge_M_bound,report.export_bound,report.initial_soc_lower, ...
    report.initial_soc_upper,report.initial_state_equality,report.power_balance, ...
    report.energy_recurrence,report.soc_lower,report.soc_upper, ...
    report.terminal_energy,report.simultaneous_charge_discharge]);
report.feasible = report.max_violation<=report.tolerance;
end
