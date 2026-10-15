function [chargeClean,dischargeClean,curtailClean,report] = normalize_storage_lp_dispatch(charge,discharge,curtail,eta_ch,eta_dis)
%NORMALIZE_STORAGE_LP_DISPATCH Remove simultaneous operation, keeping SOC/export.
% Intended for solutions that already pass the numerical feasibility check.
% Free, unpenalized curtailment and 0 < efficiencies <= 1 are essential.
charge = double(charge(:));
discharge = double(discharge(:));
curtail = double(curtail(:));
if numel(charge)~=numel(discharge)||numel(charge)~=numel(curtail)
    error('Storage LP cleanup vectors must have equal length.');
end
if any(~isfinite([charge;discharge;curtail]))|| ...
        ~isscalar(eta_ch)||~isscalar(eta_dis)|| ...
        ~isfinite(eta_ch)||~isfinite(eta_dis)|| ...
        eta_ch<=0||eta_ch>1||eta_dis<=0||eta_dis>1
    error('Invalid storage LP cleanup inputs.');
end
alpha = eta_ch*eta_dis;
chargeNonnegative = max(charge,0);
dischargeNonnegative = max(discharge,0);
curtailNonnegative = max(curtail,0);
removedCharge = min(chargeNonnegative,dischargeNonnegative/alpha);
chargeClean = max(chargeNonnegative-removedCharge,0);
dischargeClean = max(dischargeNonnegative-alpha*removedCharge,0);
curtailClean = curtailNonnegative+(1-alpha)*removedCharge;
report = struct();
report.hours_with_simultaneous_operation = nnz(chargeNonnegative>0&dischargeNonnegative>0);
report.removed_charge_MWh = sum(removedCharge);
report.removed_discharge_MWh = sum(alpha*removedCharge);
report.additional_curtailment_MWh = sum(curtailClean-curtail);
report.max_soc_update_change_MWh = max([0;abs( ...
    eta_ch*chargeClean-dischargeClean/eta_dis-(eta_ch*charge-discharge/eta_dis))]);
report.max_power_balance_change_MW = max([0;abs( ...
    dischargeClean-chargeClean-curtailClean-(discharge-charge-curtail))]);
report.max_remaining_simultaneous_power_MW = max([0;min(chargeClean,dischargeClean)]);
end
