function S = b0_compact_candidates(B)
% Rows: storage-HVAC, storage-HVDC, no-storage-HVAC, no-storage-HVDC.
n=numel(B);
S=struct();
S.row_labels={'storage_HVAC','storage_VSC_HVDC','no_storage_HVAC','no_storage_VSC_HVDC'};
S.epsilon_MW=nan(1,n);
S.stage2_npv_CNY=nan(4,n);
S.export_system_cost_CNY=nan(4,n);
S.export_system_capex_CNY=nan(4,n);
S.export_maintenance_pv_CNY=nan(4,n);
S.battery_cost_CNY=nan(4,n);
S.annual_revenue_CNY=nan(4,n);
S.storage_energy_MWh=nan(2,n);
S.feasible=false(4,n);
for k=1:n
    S.epsilon_MW(k)=B(k).epsilon_MW;
    A=B(k).storage_optimized; N=B(k).no_storage;
    if isfield(A,'HVAC') && ~isempty(fieldnames(A.HVAC))
        S.stage2_npv_CNY(1,k)=A.HVAC.stage2_npv_CNY;
        S.export_system_cost_CNY(1,k)=A.HVAC.export_system_cost_CNY;
        S.export_system_capex_CNY(1,k)=A.HVAC.export_system_capex_CNY;
        S.export_maintenance_pv_CNY(1,k)=A.HVAC.export_maintenance_pv_CNY;
        S.battery_cost_CNY(1,k)=A.HVAC.battery_cost_CNY;
        S.annual_revenue_CNY(1,k)=A.HVAC.annual_revenue_CNY;
        S.feasible(1,k)=A.HVAC.feasible;
        S.storage_energy_MWh(1,k)=A.HVAC.storage_energy_MWh;
    end
    if isfield(A,'VSC_HVDC') && ~isempty(fieldnames(A.VSC_HVDC))
        S.stage2_npv_CNY(2,k)=A.VSC_HVDC.stage2_npv_CNY;
        S.export_system_cost_CNY(2,k)=A.VSC_HVDC.export_system_cost_CNY;
        S.export_system_capex_CNY(2,k)=A.VSC_HVDC.export_system_capex_CNY;
        S.export_maintenance_pv_CNY(2,k)=A.VSC_HVDC.export_maintenance_pv_CNY;
        S.battery_cost_CNY(2,k)=A.VSC_HVDC.battery_cost_CNY;
        S.annual_revenue_CNY(2,k)=A.VSC_HVDC.annual_revenue_CNY;
        S.feasible(2,k)=A.VSC_HVDC.feasible;
    end
    if isfield(N,'HVAC') && ~isempty(fieldnames(N.HVAC))
        S.stage2_npv_CNY(3,k)=N.HVAC.stage2_npv_CNY;
        S.export_system_cost_CNY(3,k)=N.HVAC.export_system_cost_CNY;
        S.export_system_capex_CNY(3,k)=N.HVAC.export_system_capex_CNY;
        S.export_maintenance_pv_CNY(3,k)=N.HVAC.export_maintenance_pv_CNY;
        S.battery_cost_CNY(3,k)=N.HVAC.battery_cost_CNY;
        S.annual_revenue_CNY(3,k)=N.HVAC.annual_revenue_CNY;
        S.feasible(3,k)=N.HVAC.feasible;
        S.storage_energy_MWh(2,k)=0;
    end
    if isfield(N,'VSC_HVDC') && ~isempty(fieldnames(N.VSC_HVDC))
        S.stage2_npv_CNY(4,k)=N.VSC_HVDC.stage2_npv_CNY;
        S.export_system_cost_CNY(4,k)=N.VSC_HVDC.export_system_cost_CNY;
        S.export_system_capex_CNY(4,k)=N.VSC_HVDC.export_system_capex_CNY;
        S.export_maintenance_pv_CNY(4,k)=N.VSC_HVDC.export_maintenance_pv_CNY;
        S.battery_cost_CNY(4,k)=N.VSC_HVDC.battery_cost_CNY;
        S.annual_revenue_CNY(4,k)=N.VSC_HVDC.annual_revenue_CNY;
        S.feasible(4,k)=N.VSC_HVDC.feasible;
    end
end
end
