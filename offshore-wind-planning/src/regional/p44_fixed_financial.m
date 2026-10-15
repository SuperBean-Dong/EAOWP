function financial = p44_fixed_financial(r,selectedCandidate,dispatch,updatedFarmCost,updatedCollectionCost,cfg)
% Preserve B0 export/BESS costs and cable choice; recalculate shore delivery.
assert(dispatch.feasible&&dispatch.is_optimal,'A feasible, optimal/reused dispatch is required.');
assert(isscalar(updatedFarmCost)&&isfinite(updatedFarmCost)&&updatedFarmCost>=0);
assert(isscalar(updatedCollectionCost)&&isfinite(updatedCollectionCost)&&updatedCollectionCost>=0);
P=double(dispatch.farm_delivery_MW(:)); distance=double(r.distance_to_shore_km);
td=selectedCandidate.transmission_detail;
index=double(td.selected_cable_index); n=double(td.n_circuits);
assert(isfinite(distance)&&distance>=0&&n>=1&&n==floor(n));
mode=char(selectedCandidate.export_mode);
assert(strcmp(mode,char(r.export_mode)),'Fixed export mode mismatch.');
Pactual=P*1e6;
if strcmp(mode,'HVAC')
    resistance=[48.9e-3,39.1e-3,31.9e-3,27e-3]; V=220e3;
    assert(index>=1&&index<=4&&index==floor(index));
    P1=Pactual*(1-.0004);
    loss=3*(P1/(n*sqrt(3)*V)).^2*resistance(index)*distance*n;
    onshore=max(P1-loss,0)/1e3;
elseif strcmp(mode,'VSC-HVDC')
    resistance=[22.4e-3,19.2e-3,16.5e-3,14.4e-3,11.5e-3]; V=300e3;
    assert(index>=1&&index<=5&&index==floor(index));
    P1=Pactual*(1-.0143);
    loss=2*(P1/(n*V)).^2*resistance(index)*distance*n;
    onshore=max(P1-loss,0)*(1-.0143)/1e3;
else
    error('p44:UnknownExportMode','Unknown export mode: %s',mode);
end
price=double(cfg.PriceCNYPerMWh(:));
if isscalar(price); price=repmat(price,numel(P),1); end
assert(numel(price)==numel(P)&&all(isfinite(price))&&all(price>0),'Invalid financial price.');
annualRevenue=sum(onshore/1000.*price);
annualEnergy=sum(onshore);
% Exact saved values preserve B0's reporting on the no-change path.
if dispatch.reused_baseline
    original=double(selectedCandidate.onshore_power_kW(:));
    assert(numel(original)==numel(onshore)&&max(abs(original-onshore))<=max(1e-6,1e-10*max([1;abs(original)])), ...
        'Fixed loss formula failed to reproduce saved B0 shore power.');
    onshore=original; annualEnergy=selectedCandidate.annual_delivered_energy_kWh;
    assert(abs(annualRevenue-selectedCandidate.annual_revenue_CNY)<=max(1,1e-10*abs(annualRevenue)), ...
        'Saved B0 shore revenue differs from the supplied price.');
    annualRevenue=selectedCandidate.annual_revenue_CNY;
end
years=(1:cfg.LifeFarmYears).';
factor=sum((1-cfg.AnnualWindDegradationRate).^(years-1)./(1+cfg.AnnualInterestRate).^years);
exportCost=double(selectedCandidate.export_system_cost_CNY);
batteryCost=double(selectedCandidate.battery_cost_CNY);
assert(isfinite(exportCost)&&exportCost>=0&&isfinite(batteryCost)&&batteryCost>=0);
if isfield(r,'export_system_cost_CNY')
    assertClose(double(r.export_system_cost_CNY),exportCost,'Fixed export lifecycle cost');
end
if isfield(r,'battery_cost_CNY')
    assertClose(double(r.battery_cost_CNY),batteryCost,'Fixed BESS lifecycle cost');
end
if all(isfield(r,{'platform_turbine_cost','internal_cable_cost','total_cost_best','NPV_best'}))
    originalCost=double(r.platform_turbine_cost)+double(r.internal_cable_cost)+exportCost+batteryCost;
    assertClose(double(r.total_cost_best),originalCost,'Original B0 total lifecycle cost');
    assertClose(double(r.NPV_best),double(selectedCandidate.pv_revenue_CNY)-originalCost, ...
        'Original B0 NPV and selected snapshot');
end
total=updatedFarmCost+updatedCollectionCost+exportCost+batteryCost;
pvRevenue=annualRevenue*factor;
discountedEnergy=annualEnergy*factor;
financial=struct('farm_cost_CNY',updatedFarmCost,'collection_cost_CNY',updatedCollectionCost, ...
    'export_system_cost_CNY',exportCost,'battery_cost_CNY',batteryCost, ...
    'total_cost_best',total,'NPV_best',pvRevenue-total,'LCOE_min',total/max(discountedEnergy,eps), ...
    'annual_revenue_CNY',annualRevenue,'pv_revenue_CNY',pvRevenue, ...
    'annual_delivered_energy_kWh',annualEnergy,'discounted_onshore_energy_kWh',discountedEnergy, ...
    'AEP_wake_adjusted_MWh',sum(dispatch.Pwind_MW), ...
    'annual_available_energy_MWh',sum(dispatch.Pwind_MW), ...
    'annual_curtailed_energy_MWh',sum(dispatch.curtailment_MW), ...
    'annual_farm_delivery_MWh',sum(P),'annual_onshore_delivery_MWh',annualEnergy/1000, ...
    'curtailment_fraction',sum(dispatch.curtailment_MW)/max(sum(dispatch.Pwind_MW),eps), ...
    'onshore_power_kW',onshore,'export_mode',mode,'selected_cable_index',index, ...
    'n_circuits',n,'storage_energy_MWh',dispatch.storage_energy_MWh, ...
    'epsilon_MW',dispatch.epsilon_MW,'fixed_equipment_costs_preserved',true, ...
    'dispatch_objective_basis',dispatch.objective_basis);
financial.NPV_CNY=financial.NPV_best;
financial.LCOE_CNY_per_kWh=financial.LCOE_min;
financial.total_cost_CNY=total;
end

function assertClose(saved,recomputed,label)
assert(isfinite(saved)&&isfinite(recomputed)&& ...
    abs(saved-recomputed)<=max(1,1e-10*max(abs([saved,recomputed]))), ...
    'p44:BaselineSnapshotMismatch','%s differs between the input record and selected B0 snapshot.',label);
end
