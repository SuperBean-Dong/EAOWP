function q=p44_revalue_grid(r,keep,cfg)
% Reuse the literal B0 financial values for unchanged layouts.
keep=logical(keep(:)); assert(numel(keep)==r.turbine_count && any(keep));
q=baselineSummary(r,cfg);
if all(keep); return; end
key=p44_sha256([cfg.ScienceSignature '|' num2str(r.grid_index) '|' char('0'+keep.')]);
file=fullfile(cfg.CacheDir,sprintf('grid_%04d_%s.mat',r.grid_index,key));
if isfile(file)
    C=load(file,'summary','cache_signature');
    assert(strcmp(C.cache_signature,key),'Invalid repaired-grid cache.');
    q=C.summary; q.cache_hit=true; return;
end
timer=tic;
wind=p44_load_wind(r.grid_index,cfg);
generation=p44_recalculate_wind(r,keep,wind,cfg);
collection=p44_repair_collection(r,keep,cfg);
farm=p44_farm_cost(r,keep,collection.junction_original_indices,cfg);
path=fullfile(cfg.BaselineDir,'stage4_checkpoints',sprintf('grid_%04d.mat',r.grid_index));
C=load(path,'checkpoint'); selected=C.checkpoint.selected_candidate;
assert(C.checkpoint.grid_index==r.grid_index && C.checkpoint.search_complete);
assert(strcmp(C.checkpoint.model_signature,cfg.BaselineModelSignature));
dispatch=p44_fixed_dispatch(generation.Pwind_MW,r,selected,cfg.PriceCNYPerMWh,cfg);
financial=p44_fixed_financial(r,selected,dispatch,farm.platform_turbine_cost_CNY, ...
    collection.internal_cable_cost_CNY,cfg);
q.kept_original_indices=find(keep); q.turbine_count=sum(keep);
q.deleted_turbine_count=r.turbine_count-sum(keep); q.capacity_MW=generation.capacity_MW;
q.NPV_CNY=financial.NPV_CNY; q.LCOE_CNY_per_kWh=financial.LCOE_CNY_per_kWh;
q.AEP_wake_adjusted_MWh=generation.AEP_wake_adjusted_MWh;
q.annual_onshore_delivery_MWh=financial.annual_onshore_delivery_MWh;
q.farm_cost_CNY=farm.platform_turbine_cost_CNY;
q.collection_cost_CNY=collection.internal_cable_cost_CNY;
q.total_cost_CNY=financial.total_cost_CNY;
q.junction_count=collection.junction_count; q.junction_proxy_cost_CNY=farm.junction_proxy_cost_CNY;
q.baseline_reused=false; q.cache_file=file; q.cache_hit=false; q.financial_evaluated=true;
q.recompute_seconds=toc(timer); q.dispatch_status=dispatch.status;
summary=q; cache_signature=key;
temporary=[file '.tmp.mat'];
save(temporary,'summary','cache_signature','generation','collection','farm','dispatch','financial','-v7.3');
movefile(temporary,file,'f');
end
function q=baselineSummary(r,cfg)
y=(1:cfg.LifeFarmYears).'; factor=sum((1-cfg.AnnualWindDegradationRate).^(y-1)./(1+cfg.AnnualInterestRate).^y);
q=struct('grid_id',double(r.grid_index),'original_capacity_MW',double(r.wind_farm_capacity), ...
    'capacity_MW',double(r.wind_farm_capacity),'original_turbine_count',double(r.turbine_count), ...
    'turbine_count',double(r.turbine_count),'deleted_turbine_count',0, ...
    'kept_original_indices',(1:double(r.turbine_count)).', ...
    'original_NPV_CNY',double(r.NPV_best),'NPV_CNY',double(r.NPV_best), ...
    'LCOE_CNY_per_kWh',double(r.LCOE_min),'AEP_wake_adjusted_MWh',double(r.AEP_wake_adjusted_MWh), ...
    'annual_onshore_delivery_MWh',double(r.discounted_onshore_energy_kWh)/factor/1000, ...
    'farm_cost_CNY',double(r.platform_turbine_cost),'collection_cost_CNY',double(r.internal_cable_cost), ...
    'export_system_cost_CNY',double(r.export_system_cost_CNY),'battery_cost_CNY',double(r.battery_cost_CNY), ...
    'storage_energy_MWh',double(r.E_bat_best),'export_capacity_MW',double(r.optimal_epsilon), ...
    'total_cost_CNY',double(r.total_cost_best),'junction_count',0,'junction_proxy_cost_CNY',0, ...
    'baseline_reused',true,'cache_file','','cache_hit',false,'financial_evaluated',true, ...
    'recompute_seconds',0,'dispatch_status','reuse_literal_B0_values');
end
