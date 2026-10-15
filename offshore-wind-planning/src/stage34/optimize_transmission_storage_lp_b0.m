function out = optimize_transmission_storage_lp_b0(P_rated,distance,storage_par,power_pu, ...
    life_farm,annual_interest_rate,power_price,varargin)
%OPTIMIZE_TRANSMISSION_STORAGE_LP_FAST_TEST Standalone sensitivity LP search.
% Production source SHA256:
% 6d0cf906871e91efa9c63cfa338b5cb0d82114a9e7e92bdf7979218240fd4040
% MATLAB + direct Gurobi. The dispatch LP is algebraically identical to the
% production YALMIP model. Its sparse constraint matrix is cached by horizon
% and physical battery parameters; only RHS, upper bounds and costs change.
% E_ini is a FREE decision variable, as in production, not fixed to 0.5 E_max.
% The production LP already removes hourly charge/discharge binaries. Since free
% curtailment is available, simultaneous operation can be eliminated after
% solving without changing SOC, export, capacity, or the storage objective.
% Search deadlines are soft wall-time checkpoints: model compilation and
% solver termination may slightly exceed the requested remaining time.
%
% Final-formula implementation highlights:
%   1) Export-capacity candidates keep the code-side revision:
%      regular 10-MW grid + exact hourly wind peak + installed capacity.
%   2) Every epsilon retains BOTH storage-optimized and explicit E_bat=0 branches.
%   3) Eq. (41): BESS purchase CAPEX is converted to an equivalent annual cost
%      over the battery service life, annual fixed O&M is added, and that annual
%      cost is incurred in every project year and discounted to present value.
%   4) BESS platform CAPEX is an additional storage CAPEX item.
%   5) HVAC/VSC-HVDC annual maintenance is included in lifecycle cost.
%   6) Lifetime delivered-shore revenue and discounted energy include annual
%      wind-turbine generation degradation.
%
% storage_par required fields:
%   eta_ch, eta_dis, self_discharge_per_h, soc_max, soc_min, duration_h,
%   energy_capex_CNY_per_kWh, platform_capex_EUR_per_kWh,
%   battery_life_years, om_fraction_per_year, EUR2CNY

p = inputParser;
p.addParameter('ExportStepMW',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinimumExportCandidateMW',50,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('StorageVerbose',0,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('AnnualWindDegradationRate',0.0072,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
p.addParameter('EUR2CNY',8.39,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.addParameter('StorageTimeLimitSeconds',30,@(x)isnumeric(x)&&isscalar(x)&&~isnan(x)&&x>0);
p.addParameter('StorageGridTimeLimitSeconds',300,@(x)isnumeric(x)&&isscalar(x)&&~isnan(x)&&x>=0);
p.addParameter('SolverThreads',2,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==floor(x));
p.addParameter('ExportCapacityCandidatesMW',[],@(x)isnumeric(x)&&isvector(x)||isempty(x));
p.addParameter('ReuseLPBasis',false,@(x)islogical(x)&&isscalar(x)||isnumeric(x)&&isscalar(x)&&ismember(x,[0 1]));
p.addParameter('CertifiedStoragePruning',[],@(x)isempty(x)||isstruct(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

storage_par = validateStorageParameters(storage_par,opt.EUR2CNY);

power_pu = double(power_pu(:));
power_pu(~isfinite(power_pu)) = 0;
power_pu = min(max(power_pu,0),1);
T = numel(power_pu);
price = normalizePrice(power_price,T);
P_wind = P_rated.*power_pu; % MW, first-year hourly available farm power
P_peak = max(P_wind);

% Regular search grid + exact actual peak + exact installed capacity.
regular = opt.MinimumExportCandidateMW:opt.ExportStepMW:P_rated;
if P_rated < opt.MinimumExportCandidateMW
    regular = [];
end
epsilon_set = [regular(:);P_peak;P_rated];
epsilon_set = epsilon_set(isfinite(epsilon_set)&epsilon_set>0);
epsilon_set = unique(epsilon_set,'sorted').';
% Optional explicit subset is for equivalence checks only. Normal testing
% leaves this empty and uses every original production epsilon candidate.
if ~isempty(opt.ExportCapacityCandidatesMW)
    epsilon_set = unique(double(opt.ExportCapacityCandidatesMW(:)),'sorted').';
    if any(~isfinite(epsilon_set)|epsilon_set<=0)
        error('ExportCapacityCandidatesMW must contain finite positive values.');
    end
end
if isempty(epsilon_set)
    epsilon_set = max(P_rated,eps);
end
[certifiedPruning,certifiedMask] = validateCertifiedPruning(opt.CertifiedStoragePruning,epsilon_set);

candidateTemplate = struct( ...
    'epsilon_MW',NaN, ...
    'storage_case','', ...
    'export_mode','', ...
    'storage_energy_MWh',NaN, ...
    'battery_initial_capex_CNY',NaN, ...
    'battery_annualized_capex_CNY_per_year',NaN, ...
    'battery_annual_OM_CNY_per_year',NaN, ...
    'battery_cost_CNY',NaN, ...
    'farm_delivery_MW',[], ...
    'curtailment_MW',[], ...
    'annual_delivered_energy_kWh',NaN, ...
    'annual_revenue_CNY',NaN, ...
    'pv_revenue_CNY',NaN, ...
    'export_system_capex_CNY',NaN, ...
    'export_annual_maintenance_CNY',NaN, ...
    'export_maintenance_pv_CNY',NaN, ...
    'export_system_cost_CNY',NaN, ...
    'stage2_npv_CNY',-Inf, ...
    'onshore_power_kW',[], ...
    'transmission_detail',struct(), ...
    'storage_detail',struct(), ...
    'storage_evaluation_status','', ...
    'storage_dispatch_is_optimal',false, ...
    'feasible',false);
candidates = repmat(candidateTemplate,0,1);
by_epsilon = repmat(struct('epsilon_MW',NaN,'storage_optimized',struct(), ...
    'no_storage',struct()),numel(epsilon_set),1);
storage_status_by_epsilon = repmat(storageStatusTemplate(),numel(epsilon_set),1);
storageSearchTimer = tic;

for ee = 1:numel(epsilon_set)
    epsilon = epsilon_set(ee);
    by_epsilon(ee).epsilon_MW = epsilon;

    % Branch A: storage optimum under this epsilon using final Eq. (41).
    remainingGridSeconds = opt.StorageGridTimeLimitSeconds-toc(storageSearchTimer);
    if certifiedMask(ee)
        Eopt = NaN;
        PdelOpt = nan(T,1);
        storageOpt = emptyStorageDetail(T,'CERTIFIED_BELOW_NO_STORAGE_INCUMBENT');
        storageOpt.solver_status = 'CERTIFIED_BELOW_NO_STORAGE_INCUMBENT';
        storageOpt.solver_info = 'Certified shore-NPV upper bound is strictly below an explicit no-storage incumbent.';
        storageOpt.certified_pruned = true;
        storageOpt.certified_upper_bound_CNY = certifiedPruning.upper_bound_CNY(ee);
        storageOpt.rounding_tolerance_CNY = certifiedPruning.rounding_tolerance_CNY(ee);
        storageOpt.no_storage_incumbent_CNY = certifiedPruning.no_storage_incumbent_CNY;
    elseif remainingGridSeconds<=0
        Eopt = NaN;
        PdelOpt = nan(T,1);
        storageOpt = emptyStorageDetail(T,'skipped_grid_time_budget');
        storageOpt.solver_status = 'NOT_RUN_GRID_TIME_BUDGET';
        storageOpt.solver_info = 'Storage grid wall-time budget exhausted; no-storage branch retained.';
    else
        [Eopt,PdelOpt,storageOpt] = optimizeStorageDispatch(P_wind,epsilon,storage_par, ...
            life_farm,price,annual_interest_rate,opt.AnnualWindDegradationRate,opt.StorageVerbose, ...
            opt.StorageTimeLimitSeconds,remainingGridSeconds,opt.SolverThreads,opt.ReuseLPBasis);
    end
    storage_status_by_epsilon(ee) = summarizeStorageStatus(epsilon,storageOpt);
    branchOpt = evaluateBranch('storage_optimized',Eopt,PdelOpt,storageOpt, ...
        epsilon,distance,storage_par,life_farm,annual_interest_rate,price, ...
        opt.AnnualWindDegradationRate,opt.EUR2CNY);
    by_epsilon(ee).storage_optimized = branchOpt;
    candidates = appendModeCandidates(candidates,branchOpt);

    % Branch B: explicit E_bat = 0 candidate.
    Pdel0 = min(P_wind,epsilon);
    storage0 = makeNoStorageDetail(P_wind,Pdel0);
    branch0 = evaluateBranch('no_storage',0,Pdel0,storage0,epsilon,distance, ...
        storage_par,life_farm,annual_interest_rate,price, ...
        opt.AnnualWindDegradationRate,opt.EUR2CNY);
    by_epsilon(ee).no_storage = branch0;
    candidates = appendModeCandidates(candidates,branch0);
    if ~isempty(certifiedPruning)
        localNoStorage = [branch0.HVAC.stage2_npv_CNY,branch0.VSC_HVDC.stage2_npv_CNY];
        localNoStorage = localNoStorage(isfinite(localNoStorage));
        if ~isempty(localNoStorage) && ...
                certifiedPruning.upper_bound_CNY(ee)+certifiedPruning.rounding_tolerance_CNY(ee)<max(localNoStorage)-1
            error('Certified upper bound is below this epsilon''s explicit feasible no-storage NPV.');
        end
    end
end

if isempty(candidates)
    error('No transmission-storage candidates were generated.');
end
npv = [candidates.stage2_npv_CNY];
feas = [candidates.feasible]&isfinite(npv);
if ~any(feas)
    error('All HVAC/VSC-HVDC transmission-storage candidates are infeasible.');
end
validIdx = find(feas);
[~,k] = max(npv(validIdx));
bestIdx = validIdx(k);
best = candidates(bestIdx);
certifiedValidation=struct('enabled',false);
if ~isempty(certifiedPruning)
    noStorageMask=strcmp({candidates.storage_case},'no_storage')&feas;
    actualNoStorageIncumbent=max(npv(noStorageMask));
    incumbentTolerance=max(1,1e-10*abs(certifiedPruning.no_storage_incumbent_CNY));
    incumbentDifference=abs(actualNoStorageIncumbent-certifiedPruning.no_storage_incumbent_CNY);
    assert(incumbentDifference<=incumbentTolerance, ...
        'Frozen certified-pruning incumbent disagrees with the actual complete no-storage search.');
    certifiedValidation=struct('enabled',true,'pruned_count',nnz(certifiedMask), ...
        'provided_no_storage_incumbent_CNY',certifiedPruning.no_storage_incumbent_CNY, ...
        'actual_no_storage_incumbent_CNY',actualNoStorageIncumbent, ...
        'incumbent_difference_CNY',incumbentDifference,'incumbent_tolerance_CNY',incumbentTolerance, ...
        'all_upper_bounds_above_corresponding_no_storage_candidate',true);
end

out = struct();
out.actual_wind_peak_MW = P_peak;
out.installed_capacity_MW = P_rated;
out.annual_wind_degradation_rate = opt.AnnualWindDegradationRate;
out.epsilon_set_MW = epsilon_set;
out.storage_parameters = storage_par;
out.by_epsilon = by_epsilon;
out.candidates = candidates;
out.best_candidate_index = bestIdx;
out.best = best;
out.best_storage_energy_MWh = best.storage_energy_MWh;
out.best_epsilon_MW = best.epsilon_MW;
out.best_export_mode = best.export_mode;
out.best_storage_case = best.storage_case;
out.best_stage2_npv_CNY = best.stage2_npv_CNY;
out.best_onshore_power_kW = best.onshore_power_kW;
out.best_transmission_detail = best.transmission_detail;
out.best_farm_delivery_MW = best.farm_delivery_MW;
out.best_battery_cost_CNY = best.battery_cost_CNY;
out.best_export_system_cost_CNY = best.export_system_cost_CNY;
out.storage_status_by_epsilon = storage_status_by_epsilon;
out.storage_all_epsilon_attempted = all([storage_status_by_epsilon.attempted]);
out.storage_search_all_optimal = all([storage_status_by_epsilon.is_optimal]);
out.storage_search_complete = all([storage_status_by_epsilon.is_optimal]| ...
    [storage_status_by_epsilon.certified_pruned]);
out.certified_storage_pruning_count = nnz([storage_status_by_epsilon.certified_pruned]);
out.certified_storage_pruning_validation = certifiedValidation;
out.storage_search_wall_time_s = toc(storageSearchTimer);
out.storage_time_limit_seconds = opt.StorageTimeLimitSeconds;
out.storage_grid_time_limit_seconds = opt.StorageGridTimeLimitSeconds;
out.solver_threads = opt.SolverThreads;
out.reuse_lp_basis = logical(opt.ReuseLPBasis);
out.dispatch_implementation = 'direct_Gurobi_cached_sparse_LP';
out.production_source_sha256 = '6d0cf906871e91efa9c63cfa338b5cb0d82114a9e7e92bdf7979218240fd4040';
if out.storage_search_complete
    if out.certified_storage_pruning_count>0
        out.optimality_status = 'all_storage_dispatches_optimal_or_certified_pruned_on_epsilon_grid';
    else
        out.optimality_status = 'all_storage_dispatches_optimal_on_epsilon_grid';
    end
else
    out.optimality_status = 'best_known_feasible_candidate_incomplete_storage_search';
end
end

function branch = evaluateBranch(caseName,Ebat,Pdelivery,storageDetail,epsilon,distance, ...
    storage_par,life_farm,r,price,degradation,EUR2CNY)
branch = struct();
branch.case_name = caseName;
branch.storage_energy_MWh = Ebat;
branch.farm_delivery_MW = Pdelivery(:);
branch.storage_detail = storageDetail;
branch.HVAC = struct();
branch.VSC_HVDC = struct();

if ~isfinite(Ebat)||isempty(Pdelivery)||any(~isfinite(Pdelivery))
    return;
end

B = batteryLifecycleCost(Ebat,storage_par,life_farm,r);
Pmax = max(Pdelivery);

[~,~,Pac,ac] = calculateACTransmissionCost(Pmax,Pdelivery*1e6,distance,life_farm,r, ...
    price,degradation,EUR2CNY);
branch.HVAC = makeModeResult('HVAC',caseName,epsilon,Ebat,B,Pdelivery, ...
    storageDetail,Pac,ac,life_farm,r,degradation);

[~,~,Pdc,dc] = calculateDCTransmissionCost(Pmax,Pdelivery*1e6,distance,life_farm,r, ...
    price,degradation,EUR2CNY);
branch.VSC_HVDC = makeModeResult('VSC-HVDC',caseName,epsilon,Ebat,B,Pdelivery, ...
    storageDetail,Pdc,dc,life_farm,r,degradation);
end

function c = makeModeResult(modeName,caseName,epsilon,Ebat,B,Pdelivery, ...
    storageDetail,Ponshore,transDetail,life_farm,r,degradation)
c = struct();
c.epsilon_MW = epsilon;
c.storage_case = caseName;
c.export_mode = modeName;
c.storage_energy_MWh = Ebat;
c.battery_initial_capex_CNY = B.initial_capex_CNY;
c.battery_annualized_capex_CNY_per_year = B.annualized_capex_CNY_per_year;
c.battery_annual_OM_CNY_per_year = B.annual_OM_CNY_per_year;
c.battery_cost_CNY = B.lifecycle_PV_CNY;
c.farm_delivery_MW = Pdelivery(:);
if isfield(storageDetail,'curtailment_MW')
    c.curtailment_MW = storageDetail.curtailment_MW(:);
else
    c.curtailment_MW = [];
end
c.onshore_power_kW = Ponshore(:);
c.transmission_detail = transDetail;
c.storage_detail = storageDetail;
c.storage_evaluation_status = storageDetail.status;
c.storage_dispatch_is_optimal = storageDetail.is_optimal;
c.annual_delivered_energy_kWh = sum(c.onshore_power_kW);
if isfield(transDetail,'annual_revenue_CNY')
    c.annual_revenue_CNY = transDetail.annual_revenue_CNY;
else
    c.annual_revenue_CNY = NaN;
end
c.pv_revenue_CNY = presentValueDegraded(c.annual_revenue_CNY,life_farm,r,degradation);

if isfield(transDetail,'capex_10kCNY')
    c.export_system_capex_CNY = transDetail.capex_10kCNY*1e4;
else
    c.export_system_capex_CNY = Inf;
end
if isfield(transDetail,'annual_maintenance_10kCNY')
    c.export_annual_maintenance_CNY = transDetail.annual_maintenance_10kCNY*1e4;
else
    c.export_annual_maintenance_CNY = NaN;
end
if isfield(transDetail,'pv_maintenance_10kCNY')
    c.export_maintenance_pv_CNY = transDetail.pv_maintenance_10kCNY*1e4;
else
    c.export_maintenance_pv_CNY = NaN;
end
if isfield(transDetail,'lifecycle_cost_10kCNY')
    c.export_system_cost_CNY = transDetail.lifecycle_cost_10kCNY*1e4;
elseif isfield(transDetail,'total_cost_10kCNY')
    c.export_system_cost_CNY = transDetail.total_cost_10kCNY*1e4;
else
    c.export_system_cost_CNY = Inf;
end

c.stage2_npv_CNY = c.pv_revenue_CNY - c.export_system_cost_CNY - c.battery_cost_CNY;
c.feasible = isfinite(c.stage2_npv_CNY)&&isfinite(c.export_system_cost_CNY)&&isfinite(c.battery_cost_CNY);
end

function C = appendModeCandidates(C,branch)
if isfield(branch,'HVAC')&&~isempty(fieldnames(branch.HVAC))
    if isempty(C)
        C = branch.HVAC;
    else
        C(end+1,1) = orderfields(branch.HVAC,C(1)); 
    end
end
if isfield(branch,'VSC_HVDC')&&~isempty(fieldnames(branch.VSC_HVDC))
    if isempty(C)
        C = branch.VSC_HVDC;
    else
        C(end+1,1) = orderfields(branch.VSC_HVDC,C(1)); 
    end
end
end

function detail = makeNoStorageDetail(Pwind,Pdelivery)
detail = struct();
detail.feasible = true;
detail.storage_energy_MWh = 0;
detail.charge_MW = zeros(size(Pwind));
detail.discharge_MW = zeros(size(Pwind));
detail.energy_MWh = zeros(numel(Pwind)+1,1);
detail.curtailment_MW = max(Pwind(:)-Pdelivery(:),0);
detail.farm_delivery_MW = Pdelivery(:);
detail.solver_info = 'explicit E_bat=0 candidate';
detail.initial_capex_CNY = 0;
detail.annualized_capex_CNY_per_year = 0;
detail.annual_OM_CNY_per_year = 0;
detail.lifecycle_PV_CNY = 0;
detail.status = 'explicit_no_storage';
detail.solver_status = 'NOT_RUN_NO_STORAGE';
detail.is_optimal = true;
detail.accepted_feasible = true;
detail.attempted = false;
end

function [Emax,Pdelivery,detail] = optimizeStorageDispatch(Pwind,epsilon,sp,life_farm,price,r,degradation,verbose,perSolveSeconds,remainingGridSeconds,solverThreads,reuseLPBasis)
% Final Eq. (41): annualized BESS acquisition + annual fixed O&M - annual
% electricity benefit, repeated over project life and discounted to present.
attemptTimer = tic;
Pwind = double(Pwind(:));
T = numel(Pwind);
[model,ix] = cachedDispatchLP(T,sp);

% C_ann^B = E_max * unit CAPEX * CRF(L_B,r) + E_max * unit CAPEX * o_B.
unitInitialCapex_CNY_per_MWh = combinedStorageCapexPerMWh(sp);
crf = capitalRecoveryFactor(r,sp.battery_life_years);
annualizedUnitCapex = unitInitialCapex_CNY_per_MWh*crf;
annualUnitOM = unitInitialCapex_CNY_per_MWh*sp.om_fraction_per_year;
pvStorageUnitCost = (annualizedUnitCapex + annualUnitOM)*annuityPVFactor(life_farm,r);

% First-year operating benefit. Lifetime benefit includes 0.72%/yr degradation.
pvBenefitFactor = degradedPVFactor(life_farm,r,degradation);
model.obj(ix.Emax) = pvStorageUnitCost;
model.obj(ix.delivery) = -price(:)*pvBenefitFactor;
M = max(max(Pwind),epsilon);
model.ub(ix.charge) = M;
model.ub(ix.discharge) = M;
model.ub(ix.delivery) = epsilon;
model.rhs(ix.balanceRows) = -Pwind;

detail = emptyStorageDetail(T,'solver_not_started');
Emax = NaN;
Pdelivery = nan(T,1);
detail.modeling_time_s = toc(attemptTimer);
timeLimit = min(perSolveSeconds,remainingGridSeconds-detail.modeling_time_s);
detail.time_limit_seconds = max(timeLimit,0);
if timeLimit<=0
    detail.status = 'grid_time_budget_during_model_build';
    detail.solver_status = 'NOT_RUN_GRID_TIME_BUDGET';
    detail.solver_info = 'Storage grid budget exhausted during YALMIP model building.';
    detail.wall_time_s = toc(attemptTimer);
    return;
end

params = struct('Threads',solverThreads,'TimeLimit',timeLimit, ...
    'TuneTimeLimit',0,'OutputFlag',double(verbose~=0));
% Basis reuse changes only the algorithmic starting point. The model and
% objective are unchanged. Gurobi crushes the start through presolve and
% dual simplex reoptimizes when epsilon/RHS/bounds change.
persistent basisSignature basisVBasis basisCBasis
signature=[T,sp.eta_ch,sp.eta_dis,sp.self_discharge_per_h, ...
    sp.soc_min,sp.soc_max,sp.duration_h];
detail.basis_reused=false;
if reuseLPBasis
    params.Method=1; params.LPWarmStart=2;
    if isequal(signature,basisSignature) && ...
            numel(basisVBasis)==numel(model.obj) && numel(basisCBasis)==numel(model.rhs)
        model.vbasis=basisVBasis; model.cbasis=basisCBasis;
        detail.basis_reused=true;
    end
end
detail.attempted = true;
solveTimer = tic;
try
    raw = gurobi(model,params);
    if reuseLPBasis && strcmpi(raw.status,'OPTIMAL') && ...
            isfield(raw,'vbasis') && isfield(raw,'cbasis')
        basisSignature=signature; basisVBasis=raw.vbasis; basisCBasis=raw.cbasis;
    end
    if strcmpi(raw.status,'OPTIMAL'); problem=0; else; problem=3; end
    res = struct('problem',problem,'info',raw.status, ...
        'solvertime',getNumericField(raw,'runtime',toc(solveTimer)), ...
        'yalmiptime',0,'solveroutput',struct('result',raw));
catch ME
    detail.status = 'solver_exception';
    detail.solver_status = 'EXCEPTION';
    detail.solver_problem = -1;
    detail.solver_info = ME.message;
    detail.solver_time_s = toc(solveTimer);
    detail.wall_time_s = toc(attemptTimer);
    return;
end
detail.solver_time_s = getNumericField(res,'solvertime',toc(solveTimer));
detail.yalmip_time_s = getNumericField(res,'yalmiptime',NaN);
detail.solver_problem = res.problem;
detail.solver_info = res.info;
raw = extractGurobiResult(res);
detail.solver_status = getTextField(raw,'status',sprintf('YALMIP_PROBLEM_%d',res.problem));
detail.solver_runtime_s = getNumericField(raw,'runtime',detail.solver_time_s);
detail.solver_objective_CNY = getNumericField(raw,'objval',NaN);
detail.objective_bound_CNY = getNumericField(raw,'objbound',NaN);
if ~isfinite(detail.objective_bound_CNY)
    detail.objective_bound_CNY = getNumericField(raw,'objboundc',NaN);
end

% YALMIP may populate zeros when a solver returns no primal vector. Require
% an explicit primal result for non-optimal solves before accepting values.
hasSolverPrimal = isfield(raw,'x')&&~isempty(raw.x)&&all(isfinite(raw.x(:)));
canReadPrimal = (res.problem==0)||hasSolverPrimal;
excludedStatus = any(strcmpi(detail.solver_status, ...
    {'INFEASIBLE','UNBOUNDED','INF_OR_UNBD','EXCEPTION'}));
if ~canReadPrimal||excludedStatus
    detail.status = 'no_feasible_primal';
    detail.wall_time_s = toc(attemptTimer);
    return;
end

Ev = raw.x(ix.Emax);
initialEnergy = raw.x(ix.Eini);
energy = raw.x(ix.energy);
delivery = raw.x(ix.delivery);
charge = raw.x(ix.charge);
discharge = raw.x(ix.discharge);
curtail = raw.x(ix.curtail);
if any(~isfinite([Ev;initialEnergy;energy(:);delivery(:);charge(:);discharge(:);curtail(:)]))
    detail.status = 'nonfinite_primal_rejected';
    detail.wall_time_s = toc(attemptTimer);
    return;
end
detail.raw_constraint_residuals = verify_storage_lp_dispatch( ...
    Pwind,epsilon,sp,Ev,initialEnergy,energy,charge,discharge,delivery,curtail,false);
if ~detail.raw_constraint_residuals.feasible
    detail.status = 'constraint_violation_rejected';
    detail.constraint_residuals = detail.raw_constraint_residuals;
    detail.wall_time_s = toc(attemptTimer);
    return;
end
% Clipping only tiny negative flow values is followed by another complete
% physical check, so numerical cleanup cannot silently invalidate dispatch.
[charge,discharge,curtail,cleanup] = normalize_storage_lp_dispatch( ...
    charge,discharge,curtail,sp.eta_ch,sp.eta_dis);
detail.simultaneous_operation_cleanup = cleanup;
detail.constraint_residuals = verify_storage_lp_dispatch( ...
    Pwind,epsilon,sp,Ev,initialEnergy,energy,charge,discharge,delivery,curtail,true);
if ~detail.constraint_residuals.feasible
    detail.status = 'cleaned_constraint_violation_rejected';
    detail.wall_time_s = toc(attemptTimer);
    return;
end

Emax = max(Ev,0);
Pdelivery = max(delivery(:),0);
% Check again with the reported capacity, which only differs when the solver
% has returned a tiny negative capacity within its feasibility tolerance.
detail.constraint_residuals = verify_storage_lp_dispatch( ...
    Pwind,epsilon,sp,Emax,initialEnergy,energy,charge,discharge,Pdelivery,curtail,true);
if ~detail.constraint_residuals.feasible
    Emax = NaN;
    Pdelivery = nan(T,1);
    detail.status = 'reported_capacity_constraint_violation_rejected';
    detail.wall_time_s = toc(attemptTimer);
    return;
end
detail.feasible = true;
detail.accepted_feasible = true;
detail.is_optimal = res.problem==0&& ...
    (strcmpi(detail.solver_status,'OPTIMAL')||isempty(fieldnames(raw)));
if detail.is_optimal
    detail.status = 'optimal';
else
    detail.status = 'feasible_unproven';
end
B = batteryLifecycleCost(Emax,sp,life_farm,r);
detail.storage_energy_MWh = Emax;
detail.initial_energy_MWh = initialEnergy;
detail.charge_MW = charge;
detail.discharge_MW = discharge;
detail.energy_MWh = energy;
detail.curtailment_MW = curtail;
detail.farm_delivery_MW = Pdelivery;
detail.initial_capex_CNY = B.initial_capex_CNY;
detail.annualized_capex_CNY_per_year = B.annualized_capex_CNY_per_year;
detail.annual_OM_CNY_per_year = B.annual_OM_CNY_per_year;
detail.lifecycle_PV_CNY = B.lifecycle_PV_CNY;
detail.first_year_export_revenue_CNY = sum(Pdelivery.*price);
detail.lifetime_export_revenue_PV_CNY = detail.first_year_export_revenue_CNY*pvBenefitFactor;
detail.objective_CNY = B.lifecycle_PV_CNY-detail.lifetime_export_revenue_PV_CNY;
detail.numerical_cleanup_objective_change_CNY = detail.objective_CNY- ...
    (Ev*(annualizedUnitCapex+annualUnitOM)*annuityPVFactor(life_farm,r)- ...
    sum(delivery(:).*price)*pvBenefitFactor);
if detail.is_optimal
    detail.objective_bound_CNY = detail.objective_CNY;
end
if isfinite(detail.objective_bound_CNY)
    detail.absolute_optimality_gap_CNY = max(0,detail.objective_CNY-detail.objective_bound_CNY);
    detail.optimality_gap = detail.absolute_optimality_gap_CNY/max(1,abs(detail.objective_CNY));
end
detail.wall_time_s = toc(attemptTimer);
end

function detail = emptyStorageDetail(T,status)
detail = struct('feasible',false,'accepted_feasible',false,'is_optimal',false, ...
    'attempted',false,'status',status,'solver_status','NOT_RUN', ...
    'solver_problem',NaN,'solver_info','','storage_energy_MWh',NaN, ...
    'charge_MW',nan(T,1),'discharge_MW',nan(T,1),'energy_MWh',nan(T+1,1), ...
    'curtailment_MW',nan(T,1),'farm_delivery_MW',nan(T,1), ...
    'objective_CNY',NaN,'solver_objective_CNY',NaN,'objective_bound_CNY',NaN, ...
    'absolute_optimality_gap_CNY',NaN,'optimality_gap',NaN, ...
    'modeling_time_s',0,'yalmip_time_s',NaN,'solver_time_s',0, ...
    'solver_runtime_s',NaN,'wall_time_s',0,'time_limit_seconds',0, ...
    'raw_constraint_residuals',struct(),'constraint_residuals',struct(), ...
    'certified_pruned',false,'certified_upper_bound_CNY',NaN, ...
    'rounding_tolerance_CNY',NaN,'no_storage_incumbent_CNY',NaN, ...
    'simultaneous_operation_cleanup',struct());
end

function [C,mask] = validateCertifiedPruning(C,epsilonSet)
% Certificates are external, frozen upper bounds, not estimated candidates.
% Require exactly one finite bound per original epsilon. A 1e-9 MW numerical
% match tolerance accommodates decimal CSV serialization and is dominated
% by the deliberately conservative monetary rounding allowance.
mask=false(size(epsilonSet));
if isempty(C); C=[]; return; end
required={'epsilon_MW','upper_bound_CNY','rounding_tolerance_CNY','no_storage_incumbent_CNY'};
for k=1:numel(required)
    assert(isfield(C,required{k}),'CertifiedStoragePruning missing required field.');
    assert(isnumeric(C.(required{k})) && all(isfinite(C.(required{k})(:))), ...
        'CertifiedStoragePruning values must be finite numeric values.');
end
epsilon=double(C.epsilon_MW(:));
upper=double(C.upper_bound_CNY(:));
tolerance=double(C.rounding_tolerance_CNY(:));
n=numel(epsilonSet);
assert(numel(epsilon)==n && numel(unique(epsilon))==n && numel(upper)==n, ...
    'CertifiedStoragePruning must map one-to-one to the complete epsilon set.');
if isscalar(tolerance); tolerance=repmat(tolerance,n,1); end
assert(numel(tolerance)==n && all(tolerance>=0), ...
    'CertifiedStoragePruning requires nonnegative rounding tolerances for each epsilon.');
assert(isscalar(C.no_storage_incumbent_CNY),'Certified no-storage incumbent must be scalar.');
[epsilon,order]=sort(epsilon); upper=upper(order); tolerance=tolerance(order);
assert(all(abs(epsilon-epsilonSet(:))<=1e-9), ...
    'CertifiedStoragePruning epsilon values do not match the original candidate grid.');
C=struct('epsilon_MW',epsilon.','upper_bound_CNY',upper.', ...
    'rounding_tolerance_CNY',tolerance.', ...
    'no_storage_incumbent_CNY',double(C.no_storage_incumbent_CNY));
mask=(C.upper_bound_CNY+C.rounding_tolerance_CNY<C.no_storage_incumbent_CNY);
end

function [model,ix] = cachedDispatchLP(T,sp)
% Columns are E_ini, E_max, charge(T), discharge(T), energy(T+1),
% delivery(T), curtailment(T). The additional nonnegative lower bounds on
% E_ini and energy are implied by production SOC bounds and E_max>=0.
persistent cacheSignature cacheModel cacheIndex
signature = [T,sp.eta_ch,sp.eta_dis,sp.self_discharge_per_h, ...
    sp.soc_min,sp.soc_max,sp.duration_h];
if isequal(signature,cacheSignature)
    model=cacheModel; ix=cacheIndex; return;
end
ix.Eini=1; ix.Emax=2;
ix.charge=(3:T+2).';
ix.discharge=(T+3:2*T+2).';
ix.energy=(2*T+3:3*T+3).';
ix.delivery=(3*T+4:4*T+3).';
ix.curtail=(4*T+4:5*T+3).';
nv=5*T+3; nr=6*T+6;
rr=[]; cc=[]; vv=[];
    function add(entriesR,entriesC,entriesV)
        rr=[rr;entriesR(:)]; cc=[cc;entriesC(:)]; vv=[vv;entriesV(:)];
    end
q=(1:T).'; u=ones(T,1);
% Charge and discharge <= E_max / duration_h.
add(q,ix.charge,u); add(q,repmat(ix.Emax,T,1),-u/sp.duration_h);
add(T+q,ix.discharge,u); add(T+q,repmat(ix.Emax,T,1),-u/sp.duration_h);
% E_bat(1)==E_ini; lower/upper SOC bounds on E_ini explicitly retained.
add(2*T+1,ix.energy(1),1); add(2*T+1,ix.Eini,-1);
add(2*T+2,ix.Emax,sp.soc_min); add(2*T+2,ix.Eini,-1);
add(2*T+3,ix.Eini,1); add(2*T+3,ix.Emax,-sp.soc_max);
% Pwind + discharge - charge == delivery + curtailment.
ix.balanceRows=2*T+3+q;
add(ix.balanceRows,ix.discharge,u); add(ix.balanceRows,ix.charge,-u);
add(ix.balanceRows,ix.delivery,-u); add(ix.balanceRows,ix.curtail,-u);
% Every hourly SOC recurrence.
recur=3*T+3+q;
add(recur,ix.energy(2:end),u);
add(recur,ix.energy(1:end-1),-(1-sp.self_discharge_per_h)*u);
add(recur,ix.charge,-sp.eta_ch*u);
add(recur,ix.discharge,u/sp.eta_dis);
% SOC lower and upper bounds for every state, including the terminal state.
qstate=(1:T+1).'; us=ones(T+1,1);
lower=4*T+3+qstate;
upper=5*T+4+qstate;
add(lower,repmat(ix.Emax,T+1,1),sp.soc_min*us);
add(lower,ix.energy,-us);
add(upper,ix.energy,us);
add(upper,repmat(ix.Emax,T+1,1),-sp.soc_max*us);
% E_bat(end)>=E_ini. No fixed half-full initial state is imposed.
add(nr,ix.Eini,1); add(nr,ix.energy(end),-1);
model=struct();
model.A=sparse(rr,cc,vv,nr,nv);
model.rhs=zeros(nr,1);
model.sense=repmat('<',nr,1);
model.sense(2*T+1)='=';
model.sense(ix.balanceRows)='=';
model.sense(recur)='=';
model.obj=zeros(nv,1);
model.lb=zeros(nv,1);
model.ub=inf(nv,1);
model.vtype=repmat('C',nv,1);
model.modelsense='min';
model.modelname='Production_storage_LP_direct_equivalent';
cacheSignature=signature; cacheModel=model; cacheIndex=ix;
end

function s = storageStatusTemplate()
s = struct('epsilon_MW',NaN,'status','','solver_status','','solver_problem',NaN, ...
    'attempted',false,'accepted_feasible',false,'is_optimal',false, ...
    'modeling_time_s',0,'yalmip_time_s',NaN,'solver_time_s',0,'solver_runtime_s',NaN, ...
    'wall_time_s',0,'time_limit_seconds',0,'objective_CNY',NaN, ...
    'objective_bound_CNY',NaN,'absolute_optimality_gap_CNY',NaN,'optimality_gap',NaN, ...
    'max_constraint_violation',NaN,'feasibility_tolerance',NaN, ...
    'certified_pruned',false,'certified_upper_bound_CNY',NaN, ...
    'rounding_tolerance_CNY',NaN,'no_storage_incumbent_CNY',NaN);
end

function s = summarizeStorageStatus(epsilon,detail)
s = storageStatusTemplate();
s.epsilon_MW = epsilon;
names = fieldnames(s);
for k = 1:numel(names)
    if isfield(detail,names{k}); s.(names{k}) = detail.(names{k}); end
end
if isfield(detail.constraint_residuals,'max_violation')
    s.max_constraint_violation = detail.constraint_residuals.max_violation;
    s.feasibility_tolerance = detail.constraint_residuals.tolerance;
end
end

function raw = extractGurobiResult(res)
raw = struct();
if isfield(res,'solveroutput')&&isstruct(res.solveroutput)
    if isfield(res.solveroutput,'result')&&isstruct(res.solveroutput.result)
        raw = res.solveroutput.result;
    elseif isfield(res.solveroutput,'status')||isfield(res.solveroutput,'x')
        raw = res.solveroutput;
    end
end
end

function v = getNumericField(s,name,fallback)
v = fallback;
if isfield(s,name)&&isnumeric(s.(name))&&isscalar(s.(name))&&isfinite(s.(name))
    v = double(s.(name));
end
end

function v = getTextField(s,name,fallback)
v = fallback;
if isfield(s,name)&&(ischar(s.(name))||(isstring(s.(name))&&isscalar(s.(name))))
    v = char(s.(name));
end
end

function [total_cost_ac,sumpower_final_ac,power_final_ac,detail] = calculateACTransmissionCost(P_max,P_actual,distance,life_farm,r,power_price,degradation,EUR2CNY)
% HVAC: final Appendix-D cable costs, annual O&M, reactive compensation, losses.
f = 50;
% columns: V, area, R(ohm/km), C(F/km), I(A), cable cost(kEUR/km)
hvac = [ ...
    220e3, 500,  48.9e-3, 136e-9, 732, 815; ...
    220e3, 630,  39.1e-3, 151e-9, 808, 850; ...
    220e3, 800,  31.9e-3, 163e-9, 879, 975; ...
    220e3,1000,  27.0e-3, 177e-9, 942,1000];
num = size(hvac,1); T = numel(P_actual);
Ccable = Inf(num,1); Creact = Inf(num,1); nc = nan(num,1); Seff = nan(num,1);
capex = Inf(num,1); annualOM = Inf(num,1); pvOM = Inf(num,1); lifecycle = Inf(num,1);
energy = zeros(num,1); Pfinal = zeros(num,T); npv = -Inf(num,1); revenue = zeros(num,1); pvRevenue = zeros(num,1);

% Eq. (46), corrected linear coefficient 40 (10^4 CNY/MW).
Csea = 4390 + 40*P_max + 23*(P_max^0.7513);
Cland = 18000;
pvifa = annuityPVFactor(life_farm,r);
revFactor = degradedPVFactor(life_farm,r,degradation);

for i = 1:num
    V = hvac(i,1); rr = hvac(i,3); cap = hvac(i,4); I = hvac(i,5);
    Ckm_10kCNY = hvac(i,6)*EUR2CNY/10; % kEUR/km -> 10^4 CNY/km
    Q = 2*pi*f*cap*V^2*distance;
    S = sqrt(3)*V*I;
    if S<=Q/2; continue; end
    Se = sqrt(S^2-(Q/2)^2)/1e6;
    if ~isfinite(Se)||Se<=0; continue; end
    n = max(1,ceil(max(P_max,0)/Se));
    Seff(i) = Se; nc(i) = n;
    Ccable(i) = n*Ckm_10kCNY*distance;
    Qtr = 0.01*P_max;
    Qtotal = n*Q/1e6;
    % 0.015 / 0.025 MGBP/MVAr converted with GBP2CNY=8.41 -> 12.615/21.025 (10^4 CNY/MVAr).
    Creact(i) = 12.615*(Qtotal/2) + 21.025*(Qtotal/2+Qtr);

    P1 = P_actual*(1-0.0004);
    Ploss = 3*(P1/(n*sqrt(3)*V)).^2*rr*distance*n;
    Pf = max(P1-Ploss,0);
    Pfinal(i,:) = Pf(:)'/1e3; % kW
    energy(i) = sum(Pf)/1e3;  % kWh in first year

    capex(i) = Ccable(i)+Csea+Cland+Creact(i);
    % Final Table 4 annual HVAC maintenance cost (10^4 CNY/yr).
    annualOM(i) = 0.005*Csea + 0.002*Cland + 0.0003*Ccable(i);
    pvOM(i) = annualOM(i)*pvifa;
    lifecycle(i) = capex(i)+pvOM(i);

    revenue(i) = sum(Pfinal(i,:)'/1000.*power_price); % CNY first year
    pvRevenue(i) = revenue(i)*revFactor;
    npv(i) = pvRevenue(i)-lifecycle(i)*1e4;
end

if all(~isfinite(npv))
    total_cost_ac = Inf; sumpower_final_ac = 0; power_final_ac = zeros(T,1); idx = NaN;
else
    [~,idx] = max(npv);
    total_cost_ac = lifecycle(idx);
    sumpower_final_ac = energy(idx);
    power_final_ac = Pfinal(idx,:)';
end

detail = struct();
detail.selected_cable_index = idx;
detail.selected_cable_label = selectedLabel('AC',idx);
detail.n_circuits = getVal(nc,idx);
detail.effective_capacity_per_circuit_MW = getVal(Seff,idx);
detail.capex_10kCNY = getVal(capex,idx);
detail.annual_maintenance_10kCNY = getVal(annualOM,idx);
detail.pv_maintenance_10kCNY = getVal(pvOM,idx);
detail.lifecycle_cost_10kCNY = total_cost_ac;
detail.total_cost_10kCNY = total_cost_ac; % backward-compatible alias
detail.annual_delivered_energy_kWh = sumpower_final_ac;
detail.annual_revenue_CNY = getVal(revenue,idx);
detail.pv_revenue_CNY = getVal(pvRevenue,idx);
detail.npv_CNY = getVal(npv,idx);
detail.candidate_npv_CNY = npv;
detail.candidate_capex_10kCNY = capex;
detail.candidate_annual_maintenance_10kCNY = annualOM;
detail.candidate_pv_maintenance_10kCNY = pvOM;
detail.candidate_lifecycle_cost_10kCNY = lifecycle;
detail.candidate_n_circuits = nc;
detail.candidate_effective_capacity_MW = Seff;
end

function [total_cost_dc,sumpower_final_dc,power_final_dc,detail] = calculateDCTransmissionCost(P_max,P_actual,distance,life_farm,r,power_price,degradation,EUR2CNY)
% VSC-HVDC: final Appendix-D cable costs, annual O&M, converter/line losses.
% columns: V, area, R, I, rated MW, cable cost(kEUR/km)
hvdc = [ ...
    300e3,1000,22.4e-3,1644,493.2,670; ...
    300e3,1200,19.2e-3,1791,537.3,730; ...
    300e3,1400,16.5e-3,1962,588.6,785; ...
    300e3,1600,14.4e-3,2123,636.9,840; ...
    300e3,2000,11.5e-3,2407,722.1,900];
num = size(hvdc,1); T = numel(P_actual);
Ccable = Inf(num,1); nc = nan(num,1); energy = zeros(num,1); Pfinal = zeros(num,T);
capex = Inf(num,1); annualOM = Inf(num,1); pvOM = Inf(num,1); lifecycle = Inf(num,1);
npv = -Inf(num,1); revenue = zeros(num,1); pvRevenue = zeros(num,1);

Csea = 21950 + 97*P_max;
Cland = 72*P_max;
pvifa = annuityPVFactor(life_farm,r);
revFactor = degradedPVFactor(life_farm,r,degradation);

for i = 1:num
    Sc = hvdc(i,5);
    n = max(1,ceil(max(P_max,0)/Sc)); nc(i) = n;
    Ckm_10kCNY = hvdc(i,6)*EUR2CNY/10; % kEUR/km -> 10^4 CNY/km
    Ccable(i) = n*Ckm_10kCNY*distance;

    loss = 0.0143;
    P1 = P_actual*(1-loss);
    V = hvdc(i,1); rr = hvdc(i,3);
    Ploss = 2*(P1/(n*V)).^2*rr*distance*n;
    Pf = max(P1-Ploss,0)*(1-loss);
    Pfinal(i,:) = Pf(:)'/1e3; % kW
    energy(i) = sum(Pf)/1e3;  % kWh first year

    capex(i) = Ccable(i)+Csea+Cland;
    % Final Table 4 annual VSC-HVDC maintenance cost (10^4 CNY/yr).
    annualOM(i) = 0.015*Csea + 0.015*Cland + 0.0003*Ccable(i);
    pvOM(i) = annualOM(i)*pvifa;
    lifecycle(i) = capex(i)+pvOM(i);

    revenue(i) = sum(Pfinal(i,:)'/1000.*power_price);
    pvRevenue(i) = revenue(i)*revFactor;
    npv(i) = pvRevenue(i)-lifecycle(i)*1e4;
end

if all(~isfinite(npv))
    total_cost_dc = Inf; sumpower_final_dc = 0; power_final_dc = zeros(T,1); idx = NaN;
else
    [~,idx] = max(npv);
    total_cost_dc = lifecycle(idx);
    sumpower_final_dc = energy(idx);
    power_final_dc = Pfinal(idx,:)';
end

detail = struct();
detail.selected_cable_index = idx;
detail.selected_cable_label = selectedLabel('DC',idx);
detail.n_circuits = getVal(nc,idx);
detail.rated_capacity_per_circuit_MW = getVal(hvdc(:,5),idx);
detail.capex_10kCNY = getVal(capex,idx);
detail.annual_maintenance_10kCNY = getVal(annualOM,idx);
detail.pv_maintenance_10kCNY = getVal(pvOM,idx);
detail.lifecycle_cost_10kCNY = total_cost_dc;
detail.total_cost_10kCNY = total_cost_dc; % backward-compatible alias
detail.annual_delivered_energy_kWh = sumpower_final_dc;
detail.annual_revenue_CNY = getVal(revenue,idx);
detail.pv_revenue_CNY = getVal(pvRevenue,idx);
detail.npv_CNY = getVal(npv,idx);
detail.candidate_npv_CNY = npv;
detail.candidate_capex_10kCNY = capex;
detail.candidate_annual_maintenance_10kCNY = annualOM;
detail.candidate_pv_maintenance_10kCNY = pvOM;
detail.candidate_lifecycle_cost_10kCNY = lifecycle;
detail.candidate_n_circuits = nc;
end

function sp = validateStorageParameters(sp,EUR2CNY)
required = {'eta_ch','eta_dis','self_discharge_per_h','soc_max','soc_min','duration_h', ...
    'energy_capex_CNY_per_kWh','platform_capex_EUR_per_kWh','battery_life_years','om_fraction_per_year'};
if ~isstruct(sp)
    error('storage_par must be a struct in the final parameter-corrected code.');
end
for k = 1:numel(required)
    if ~isfield(sp,required{k})
        error('storage_par missing field: %s',required{k});
    end
end
if ~isfield(sp,'EUR2CNY')||isempty(sp.EUR2CNY)
    sp.EUR2CNY = EUR2CNY;
end
if sp.eta_ch<=0||sp.eta_ch>1||sp.eta_dis<=0||sp.eta_dis>1
    error('Storage efficiencies must be in (0,1].');
end
if sp.soc_min<0||sp.soc_max>1||sp.soc_min>=sp.soc_max
    error('Require 0 <= soc_min < soc_max <= 1.');
end
if sp.duration_h<=0||sp.battery_life_years<=0||sp.om_fraction_per_year<0
    error('Invalid storage duration/lifetime/O&M parameter.');
end
end

function unit = combinedStorageCapexPerMWh(sp)
% Energy-system CAPEX + dedicated offshore BESS platform CAPEX.
unitEnergy = sp.energy_capex_CNY_per_kWh*1000;
unitPlatform = sp.platform_capex_EUR_per_kWh*sp.EUR2CNY*1000;
unit = unitEnergy + unitPlatform;
end

function B = batteryLifecycleCost(E_MWh,sp,projectLife,r)
unitCapex = combinedStorageCapexPerMWh(sp);
initialCapex = E_MWh*unitCapex;
crf = capitalRecoveryFactor(r,sp.battery_life_years);
annualizedCapex = initialCapex*crf;
annualOM = initialCapex*sp.om_fraction_per_year;
annualTotal = annualizedCapex+annualOM;
lifecyclePV = annualTotal*annuityPVFactor(projectLife,r);
B = struct();
B.initial_capex_CNY = initialCapex;
B.annualized_capex_CNY_per_year = annualizedCapex;
B.annual_OM_CNY_per_year = annualOM;
B.annual_total_CNY_per_year = annualTotal;
B.lifecycle_PV_CNY = lifecyclePV;
B.unit_initial_capex_CNY_per_MWh = unitCapex;
B.capital_recovery_factor = crf;
end

function f = capitalRecoveryFactor(r,n)
if r==0
    f = 1/n;
else
    f = r*(1+r)^n/((1+r)^n-1);
end
end

function f = annuityPVFactor(L,r)
if r==0
    f = L;
else
    f = (1-(1+r)^(-L))/r;
end
end

function f = degradedPVFactor(L,r,d)
y = (1:L).';
f = sum((1-d).^(y-1)./(1+r).^y);
end

function pv = presentValueDegraded(annual_value,L,r,d)
pv = annual_value*degradedPVFactor(L,r,d);
end

function price = normalizePrice(power_price,T)
price = double(power_price(:));
if isscalar(price)
    price = repmat(price,T,1);
elseif numel(price)~=T
    error('power_price must be scalar or have the same length as power_pu.');
end
if any(~isfinite(price))
    error('power_price contains non-finite values.');
end
end

function v = getVal(vec,idx)
if isempty(idx)||~isfinite(idx); v=NaN; else; v=vec(idx); end
end

function s = selectedLabel(prefix,idx)
if isempty(idx)||~isfinite(idx); s='infeasible'; else; s=sprintf('%s-%d',prefix,idx); end
end
