function [selected,info] = p44_select(capacityMW,npvCNY,retainedMask,targetMW,toleranceMW,cfg)
% Maximize total project NPV in [targetMW,targetMW+toleranceMW].
% Retained projects are fixed to one, including projects whose corrected NPV
% became nonpositive. Only positive-NPV unretained projects may be added.
% Direct Gurobi MILP: binary variables and two capacity rows; no dispatch or
% layout optimization is called. selected=[] means no verified incumbent.
if nargin<6 || isempty(cfg); cfg=struct; end
cfg=defaults(cfg);
capacityMW=double(capacityMW(:)); npvCNY=double(npvCNY(:));
rawRetained=retainedMask(:);
n=numel(capacityMW);
assert(numel(npvCNY)==n && numel(rawRetained)==n,'p44_select:SizeMismatch', ...
    'Capacity, NPV and retained-mask arrays must have the same length.');
assert(all(isfinite(capacityMW) & capacityMW>0),'p44_select:InvalidCapacity', ...
    'Candidate capacities must be finite and positive.');
assert(all(isfinite(npvCNY)),'p44_select:InvalidNPV','Candidate NPVs must be finite.');
assert(all(rawRetained==0 | rawRetained==1),'p44_select:InvalidRetainedMask', ...
    'Retained mask must contain only zero or one.');
assert(isnumeric(targetMW) && isscalar(targetMW) && isfinite(targetMW) && targetMW>=0, ...
    'p44_select:InvalidTarget','Capacity target must be a finite nonnegative scalar.');
assert(isnumeric(toleranceMW) && isscalar(toleranceMW) && isfinite(toleranceMW) && toleranceMW>=0, ...
    'p44_select:InvalidTolerance','Capacity-band tolerance must be a finite nonnegative scalar.');
retainedMask=logical(rawRetained);
capacityUnitMW=12; moneyScaleCNY=1e6;
units=round(capacityMW/capacityUnitMW);
assert(all(abs(capacityMW/capacityUnitMW-units)<=1e-8), ...
    'p44_select:NonDiscreteCapacity','Every project capacity must be an integer multiple of 12 MW.');
lowerUnits=ceil(targetMW/capacityUnitMW-1e-10);
upperUnits=floor((targetMW+toleranceMW)/capacityUnitMW+1e-10);
eligible=retainedMask | npvCNY>0;
retainedUnits=sum(units(retainedMask));
selected=[];
info=struct('status','NOT_SOLVED','solver_status','NOT_SOLVED', ...
    'stop_reason','','message','','has_feasible_solution',false, ...
    'proven_infeasible',false,'gap_target_met',false, ...
    'exact_optimum_certified',false,'mip_gap',NaN, ...
    'requested_mip_gap',cfg.SelectionMIPGap, ...
    'time_limit_seconds',cfg.SelectionTimeLimitSeconds, ...
    'runtime_seconds',0,'wall_time_seconds',0, ...
    'objective_CNY',NaN,'bound_CNY',NaN,'capacity_MW',NaN, ...
    'target_MW',targetMW,'tolerance_MW',toleranceMW, ...
    'capacity_unit_MW',capacityUnitMW,'lower_capacity_units',lowerUnits, ...
    'upper_capacity_units',upperUnits,'retained_capacity_MW',sum(capacityMW(retainedMask)), ...
    'candidate_count',n,'eligible_candidate_count',sum(eligible), ...
    'fixed_retained_count',sum(retainedMask),'selected_count',0, ...
    'new_selected_count',0,'warm_start_used',false,'node_count',NaN, ...
    'capacity_band_satisfied',false,'retained_constraints_satisfied',false);
if lowerUnits>upperUnits
    info=precheck_failure(info,'capacity_band_contains_no_12_MW_multiple'); return
elseif retainedUnits>upperUnits
    info=precheck_failure(info,'retained_capacity_exceeds_upper_bound'); return
elseif sum(units(eligible))<lowerUnits
    info=precheck_failure(info,'insufficient_positive_NPV_candidate_capacity'); return
end
if ~any(eligible & ~retainedMask)
    selected=retainedMask;
    info.status='OPTIMAL'; info.solver_status='NOT_CALLED';
    info.stop_reason='all_eligible_projects_are_fixed';
    info.has_feasible_solution=true; info.gap_target_met=true;
    info.exact_optimum_certified=true; info.mip_gap=0;
    info.objective_CNY=sum(npvCNY(selected)); info.bound_CNY=info.objective_CNY;
    info=incumbent_info(info,selected,capacityMW,retainedMask); return
end
model=struct;
model.A=sparse([units.';units.']);
model.obj=npvCNY/moneyScaleCNY;
model.rhs=[upperUnits;lowerUnits]; model.sense='<>';
model.modelsense='max'; model.vtype=repmat('B',1,n);
model.lb=double(retainedMask); model.ub=double(eligible);
if cfg.SelectionUseWarmStart
    % Ranking supplies a feasible starting point only. It never restricts the
    % exact search or changes the total-NPV objective.
    startMask=retainedMask; total=retainedUnits;
    ids=find(eligible & ~retainedMask);
    density=npvCNY(ids)./capacityMW(ids);
    [~,order]=sortrows([-density,double(ids)],[1 2]);
    for k=reshape(ids(order),1,[])
        if total>=lowerUnits; break; end
        if total+units(k)<=upperUnits
            startMask(k)=true; total=total+units(k);
        end
    end
    if total>=lowerUnits && total<=upperUnits
        model.start=double(startMask); info.warm_start_used=true;
    end
end
params=struct('OutputFlag',0,'Threads',cfg.SolverThreads, ...
    'TimeLimit',cfg.SelectionTimeLimitSeconds,'MIPGap',cfg.SelectionMIPGap, ...
    'Seed',0);
timer=tic;
try
    result=gurobi(model,params);
catch err
    info.status='SOLVER_ERROR'; info.solver_status='SOLVER_ERROR';
    info.stop_reason='solver_error'; info.message=err.message;
    info.wall_time_seconds=toc(timer); return
end
info.wall_time_seconds=toc(timer);
info.status=char(result.status); info.solver_status=info.status;
if isfield(result,'runtime'); info.runtime_seconds=result.runtime; end
if isfield(result,'nodecount'); info.node_count=result.nodecount; end
if isfield(result,'objbound'); info.bound_CNY=result.objbound*moneyScaleCNY; end
if isfield(result,'mipgap'); info.mip_gap=result.mipgap; end
info.proven_infeasible=strcmpi(info.solver_status,'INFEASIBLE');
if isfield(result,'x') && numel(result.x)==n && all(isfinite(result.x))
    trial=round(result.x(:));
    integral=all(abs(result.x(:)-trial)<=1e-5) && all(trial==0 | trial==1);
    capUnits=sum(units(logical(trial)));
    feasible=integral && capUnits>=lowerUnits && capUnits<=upperUnits && ...
        all(trial(retainedMask)==1) && all(trial(~eligible)==0);
    if feasible
        selected=logical(trial); info.has_feasible_solution=true;
        info.objective_CNY=sum(npvCNY(selected));
        info=incumbent_info(info,selected,capacityMW,retainedMask);
        if ~isfinite(info.mip_gap) && isfinite(info.bound_CNY)
            if info.objective_CNY==0 && info.bound_CNY==0
                info.mip_gap=0;
            elseif info.objective_CNY~=0
                info.mip_gap=abs(info.bound_CNY-info.objective_CNY)/abs(info.objective_CNY);
            end
        end
        info.gap_target_met=isfinite(info.mip_gap) && ...
            info.mip_gap<=cfg.SelectionMIPGap+1e-10;
        moneyTol=max(1e-3,1e-12*abs(info.objective_CNY));
        info.exact_optimum_certified=isfinite(info.bound_CNY) && ...
            abs(info.bound_CNY-info.objective_CNY)<=moneyTol;
    else
        info.status='INVALID_INCUMBENT'; info.stop_reason='incumbent_failed_exact_feasibility_check';
        info.message='Returned binary solution failed capacity, retained, or positive-NPV eligibility checks.';
        return
    end
end
if info.has_feasible_solution
    if strcmpi(info.solver_status,'OPTIMAL')
        info.stop_reason='solver_optimal_within_requested_gap';
    else
        info.stop_reason=[lower(info.solver_status) '_with_feasible_solution'];
    end
elseif info.proven_infeasible
    info.stop_reason='solver_proven_infeasible';
else
    % TIME_LIMIT without x is an unresolved search, never an infeasibility proof.
    info.stop_reason=[lower(info.solver_status) '_without_feasible_solution'];
end
end

function cfg=defaults(cfg)
assert(isstruct(cfg) && isscalar(cfg),'p44_select:InvalidConfig','cfg must be a scalar struct.');
names={'SolverThreads','SelectionTimeLimitSeconds','SelectionMIPGap','SelectionUseWarmStart'};
values={2,60,.001,true};
for k=1:numel(names)
    if ~isfield(cfg,names{k}) || isempty(cfg.(names{k})); cfg.(names{k})=values{k}; end
end
assert(isnumeric(cfg.SolverThreads) && isscalar(cfg.SolverThreads) && ...
    isfinite(cfg.SolverThreads) && cfg.SolverThreads>=1 && cfg.SolverThreads==round(cfg.SolverThreads), ...
    'p44_select:InvalidThreads','SolverThreads must be a positive integer.');
assert(isnumeric(cfg.SelectionTimeLimitSeconds) && isscalar(cfg.SelectionTimeLimitSeconds) && ...
    isfinite(cfg.SelectionTimeLimitSeconds) && cfg.SelectionTimeLimitSeconds>0, ...
    'p44_select:InvalidTimeLimit','SelectionTimeLimitSeconds must be finite and positive.');
assert(isnumeric(cfg.SelectionMIPGap) && isscalar(cfg.SelectionMIPGap) && ...
    isfinite(cfg.SelectionMIPGap) && cfg.SelectionMIPGap>=0, ...
    'p44_select:InvalidGap','SelectionMIPGap must be finite and nonnegative.');
assert(isscalar(cfg.SelectionUseWarmStart) && ...
    (cfg.SelectionUseWarmStart==0 || cfg.SelectionUseWarmStart==1), ...
    'p44_select:InvalidWarmStart','SelectionUseWarmStart must be logical zero or one.');
end

function info=precheck_failure(info,reason)
info.status='INFEASIBLE'; info.solver_status='NOT_CALLED';
info.stop_reason=reason; info.proven_infeasible=true;
end

function info=incumbent_info(info,selected,capacityMW,retainedMask)
info.capacity_MW=sum(capacityMW(selected));
info.selected_count=sum(selected); info.new_selected_count=sum(selected & ~retainedMask);
info.capacity_band_satisfied=info.capacity_MW>=info.target_MW-1e-6 && ...
    info.capacity_MW<=info.target_MW+info.tolerance_MW+1e-6;
info.retained_constraints_satisfied=all(selected(retainedMask));
end
