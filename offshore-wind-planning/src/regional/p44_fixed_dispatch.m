function dispatch = p44_fixed_dispatch(Pwind_MW,r,selectedCandidate,price,cfg)
% One fixed-equipment dispatch; the objective follows the original B0 LP.
% It maximizes farm-side hourly revenue. Shore losses are evaluated afterward.
Pwind_MW=double(Pwind_MW(:)); T=numel(Pwind_MW);
assert(T>0&&all(isfinite(Pwind_MW))&&all(Pwind_MW>=0),'Invalid available power.');
price=double(price(:));
if isscalar(price); price=repmat(price,T,1); end
assert(numel(price)==T&&all(isfinite(price))&&all(price>0),'Invalid hourly price.');
sp=cfg.StorageParameters;
E=double(r.E_bat_best); epsilon=double(r.optimal_epsilon);
assert(isfinite(E)&&E>=0&&isfinite(epsilon)&&epsilon>0,'Invalid fixed equipment.');
assert(abs(E-selectedCandidate.storage_energy_MWh)<=max(1e-7,1e-9*E), ...
    'Battery capacity differs from selected B0 equipment.');
assert(abs(epsilon-selectedCandidate.epsilon_MW)<=max(1e-7,1e-9*epsilon), ...
    'Export capacity differs from selected B0 equipment.');
assert(strcmp(char(r.export_mode),char(selectedCandidate.export_mode)), ...
    'Export mode differs from selected B0 equipment.');
td=selectedCandidate.transmission_detail;
if strcmp(char(r.export_mode),'HVAC')
    hardwareCapacity=td.n_circuits*td.effective_capacity_per_circuit_MW;
elseif strcmp(char(r.export_mode),'VSC-HVDC')
    hardwareCapacity=td.n_circuits*td.rated_capacity_per_circuit_MW;
else
    error('p44:UnknownExportMode','Unknown export mode: %s',char(r.export_mode));
end
assert(isfinite(hardwareCapacity)&&hardwareCapacity>0,'Invalid cable circuit capacity.');
dispatchLimit=min(epsilon,double(hardwareCapacity));
dispatch=struct('Pwind_MW',Pwind_MW,'storage_energy_MWh',E, ...
    'epsilon_MW',epsilon,'dispatch_export_limit_MW',dispatchLimit, ...
    'fixed_hardware_capacity_MW',hardwareCapacity,'reused_baseline',false, ...
    'is_optimal',false,'feasible',false,'attempted',false,'status','not_run', ...
    'solver_status','NOT_RUN','wall_time_s',0,'objective_basis', ...
    'original_B0_farm_side_hourly_revenue_LP_shore_losses_evaluated_afterward');
timer=tic;

% Preserve the saved dispatch exactly when the hourly available power is unchanged.
sd=selectedCandidate.storage_detail;
baselineWind=selectedCandidate.farm_delivery_MW(:)+sd.curtailment_MW(:)+ ...
    sd.charge_MW(:)-sd.discharge_MW(:);
assert(numel(baselineWind)==T&&all(isfinite(baselineWind)), ...
    'Original B0 hourly horizon differs from available power.');
dispatch.baseline_constraint_residuals=p44_verify_storage_dispatch( ...
    baselineWind,dispatchLimit,sp,E,sd.energy_MWh(1),sd.energy_MWh(:), ...
    sd.charge_MW(:),sd.discharge_MW(:),selectedCandidate.farm_delivery_MW(:), ...
    sd.curtailment_MW(:),true);
assert(dispatch.baseline_constraint_residuals.feasible, ...
    'Selected B0 snapshot is incompatible with the fixed equipment/physical parameters.');
tol=1e-8+1e-10*max([1;abs(Pwind_MW)]);
if numel(baselineWind)==T&&max(abs(baselineWind-Pwind_MW))<=tol
    dispatch.farm_delivery_MW=selectedCandidate.farm_delivery_MW(:);
    dispatch.charge_MW=sd.charge_MW(:); dispatch.discharge_MW=sd.discharge_MW(:);
    dispatch.energy_MWh=sd.energy_MWh(:); dispatch.curtailment_MW=sd.curtailment_MW(:);
    dispatch.initial_energy_MWh=dispatch.energy_MWh(1);
    dispatch.reused_baseline=true; dispatch.is_optimal=logical(sd.is_optimal);
    dispatch.status='reused_saved_B0_dispatch'; dispatch.solver_status='NOT_RUN_B0_REUSE';
elseif E==0
    dispatch.farm_delivery_MW=min(Pwind_MW,dispatchLimit);
    dispatch.charge_MW=zeros(T,1); dispatch.discharge_MW=zeros(T,1);
    dispatch.energy_MWh=zeros(T+1,1);
    dispatch.curtailment_MW=Pwind_MW-dispatch.farm_delivery_MW;
    dispatch.initial_energy_MWh=0; dispatch.is_optimal=true;
    dispatch.status='explicit_fixed_zero_storage'; dispatch.solver_status='NOT_RUN_NO_STORAGE';
else
    [model,ix]=p44_cached_fixed_dispatch_lp(T,sp);
    model.lb(ix.Emax)=E; model.ub(ix.Emax)=E;
    model.ub(ix.delivery)=dispatchLimit;
    M=max(max(Pwind_MW),dispatchLimit);
    model.ub(ix.charge)=M; model.ub(ix.discharge)=M;
    model.rhs(ix.balanceRows)=-Pwind_MW;
    years=(1:cfg.LifeFarmYears).';
    benefit=sum((1-cfg.AnnualWindDegradationRate).^(years-1)./ ...
        (1+cfg.AnnualInterestRate).^years);
    model.obj(ix.delivery)=-price*benefit;
    assert(exist('gurobi','file')~=0,'Gurobi MATLAB interface is required for fixed positive BESS.');
    params=struct('Threads',cfg.SolverThreads,'TimeLimit',cfg.StorageTimeLimitSeconds, ...
        'TuneTimeLimit',0,'OutputFlag',0);
    dispatch.attempted=true;
    raw=gurobi(model,params);
    dispatch.solver_status=char(raw.status);
    assert(strcmpi(raw.status,'OPTIMAL'), ...
        'p44:FixedDispatchNotOptimal','Fixed BESS dispatch did not reach OPTIMAL (%s).',char(raw.status));
    assert(isfield(raw,'x')&&numel(raw.x)==numel(model.obj)&&all(isfinite(raw.x)), ...
        'Missing finite fixed-dispatch primal solution.');
    dispatch.initial_energy_MWh=double(raw.x(ix.Eini));
    dispatch.energy_MWh=double(raw.x(ix.energy));
    delivery=double(raw.x(ix.delivery));
    charge=double(raw.x(ix.charge)); discharge=double(raw.x(ix.discharge));
    curtail=double(raw.x(ix.curtail));
    rawReport=p44_verify_storage_dispatch(Pwind_MW,dispatchLimit,sp,E, ...
        dispatch.initial_energy_MWh,dispatch.energy_MWh,charge,discharge,delivery,curtail,false);
    assert(rawReport.feasible,'Fixed-dispatch raw physical constraints failed.');
    [charge,discharge,curtail,cleanup]=p44_normalize_storage_dispatch( ...
        charge,discharge,curtail,sp.eta_ch,sp.eta_dis);
    dispatch.farm_delivery_MW=max(delivery(:),0);
    dispatch.charge_MW=charge; dispatch.discharge_MW=discharge;
    dispatch.curtailment_MW=curtail;
    dispatch.simultaneous_operation_cleanup=cleanup;
    dispatch.raw_constraint_residuals=rawReport;
    dispatch.is_optimal=true; dispatch.status='optimal_fixed_equipment_LP';
    if isfield(raw,'objval'); dispatch.solver_objective_CNY=double(raw.objval); end
    if isfield(raw,'runtime'); dispatch.solver_runtime_s=double(raw.runtime); end
end
dispatch.constraint_residuals=p44_verify_storage_dispatch(Pwind_MW,dispatchLimit,sp,E, ...
    dispatch.initial_energy_MWh,dispatch.energy_MWh,dispatch.charge_MW,dispatch.discharge_MW, ...
    dispatch.farm_delivery_MW,dispatch.curtailment_MW,true);
assert(dispatch.constraint_residuals.feasible,'Reported fixed-dispatch physical constraints failed.');
dispatch.feasible=true;
dispatch.first_year_farm_revenue_CNY=sum(dispatch.farm_delivery_MW.*price);
dispatch.wall_time_s=toc(timer);
end
