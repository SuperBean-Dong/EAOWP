function report = release_stage34_run(cfg)
%RELEASE_STAGE34_RUN Portable orchestration around unchanged scientific cores.
% Core algorithms are copied verbatim. This file replaces the former local
% five-machine assignments, fixed rerun lists, and archived-result imports.
required = {'rootDir','gisFile','stage12File','priceFile','baselineFile', ...
    'startIndex','endIndex','gridIDs','resume','stopOnError', ...
    'collectionSolver','ossBlockTimeLimitSeconds','ossGridTimeLimitSeconds', ...
    'ossMIPGap','ossEnablePruning','solverThreads','storageTimeLimitSeconds', ...
    'storageGridTimeLimitSeconds','storageMaxPasses','saveStage34Every'};
assert(isstruct(cfg)&&isscalar(cfg),'stage34:InvalidConfig','Use release_config(overrides).');
missing = required(~isfield(cfg,required));
assert(isempty(missing),'stage34:MissingConfig','Missing configuration: %s.',strjoin(missing,', '));
validateConfig(cfg);
assert(isfile(cfg.stage12File),'stage34:MissingStage12','Missing Stage-I file: %s',cfg.stage12File);
assert(isfile(cfg.gisFile),'stage34:MissingGIS','Missing GIS input: %s',cfg.gisFile);
assert(isfile(cfg.priceFile),'stage34:MissingPrice','Missing hourly price file: %s',cfg.priceFile);
assert(exist('yalmip','file')~=0&&exist('sdpvar','file')~=0, ...
    'stage34:MissingYALMIP','Install YALMIP separately and add it to the MATLAB path.');
assert(exist('gurobi','file')~=0,'stage34:MissingGurobi', ...
    'Install and license Gurobi separately, including its MATLAB interface.');

S = shaperead(cfg.gisFile);
assert(~isempty(S)&&all(isfield(S,{'X','Y','X_UTM','Y_UTM'})), ...
    'stage34:GISSchema','GIS needs X/Y boundary rings and X_UTM/Y_UTM landing-point coordinates.');
A = load(cfg.stage12File,'stage12_results','stage12_grid_indices','stage12_meta');
assert(isfield(A,'stage12_results')&&isstruct(A.stage12_results)&&~isempty(A.stage12_results), ...
    'stage34:Stage12Schema','stage12_results must be a nonempty structure array.');
stage12 = A.stage12_results(:);
assert(isfield(stage12,'grid_index')&&isfield(stage12,'status'), ...
    'stage34:Stage12Schema','Stage-I records need grid_index and status.');
recordIDs = double([stage12.grid_index]);
validPositions = find(isfinite(recordIDs));
validIDs = recordIDs(validPositions);
assert(all(validIDs==floor(validIDs))&&all(validIDs>=1&validIDs<=numel(S)) ...
    &&numel(unique(validIDs))==numel(validIDs),'stage34:GridIdentity', ...
    'Stage-I grid_index values must be unique, one-based original GIS indices.');
positionByID = zeros(numel(S),1); positionByID(validIDs)=validPositions;
if isfield(A,'stage12_grid_indices')
    availableIDs = double(A.stage12_grid_indices(:));
else
    availableIDs = validIDs(:);
end
assert(all(isfinite(availableIDs))&&all(availableIDs==floor(availableIDs)) ...
    &&all(availableIDs>=1&availableIDs<=numel(S)) ...
    &&numel(unique(availableIDs))==numel(availableIDs),'stage34:Stage12Indices', ...
    'stage12_grid_indices must be unique original GIS indices.');
requested = availableIDs(availableIDs>=cfg.startIndex&availableIDs<=cfg.endIndex);
if ~isempty(cfg.gridIDs)
    assert(all(ismember(cfg.gridIDs(:),availableIDs)),'stage34:UnknownGrid', ...
        'gridIDs contains a grid absent from this Stage-I run.');
    requested = requested(ismember(requested,cfg.gridIDs(:)));
end
requested = sort(requested(:));
assert(~isempty(requested),'stage34:EmptySelection','No Stage-I records match the requested range/gridIDs.');
assert(all(positionByID(requested)>0),'stage34:MissingStage12Record','A requested Stage-I record is absent.');
for grid=requested.'
    validateStage12(stage12(positionByID(grid)),grid);
end
price_CNY_per_MWh = loadPrice(cfg.priceFile);

outDir = fileparts(cfg.baselineFile);
if isempty(outDir); outDir=pwd; end
if ~isfolder(outDir); mkdir(outDir); end
checkpointDir = fullfile(outDir,'checkpoints');
if ~isfolder(checkpointDir); mkdir(checkpointDir); end
science = scientificParameters();
[signature,identity] = runIdentity(cfg,science);
baseline_meta = struct('schema_version','PORTABLE_BASELINE_V1', ...
    'model_signature',signature,'collection_solver',char(cfg.collectionSolver), ...
    'scientific_parameters',science,'input_identity',identity, ...
    'source_stage12_file',cfg.stage12File,'stage12_source_sha256',identity.stage12_sha256, ...
    'source_price_file',cfg.priceFile,'source_price_sha256',identity.price_sha256, ...
    'data_year',2023,'expected_grid_count',numel(S), ...
    'requested_grid_count',numel(requested),'requested_grid_ids',requested, ...
    'completed_grid_count',0,'all_complete',false,'updated_utc',utcTimestamp(), ...
    'checkpoint_directory',checkpointDir, ...
    'method_note','Fresh Stage3 then Stage4; finite-budget quality reported per grid; historical merged B0 is not imported.');
results = repmat(emptyRecord(),numel(S),1);
for grid=validIDs
    results(grid)=recordFromStage12(stage12(positionByID(grid)));
end
if cfg.resume&&isfile(cfg.baselineFile)
    old=load(cfg.baselineFile,'results','baseline_meta');
    assert(isfield(old,'baseline_meta')&&strcmp(old.baseline_meta.model_signature,signature), ...
        'stage34:StaleMaster','Existing baseline inputs/strategy/budgets differ; use a new output location or resume=false.');
    assert(numel(old.results)==numel(S),'stage34:MasterGridCount','Existing baseline has a different GIS grid count.');
    results=old.results(:);
end
cfg_stage34=cfg;
report=struct('requested_grid_count',numel(requested),'computed_grid_count',0, ...
    'resumed_grid_count',0,'pending_grid_count',0,'error_grid_count',0, ...
    'baseline_file',cfg.baselineFile,'checkpoint_directory',checkpointDir, ...
    'model_signature',signature);
fprintf('Stage3-4: %d grids; collection=%s; current storage purchase cost %.2f CNY/kWh.\n', ...
    numel(requested),cfg.collectionSolver,science.storage.energy_capex_CNY_per_kWh);
for k=1:numel(requested)
    grid=requested(k); s12=stage12(positionByID(grid));
    rec=recordFromStage12(s12);
    checkpoint=[];
    gridFile=fullfile(checkpointDir,sprintf('grid_%05d.mat',grid));
    gridTimer=tic;
    try
        if cfg.resume&&isfile(gridFile)
            saved=load(gridFile,'record','checkpoint','run_signature');
            assert(isfield(saved,'run_signature')&&strcmp(saved.run_signature,signature), ...
                'stage34:StaleCheckpoint','Checkpoint grid %d has different inputs/strategy/budgets.',grid);
            assert(saved.record.grid_index==grid,'stage34:CheckpointGrid','Checkpoint grid ID differs.');
            rec=saved.record;
            checkpoint=saved.checkpoint;
            if strcmp(rec.status,'ok')&&~isempty(checkpoint)&&checkpoint.candidate_coverage_complete
                results(grid)=rec;
                report.resumed_grid_count=report.resumed_grid_count+1;
                continue;
            end
        end
        if ~isfield(rec,'module3_has_feasible_solution')||~rec.module3_has_feasible_solution
            t=tic;
            solver=str2func(['optimize_oss_collection_system_' char(cfg.collectionSolver)]);
            oss=solver(science.turbine_rated_power_MW,science.collection_capacity_MW, ...
                science.collection_cost_CNY_per_m,science.shore_proxy_CNY_per_m, ...
                S(grid).X,S(grid).Y,S(grid).X_UTM,S(grid).Y_UTM,s12.turbine_xy_m, ...
                'BlockTimeLimitSeconds',cfg.ossBlockTimeLimitSeconds, ...
                'GridTimeLimitSeconds',cfg.ossGridTimeLimitSeconds,'MIPGap',cfg.ossMIPGap, ...
                'SolverThreads',cfg.solverThreads,'EnablePruning',cfg.ossEnablePruning);
            rec.module3_wall_time_s=toc(t);
            rec.beta=oss.oss_xy_m;
            rec.internal_cable_cost=oss.internal_cable_cost_CNY;
            rec.internal_collection_objective_cost=oss.search_objective_cost_CNY;
            rec.internal_cable_gap=oss.final_gap;
            rec.internal_cable_time_s=oss.total_optimization_time_s;
            rec.collection_edge_detail=oss.edge_detail;
            rec.distance_to_shore_km=oss.oss_to_shore_distance_km;
            rec.module3_optimality_status=oss.optimality_status;
            rec.module3_diagnostics=compactDiagnostics(oss.diagnostics);
            rec.module3_has_feasible_solution=oss.has_feasible_solution;
            rec.collection_solver=char(cfg.collectionSolver);
            assert(rec.module3_has_feasible_solution,'stage34:NoCollectionIncumbent','Stage3 has no feasible incumbent.');
            clear oss
        end
        powerPU=double(s12.wake_adjusted_power_pu_2023(:));
        D=struct('dataset_id',signature,'grid_ids',grid,'capacity_MW',rec.wind_farm_capacity, ...
            'distance_to_shore_km',rec.distance_to_shore_km, ...
            'platform_turbine_cost_CNY',rec.platform_turbine_cost, ...
            'internal_cable_cost_CNY',rec.internal_cable_cost,'module3_gap',rec.internal_cable_gap);
        coreCfg=struct('ModelSignature',signature,'StorageParameters',science.storage);
        if isempty(checkpoint)
            noStorageOut=callStorage(rec,powerPU,price_CNY_per_MWh,science,cfg,[],0);
            checkpoint=b0_new_checkpoint(D,1,coreCfg,noStorageOut);
            checkpoint.search_complete=false;
            checkpoint.candidate_coverage_complete=false;
            clear noStorageOut
        else
            assert(strcmp(checkpoint.model_signature,signature)&&checkpoint.grid_index==grid, ...
                'stage34:StorageCheckpoint','Storage checkpoint identity differs.');
        end
        rec.status='pending_stage4';
        writeGrid(gridFile,rec,checkpoint,signature);
        t4=tic;
        for pass=1:cfg.storageMaxPasses
            complete=[checkpoint.status_by_epsilon.is_optimal]&[checkpoint.status_by_epsilon.accepted_feasible];
            if all(complete); break; end
            checkpoint.pass_count=checkpoint.pass_count+1;
            passTimer=tic;
            for epsilonIndex=find(~complete)
                remaining=cfg.storageGridTimeLimitSeconds-toc(passTimer);
                if remaining<=0; break; end
                epsilon=checkpoint.epsilon_set_MW(epsilonIndex);
                attempt=tic;
                out=callStorage(rec,powerPU,price_CNY_per_MWh,science,cfg,epsilon,remaining);
                checkpoint=b0_update_checkpoint(checkpoint,out,epsilonIndex);
                checkpoint.wall_time_seconds=checkpoint.wall_time_seconds+toc(attempt);
                rec=applyCandidate(rec,checkpoint,science,cfg);
                writeGrid(gridFile,rec,checkpoint,signature);
                clear out
            end
        end
        rec.module4_wall_time_s=finiteOrZero(rec.module4_wall_time_s)+toc(t4);
        rec=applyCandidate(rec,checkpoint,science,cfg);
        rec.error_message='';
        rec.grid_total_time_s=finiteOrZero(rec.grid_total_time_s)+toc(gridTimer);
        writeGrid(gridFile,rec,checkpoint,signature);
        report.computed_grid_count=report.computed_grid_count+1;
        if ~strcmp(rec.status,'ok'); report.pending_grid_count=report.pending_grid_count+1; end
        fprintf('Grid %d: %s; Stage3 gap %.4g; Stage4 complete %d.\n', ...
            grid,rec.status,rec.internal_cable_gap,rec.module4_search_complete);
    catch ME
        % Input/strategy mismatches must not overwrite a different run's data.
        if startsWith(ME.identifier,'stage34:Stale')||strcmp(ME.identifier,'stage34:CheckpointGrid') ...
                ||strcmp(ME.identifier,'stage34:StorageCheckpoint')
            rethrow(ME);
        end
        rec.status='error'; rec.error_message=getReport(ME,'extended','hyperlinks','off');
        rec.grid_total_time_s=finiteOrZero(rec.grid_total_time_s)+toc(gridTimer);
        writeGrid(gridFile,rec,checkpoint,signature);
        report.error_grid_count=report.error_grid_count+1;
        fprintf(2,'Grid %d failed: %s\n',grid,ME.message);
        if cfg.stopOnError
            results(grid)=rec;
            baseline_meta=updateMetadata(baseline_meta,results);
            writeMaster(cfg.baselineFile,results,baseline_meta,price_CNY_per_MWh,cfg_stage34);
            rethrow(ME);
        end
    end
    results(grid)=rec;
    if mod(k,cfg.saveStage34Every)==0
        writeSummary(fullfile(outDir,'baseline_results_summary.csv'),results);
    end
end
baseline_meta=updateMetadata(baseline_meta,results);
writeMaster(cfg.baselineFile,results,baseline_meta,price_CNY_per_MWh,cfg_stage34);
writeSummary(fullfile(outDir,'baseline_results_summary.csv'),results);
report.completed_grid_count=baseline_meta.completed_grid_count;
report.expected_grid_count=baseline_meta.expected_grid_count;
report.all_complete=baseline_meta.all_complete;
fprintf('Baseline saved: %s (%d/%d grids complete).\n',cfg.baselineFile, ...
    report.completed_grid_count,report.expected_grid_count);
end

function validateConfig(cfg)
assert(ismember(char(cfg.collectionSolver),{'fast','solver_init'}), ...
    'stage34:CollectionSolver','collectionSolver must be fast or solver_init.');
validateattributes(cfg.startIndex,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(cfg.endIndex,{'numeric'},{'scalar','positive'});
if isfinite(cfg.endIndex); validateattributes(cfg.endIndex,{'numeric'},{'integer'}); end
assert(cfg.endIndex>=cfg.startIndex,'stage34:Range','endIndex precedes startIndex.');
if ~isempty(cfg.gridIDs)
    validateattributes(cfg.gridIDs,{'numeric'},{'vector','integer','positive','finite'});
end
for name={'ossBlockTimeLimitSeconds','ossGridTimeLimitSeconds','storageTimeLimitSeconds','storageGridTimeLimitSeconds'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','positive'});
end
for name={'solverThreads','storageMaxPasses','saveStage34Every'}
    validateattributes(cfg.(name{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
validateattributes(cfg.ossMIPGap,{'numeric'},{'scalar','finite','>=',0,'<',1});
for name={'resume','stopOnError','ossEnablePruning'}
    validateattributes(cfg.(name{1}),{'logical','numeric'},{'scalar','binary'});
end
assert(~strcmpi(char(cfg.baselineFile),char(cfg.stage12File))&& ...
    ~strcmpi(char(cfg.baselineFile),char(cfg.priceFile)), ...
    'stage34:OutputCollision','Baseline output must differ from source inputs.');
end

function s=scientificParameters()
s=struct('turbine_rated_power_MW',12,'availability_factor',.95, ...
    'life_years',25,'discount_rate',.03,'annual_degradation_rate',.0072, ...
    'EUR2CNY',8.39,'collection_capacity_MW',[40,50,63], ...
    'collection_cost_CNY_per_m',[423,474,554]*8.39, ...
    'shore_proxy_CNY_per_m',815*8.39,'export_step_MW',10,'minimum_export_MW',50);
s.storage=struct('eta_ch',.954,'eta_dis',.954,'self_discharge_per_h',.00054, ...
    'soc_min',.1,'soc_max',1,'duration_h',2,'energy_capex_CNY_per_kWh',1043.82, ...
    'platform_capex_EUR_per_kWh',0,'battery_life_years',15,'om_fraction_per_year',.04,'EUR2CNY',8.39);
end

function [signature,identity]=runIdentity(cfg,science)
identity=struct('stage12_sha256',b0_file_sha256(cfg.stage12File), ...
    'price_sha256',b0_file_sha256(cfg.priceFile),'gis_sha256',b0_file_sha256(cfg.gisFile));
[folder,name]=fileparts(cfg.gisFile);
for ext={'.shx','.dbf','.prj'}
    p=fullfile(folder,[name ext{1}]);
    if isfile(p); identity.(['gis' strrep(ext{1},'.','_') '_sha256'])=b0_file_sha256(p); end
end
codeDir=fileparts(mfilename('fullpath'));
files=dir(fullfile(codeDir,'*.m'));
priv=dir(fullfile(codeDir,'private','*.m')); files=[files;priv];
hashes=struct();
for k=1:numel(files)
    hashes.(erase(files(k).name,'.m'))=b0_file_sha256(fullfile(files(k).folder,files(k).name));
end
runtime=struct('collection_solver',char(cfg.collectionSolver), ...
    'oss_block_seconds',cfg.ossBlockTimeLimitSeconds,'oss_grid_seconds',cfg.ossGridTimeLimitSeconds, ...
    'oss_gap',cfg.ossMIPGap,'oss_pruning',logical(cfg.ossEnablePruning), ...
    'solver_threads',cfg.solverThreads,'storage_solve_seconds',cfg.storageTimeLimitSeconds, ...
    'storage_pass_seconds',cfg.storageGridTimeLimitSeconds,'storage_max_passes',cfg.storageMaxPasses);
payload=jsonencode(struct('schema','PORTABLE_STAGE34_V1','inputs',identity, ...
    'scientific_parameters',science,'runtime',runtime,'source_code',hashes));
temporary=[tempname '.json']; fid=fopen(temporary,'w'); assert(fid>=0);
fwrite(fid,unicode2native(payload,'UTF-8'),'uint8'); fclose(fid);
guard=onCleanup(@()delete(temporary));
signature=b0_file_sha256(temporary);
end

function validateStage12(r,grid)
required={'data_year','wind_farm_capacity','turbine_xy_m','wake_adjusted_power_pu_2023', ...
    'optimal_platform','platform_turbine_cost','platform_detail','turbine_count', ...
    'wind_data_index','availability_factor','annual_wind_degradation_rate'};
assert(all(isfield(r,required))&&strcmp(r.status,'ok'),'stage34:Stage12NotReady', ...
    'Grid %d is not a complete successful Stage-I record.',grid);
assert(r.grid_index==grid&&r.data_year==2023,'stage34:Stage12Year','The baseline driver requires original 2023 grid identities.');
assert(r.wind_farm_capacity==12*size(r.turbine_xy_m,1)&&r.turbine_count==size(r.turbine_xy_m,1), ...
    'stage34:Stage12Capacity','Grid %d capacity/turbine coordinates do not agree.',grid);
u=double(r.wake_adjusted_power_pu_2023(:));
assert(numel(u)==8760&&all(isfinite(u))&&all(u>=0&u<=1), ...
    'stage34:Stage12Power','Grid %d needs 8760 finite wake/availability-adjusted p.u. values in [0,1].',grid);
assert(abs(r.availability_factor-.95)<1e-12&&isfinite(r.platform_turbine_cost), ...
    'stage34:Stage12Parameters','Stage-I availability or platform lifecycle cost is invalid.');
assert(abs(r.annual_wind_degradation_rate-.0072)<1e-12, ...
    'stage34:Stage12Parameters','Stage-I annual degradation differs from the baseline.');
end

function p=loadPrice(file)
W=whos('-file',file); names={W.name};
preferred={'price_CNY_per_MWh','Price_yuan_per_MWh','power_price','price'}; selected='';
for k=1:numel(preferred)
    hit=find(strcmpi(names,preferred{k}),1);
    if ~isempty(hit); selected=names{hit}; break; end
end
if isempty(selected)
    candidates=find(arrayfun(@(w)ismember(w.class,{'double','single'})&&prod(w.size)>=8760,W));
    assert(numel(candidates)==1,'stage34:PriceSchema', ...
        'Price MAT must contain a recognised price variable or one unambiguous numeric vector with at least 8760 values.');
    selected=W(candidates).name;
end
A=load(file,selected); p=double(A.(selected)(:));
assert(numel(p)>=8760,'stage34:PriceVector', ...
    'Provide at least 8760 first-year hourly prices in CNY/MWh.');
p=p(1:8760);
assert(all(isfinite(p)),'stage34:PriceVector', ...
    'The first 8760 hourly prices must be finite CNY/MWh values.');
end

function out=callStorage(rec,powerPU,price,science,cfg,epsilonSet,gridBudget)
out=optimize_transmission_storage_lp_b0(rec.wind_farm_capacity,rec.distance_to_shore_km, ...
    science.storage,powerPU,science.life_years,science.discount_rate,price, ...
    'AnnualWindDegradationRate',science.annual_degradation_rate,'EUR2CNY',science.EUR2CNY, ...
    'ExportStepMW',science.export_step_MW,'MinimumExportCandidateMW',science.minimum_export_MW, ...
    'StorageTimeLimitSeconds',cfg.storageTimeLimitSeconds,'StorageGridTimeLimitSeconds',gridBudget, ...
    'SolverThreads',cfg.solverThreads,'ExportCapacityCandidatesMW',epsilonSet,'ReuseLPBasis',false);
end

function rec=applyCandidate(rec,C,science,cfg)
b=C.selected_candidate;
rec.selected_candidate=b;
rec.actual_wind_peak_MW=rec.wind_farm_capacity*max(double(rec.wake_adjusted_power_pu_2023(:)));
rec.epsilon_set_MW=C.epsilon_set_MW;
rec.candidate_summary=C.candidate_summary;
rec.E_bat_best=b.storage_energy_MWh; rec.optimal_epsilon=b.epsilon_MW;
rec.storage_case=b.storage_case; rec.export_mode=b.export_mode;
rec.selected_transmission_detail=b.transmission_detail;
for f={'battery_initial_capex_CNY','battery_annualized_capex_CNY_per_year', ...
        'battery_annual_OM_CNY_per_year','battery_cost_CNY','export_system_capex_CNY', ...
        'export_annual_maintenance_CNY','export_maintenance_pv_CNY','export_system_cost_CNY'}
    rec.(f{1})=b.(f{1});
end
rec.stage2_npv_CNY=b.stage2_npv_CNY;
rec.total_cost_best=rec.platform_turbine_cost+rec.internal_cable_cost+b.battery_cost_CNY+b.export_system_cost_CNY;
rec.NPV_best=b.pv_revenue_CNY-rec.total_cost_best;
y=(1:science.life_years).';
rec.discounted_onshore_energy_kWh=b.annual_delivered_energy_kWh* ...
    sum((1-science.annual_degradation_rate).^(y-1)./(1+science.discount_rate).^y);
rec.LCOE_min=rec.total_cost_best/max(rec.discounted_onshore_energy_kWh,eps);
rec.module4_storage_status=C.status_by_epsilon;
rec.module4_search_complete=C.candidate_coverage_complete;
rec.module4_optimality_status='incomplete_storage_candidate_search';
if C.candidate_coverage_complete
    rec.status='ok'; rec.module4_optimality_status='all_storage_LPs_optimal_on_epsilon_grid';
else
    rec.status='pending_stage4';
end
rec.optimization_complete=isfinite(rec.internal_cable_gap)&& ...
    rec.internal_cable_gap<=cfg.ossMIPGap+1e-8&&C.candidate_coverage_complete;
if rec.optimization_complete
    rec.solution_quality='collection_gap_target_met_and_storage_LPs_optimal';
else
    rec.solution_quality='feasible_with_unproven_or_incomplete_search';
end
end

function r=emptyRecord()
r=struct('grid_index',NaN,'status','not_requested','stage12_status','', ...
    'data_year',2023,'availability_factor',NaN,'annual_wind_degradation_rate',.0072, ...
    'wind_data_index',[],'wind_farm_capacity',NaN,'turbine_count',NaN,'turbine_xy_m',[], ...
    'raw_power_pu_2023',[],'wake_adjusted_power_pu_2023',[], ...
    'AEP_raw_MWh',NaN,'AEP_wake_adjusted_MWh',NaN,'wake_loss_fraction',NaN, ...
    'wake_loss_model_fraction',NaN,'pareto_capacity_MW',[],'pareto_wake_loss_fraction',[], ...
    'pareto_AEP_MWh',[],'selected_pareto_index',NaN,'optimal_platform','', ...
    'platform_turbine_cost',NaN,'platform_detail',struct(), ...
    'beta',[],'internal_cable_cost',NaN,'internal_collection_objective_cost',NaN, ...
    'internal_cable_gap',NaN,'internal_cable_time_s',NaN,'collection_edge_detail',struct([]), ...
    'distance_to_shore_km',NaN,'actual_wind_peak_MW',NaN,'epsilon_set_MW',[], ...
    'E_bat_best',NaN,'optimal_epsilon',NaN,'storage_case','','export_mode','', ...
    'selected_transmission_detail',struct(),'selected_candidate',struct(),'candidate_summary',struct(), ...
    'battery_initial_capex_CNY',NaN,'battery_annualized_capex_CNY_per_year',NaN, ...
    'battery_annual_OM_CNY_per_year',NaN,'battery_cost_CNY',NaN, ...
    'export_system_capex_CNY',NaN,'export_annual_maintenance_CNY',NaN, ...
    'export_maintenance_pv_CNY',NaN,'export_system_cost_CNY',NaN,'stage2_npv_CNY',NaN, ...
    'NPV_best',NaN,'total_cost_best',NaN,'discounted_onshore_energy_kWh',NaN,'LCOE_min',NaN, ...
    'module3_wall_time_s',NaN,'module4_wall_time_s',NaN,'grid_total_time_s',NaN, ...
    'collection_solver','','module3_has_feasible_solution',false,'module3_optimality_status','', ...
    'module3_diagnostics',struct(),'module4_optimality_status','', ...
    'module4_storage_status',struct([]),'module4_search_complete',false, ...
    'optimization_complete',false,'solution_quality','','error_message','');
end

function r=recordFromStage12(s)
r=emptyRecord();
names=intersect(fieldnames(s),fieldnames(r));
for k=1:numel(names); r.(names{k})=s.(names{k}); end
r.stage12_status=s.status; r.status='not_requested'; r.error_message='';
end

function d=compactDiagnostics(d)
if isfield(d,'block_records')&&~isempty(d.block_records)
    heavy={'x','y','P_MW','t_oss_to_node_m'}; heavy=heavy(isfield(d.block_records,heavy));
    if ~isempty(heavy); d.block_records=rmfield(d.block_records,heavy); end
end
end

function writeGrid(file,record,checkpoint,run_signature)
temporary=[file '.tmp.mat'];
save(temporary,'record','checkpoint','run_signature','-v7.3');
[ok,message]=movefile(temporary,file,'f'); assert(ok,'stage34:CheckpointWrite','%s',message);
end

function m=updateMetadata(m,results)
m.completed_grid_count=nnz(strcmp({results.status},'ok'));
m.pending_grid_count=nnz(strcmp({results.status},'pending_stage4'));
m.error_grid_count=nnz(strcmp({results.status},'error'));
m.all_complete=m.completed_grid_count==m.expected_grid_count;
m.updated_utc=utcTimestamp();
end

function writeMaster(file,results,baseline_meta,price_CNY_per_MWh,cfg_stage34)
temporary=[file '.tmp.mat'];
save(temporary,'results','baseline_meta','price_CNY_per_MWh','cfg_stage34','-v7.3');
[ok,message]=movefile(temporary,file,'f'); assert(ok,'stage34:MasterWrite','%s',message);
end

function writeSummary(file,results)
T=table([results.grid_index].',string({results.status}).',[results.wind_farm_capacity].', ...
    [results.NPV_best].',[results.LCOE_min].',[results.internal_cable_gap].', ...
    [results.module4_search_complete].',[results.optimization_complete].', ...
    'VariableNames',{'grid_index','status','wind_farm_capacity_MW','NPV_best_CNY', ...
    'LCOE_selected_CNY_per_kWh','collection_joint_objective_gap', ...
    'stage4_complete','collection_gap_and_stage4_targets_met'});
writetable(T,file,'Encoding','UTF-8');
end

function x=finiteOrZero(x)
if ~isfinite(x); x=0; end
end

function value=utcTimestamp()
value=char(datetime('now','TimeZone','UTC','Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
end
