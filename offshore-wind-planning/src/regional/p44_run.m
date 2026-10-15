function reports=p44_run(overrides)
% Initial MILP -> boundary correction -> fixed-design revaluation -> supplement.
if nargin<1; overrides=struct(); end
cfg=p44_config(overrides); D=load(cfg.InputFile);
assert(strcmp(D.schema,'P44_B0_INPUT_V1'));
meta=jsondecode(fileread(fullfile(cfg.BaselineDir,'baseline_metadata.json')));
manifest=jsondecode(fileread(cfg.ManifestFile));
merge=jsondecode(fileread(fullfile(cfg.BaselineDir,'merge_status.json')));
assert(merge.num_pending==0&&merge.all_stage4_complete&&merge.all_candidate_tables_complete,'B0 merge is incomplete.');
assert(strcmp(meta.model_signature,D.model_signature),'Prepared B0 input is stale.');
assert(strcmpi(cfg.BaselineDir,D.baseline_dir),'Prepare inputs for the requested B0 directory.');
source=dir(D.baseline_file); assert(~isempty(source)&&source.bytes==D.baseline_bytes,'B0 source changed.');
assert(strcmp(p44_sha256(fullfile(cfg.BaselineDir,'baseline_results_summary.csv'),true), ...
    manifest.source_summary_sha256),'B0 source summary changed; prepare inputs again.');
price=load(cfg.PriceFile,'price_CNY_per_MWh'); cfg.PriceCNYPerMWh=double(price.price_CNY_per_MWh(:));
assert(numel(cfg.PriceCNYPerMWh)==8760&&all(isfinite(cfg.PriceCNYPerMWh))&&all(cfg.PriceCNYPerMWh>0));
frozen=load(manifest.frozen_price_source,'price_CNY_per_MWh');
assert(strcmp(p44_sha256(manifest.frozen_price_source,true),manifest.frozen_price_source_sha256));
assert(isequal(cfg.PriceCNYPerMWh,double(frozen.price_CNY_per_MWh(:))),'Price differs from the frozen B0 vector.');
cfg.BaselineModelSignature=meta.model_signature;
science=struct('schema','P44_FIXED_DESIGN_V1','baseline',meta.model_signature, ...
    'input_hash',p44_sha256(cfg.InputFile,true),'price_hash',p44_sha256(cfg.PriceFile,true), ...
    'wind_file',cfg.WindFile,'wind_stat',dir(cfg.WindFile), ...
    'life',cfg.LifeFarmYears,'discount',cfg.AnnualInterestRate,'degradation',cfg.AnnualWindDegradationRate, ...
    'spacing',cfg.MinimumSpacingM,'spacing_tolerance',cfg.SpacingNumericalToleranceM, ...
    'cost_policy','fixed full development; junction original support-facility proxy', ...
    'zero_grid_policy',cfg.ZeroGridPolicy);
code=dir(fullfile(cfg.SourceCodeDir,'*.m'));
for k=1:numel(code); science.code_hashes.(erase(code(k).name,'.m'))=p44_sha256(fullfile(code(k).folder,code(k).name),true); end
priv=dir(fullfile(cfg.SourceCodeDir,'private','*.m'));
for k=1:numel(priv); science.code_hashes.(erase(priv(k).name,'.m'))=p44_sha256(fullfile(priv(k).folder,priv(k).name),true); end
cfg.ScienceSignature=p44_sha256(jsonencode(science));
assert(~startsWith(lower(cfg.ResultDir),lower(cfg.BaselineDir))&&~startsWith(lower(cfg.CacheDir),lower(cfg.BaselineDir)), ...
    'Result/cache directories must be separate from source B0.');
if ~isfolder(cfg.ResultDir); mkdir(cfg.ResultDir); end
if ~isfolder(cfg.CacheDir); mkdir(cfg.CacheDir); end
nPoints=size(D.point_xy_m,1); pairs=double(D.conflict_pairs);
adjacency=sparse([pairs(:,1);pairs(:,2)],[pairs(:,2);pairs(:,1)],1,nPoints,nPoints);
reports=cell(numel(cfg.TargetCapacityMW),1);
for t=1:numel(cfg.TargetCapacityMW)
    target=cfg.TargetCapacityMW(t);
    folder=fullfile(cfg.ResultDir,sprintf('target_%gGW',target/1000));
    if ~isfolder(folder); mkdir(folder); end
    stateFile=fullfile(folder,'resume_state.mat');
    runSignature=p44_sha256([cfg.ScienceSignature '|' num2str(target,17) '|' num2str(cfg.CapacityToleranceMW,17)]);
    if cfg.Resume&&isfile(stateFile)
        old=load(stateFile,'state'); state=old.state;
        assert(strcmp(state.signature,runSignature),'Resume inputs/settings changed; choose a new ResultDir or Resume=false.');
        if strcmp(state.status,'selection_unresolved') || ...
                (strcmp(state.status,'iteration_limit_reached') && cfg.MaxIterations>state.iteration)
            state.status='running'; state.phase='select';
        elseif strcmp(state.status,'capacity_target_unreachable_geometry_only') && cfg.EvaluateFailedPortfolio
            state.status='running'; state.phase='evaluate';
        end
    else
        state=newState(D,target,runSignature);
    end
    fprintf('\nP44 %g GW: phase=%s, retained=%d, capacity=%.3f GW\n',target/1000,state.phase,sum(state.retained),sum(state.capacity(state.retained))/1000);
    try
        while strcmp(state.status,'running')
            if strcmp(state.phase,'select')
                if state.iteration>=cfg.MaxIterations
                    state.status='iteration_limit_reached'; break;
                end
                state.iteration=state.iteration+1;
                [chosen,selection]=p44_select(state.capacity,state.npv,state.retained,target,cfg.CapacityToleranceMW,cfg);
                state.selection_history{end+1}=selection;
                if ~selection.has_feasible_solution
                    if selection.proven_infeasible
                        state.status='no_feasible_supplement_under_retention_rule';
                    else
                        state.status='selection_unresolved';
                    end
                    break;
                end
                new=find(chosen & ~state.retained);
                if isempty(new)
                    state.phase='check'; continue;
                end
                if isempty(state.initial_selected_ids)
                    state.initial_selected_ids=D.grid_ids(chosen);
                    state.initial_selected_capacity_MW=sum(state.capacity(chosen));
                    state.initial_selection_NPV_CNY=sum(state.npv(chosen));
                end
                [state.keepGlobal,keepByGrid,repair]=p44_spatial_repair(D,adjacency,state.keepGlobal,new);
                state.repair_history{end+1}=repair;
                zero=repair.zero_survivor_grid_indices;
                state.excluded(zero)=true; state.npv(zero)=0;
                new=setdiff(new,zero,'stable');
                state.proposed=state.retained; state.proposed(new)=true;
                state.pending=new(:);
                corrected=state.capacity;
                for g=new(:).'; corrected(g)=sum(keepByGrid{g})*12; end
                remainder=~state.proposed & ~state.excluded;
                capacityUpper=sum(corrected(state.proposed))+sum(state.capacity(remainder));
                state.remaining_capacity_upper_bound_MW=capacityUpper;
                fprintf('  round %d: new=%d, deleted=%d, corrected=%.3f GW, remaining upper=%.3f GW\n', ...
                    state.iteration,numel(new),repair.deleted_total,sum(corrected(state.proposed))/1000,capacityUpper/1000);
                if capacityUpper<target-1e-6 && ~cfg.EvaluateFailedPortfolio
                    % Capacity alone proves the target unreachable for this retained
                    % layout; expensive financial revaluation cannot change it.
                    for g=new(:).'
                        state.rows{g}=geometryOnlySummary(D.records{g},keepByGrid{g});
                        state.capacity(g)=corrected(g);
                    end
                    state.retained=state.proposed;
                    state.status='capacity_target_unreachable_geometry_only';
                    state.phase='finished'; saveState(stateFile,state); break;
                end
                state.phase='evaluate'; saveState(stateFile,state);
            end
            if strcmp(state.phase,'evaluate')
                while ~isempty(state.pending)
                    count=min(numel(state.pending),max(1,cfg.RecomputeWorkers));
                    batch=state.pending(1:count); rows=cell(count,1);
                    batchRecords=D.records(batch);
                    masks=cell(count,1);
                    for b=1:count
                        g=batch(b); pts=double(D.point_start(g)):double(D.point_start(g+1))-1;
                        masks{b}=state.keepGlobal(pts);
                    end
                    if cfg.RecomputeWorkers>1
                        assert(license('test','Distrib_Computing_Toolbox'),'Parallel Computing Toolbox required; set RecomputeWorkers=1.');
                        pool=gcp('nocreate');
                        if isempty(pool); parpool('local',cfg.RecomputeWorkers); end
                        parfor b=1:count
                            rows{b}=p44_revalue_grid(batchRecords{b},masks{b},cfg);
                        end
                    else
                        rows{1}=p44_revalue_grid(batchRecords{1},masks{1},cfg);
                    end
                    for b=1:count
                        g=batch(b); state.rows{g}=rows{b};
                        state.capacity(g)=rows{b}.capacity_MW; state.npv(g)=rows{b}.NPV_CNY;
                        state.retained(g)=true;
                    end
                    state.pending(1:count)=[]; saveState(stateFile,state);
                    if mod(numel(state.pending),10)<count || isempty(state.pending)
                        fprintf('    financial grids remaining=%d; repaired cache entries reused where available\n',numel(state.pending));
                    end
                end
                state.phase='check'; saveState(stateFile,state);
            end
            if strcmp(state.phase,'check')
                violations=sum(state.keepGlobal(pairs(:,1)) & state.keepGlobal(pairs(:,2)));
                assert(violations==0,'Cross-grid spacing repair left conflicts.');
                capacity=sum(state.capacity(state.retained));
                if capacity>=target-1e-6&&capacity<=target+cfg.CapacityToleranceMW+1e-6
                    state.status='capacity_and_spacing_satisfied'; state.phase='finished';
                else
                    assert(capacity<target,'Corrected capacity unexpectedly exceeds upper band.');
                    state.phase='select';
                end
                saveState(stateFile,state);
            end
        end
    catch err
        state.last_error=struct('identifier',err.identifier,'message',err.message,'grid_ids_pending',D.grid_ids(state.pending));
        saveState(stateFile,state);
        failed=struct('status','revaluation_error_resume_available','message',err.message, ...
            'state_file',stateFile,'target_GW',target/1000);
        writeJSON(fullfile(folder,'run_error.json'),failed);
        rethrow(err);
    end
    saveState(stateFile,state);
    reports{t}=p44_export_results(D,state,cfg,folder,science);
    if isfile(fullfile(folder,'run_error.json'))
        writeJSON(fullfile(folder,'run_error.json'),struct('resolved',true,'current_status',state.status));
    end
    fprintf('  %s: %.3f GW; shortfall %.3f GW\n',reports{t}.status,reports{t}.actual_capacity_GW,reports{t}.capacity_shortfall_GW);
end
end
function state=newState(D,target,signature)
n=numel(D.records);
state=struct('schema','P44_RUN_STATE_V1','signature',signature,'target_MW',target, ...
    'status','running','phase','select','iteration',0,'retained',false(n,1), ...
    'excluded',false(n,1),'proposed',false(n,1),'pending',[], ...
    'keepGlobal',false(size(D.point_xy_m,1),1),'rows',{cell(n,1)}, ...
    'capacity',cellfun(@(r)double(r.wind_farm_capacity),D.records), ...
    'npv',cellfun(@(r)double(r.NPV_best),D.records), ...
    'selection_history',{{}},'repair_history',{{}},'initial_selected_ids',[], ...
    'initial_selected_capacity_MW',NaN,'initial_selection_NPV_CNY',NaN, ...
    'remaining_capacity_upper_bound_MW',NaN,'last_error',struct());
end
function q=geometryOnlySummary(r,keep)
q=struct('grid_id',r.grid_index,'original_capacity_MW',r.wind_farm_capacity, ...
    'capacity_MW',sum(keep)*12,'original_turbine_count',r.turbine_count, ...
    'turbine_count',sum(keep),'deleted_turbine_count',r.turbine_count-sum(keep), ...
    'kept_original_indices',find(keep),'original_NPV_CNY',r.NPV_best, ...
    'NPV_CNY',NaN,'LCOE_CNY_per_kWh',NaN,'AEP_wake_adjusted_MWh',NaN, ...
    'annual_onshore_delivery_MWh',NaN,'farm_cost_CNY',NaN,'collection_cost_CNY',NaN, ...
    'total_cost_CNY',NaN,'junction_count',NaN,'junction_proxy_cost_CNY',NaN, ...
    'baseline_reused',false,'cache_file','','cache_hit',false,'financial_evaluated',false, ...
    'recompute_seconds',0,'dispatch_status','not_evaluated_capacity_precheck_failed');
end
function saveState(file,state)
temporary=[file '.tmp.mat']; save(temporary,'state','-v7.3'); movefile(temporary,file,'f');
end
function writeJSON(file,data)
fid=fopen(file,'w','n','UTF-8'); assert(fid>0); c=onCleanup(@()fclose(fid));
fwrite(fid,unicode2native(jsonencode(data,'PrettyPrint',true),'UTF-8'),'uint8');
end
