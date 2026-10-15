function report = run_stage1(cfg)
%RUN_STAGE1 Run the paper's baseline 7D turbine and platform design modules.
%   REPORT = RUN_STAGE1() uses RELEASE_CONFIG().
%   REPORT = RUN_STAGE1(CFG) accepts a complete configuration returned by
%   RELEASE_CONFIG. Inputs are supplied by the user; no research data are
%   included in this source-only package. See docs/STAGE1_INPUTS.md.
%
% Scientific functions and numerical parameters retain the archived model.
% This entry point replaces local paths and the old machine-specific slice.
% REPORT contains completion counts and the saved master-file path.

if nargin < 1 || isempty(cfg)
    cfg = release_config();
end
assert(isstruct(cfg) && isscalar(cfg),'run_stage1:InvalidConfig', ...
    'cfg must be a scalar configuration structure from release_config.');
requiredConfig = {'rootDir','dataDir','outputDir','gisFile', ...
    'speed140File','direction140File','compactWindFile','stage12File', ...
    'startIndex','endIndex','chunkSize','resume','useParallelMicrositing', ...
    'numWorkers','stopOnError','rebuildCompactWind'};
missingConfig = requiredConfig(~isfield(cfg,requiredConfig));
assert(isempty(missingConfig),'run_stage1:MissingConfig', ...
    'Configuration fields missing: %s. Use release_config(overrides).', ...
    strjoin(missingConfig,', '));

codeDir = fileparts(mfilename('fullpath'));
previousPath = path;
pathGuard = onCleanup(@()path(previousPath)); %#ok<NASGU>
addpath(fullfile(codeDir,'src','stage1','functions'),'-begin');

% These aliases preserve the archived driver's internal configuration and
% checkpoint schema. The public paths are set only in release_config.
cfg.shpFile = cfg.gisFile;
cfg.stage12MasterFile = cfg.stage12File;
cfg.stage12Dir = fullfile(fileparts(cfg.stage12File),'chunks');

%% ------------------------ Current model parameters ---------------------
turbine.rated_power_MW = 12;
turbine.rotor_diameter_m = 242;
turbine.minimum_spacing_D = 7;
turbine.availability_factor = 0.95;
turbine.annual_degradation_rate = 0.0072;

lb = 118;
lc = 5.77;
dt = 6;

life_farm = 25;
annual_interest_rate = 0.03;
GBP2CNY = 8.41;
USD2CNY = 7.17;
EUR2CNY = 8.39;

platformPar.P_rated = turbine.rated_power_MW;
platformPar.turbine_unit_cost = 2910;
platformPar.steel_price = 4200;
platformPar.GBP2CNY = GBP2CNY;
platformPar.USD2CNY = USD2CNY;
platformPar.EUR2CNY = EUR2CNY;
platformPar.lb = lb;
platformPar.lc = lc;
platformPar.db = turbine.rotor_diameter_m;
platformPar.dt = dt;

%% ------------------------ Load static data -----------------------------
assert(isfile(cfg.shpFile),'Cannot find shapefile: %s',cfg.shpFile);
S = shaperead(cfg.shpFile);

if cfg.rebuildCompactWind || ~isfile(cfg.compactWindFile)
    fprintf('Compact 2023 wind database not found (or rebuild requested). Building once...\n');
    build_2023_compact_wind_db(cfg.shpFile,cfg.speed140File, ...
        cfg.direction140File,cfg.compactWindFile, ...
        turbine.rated_power_MW,turbine.availability_factor);
end

fprintf('Loading compact 2023 wind database fully into RAM...\n');
W2023 = load(cfg.compactWindFile,'wind2023');
windDB = W2023.wind2023;
clear W2023
validateCompactWindDB(windDB,S,turbine);

if ~isfolder(cfg.stage12Dir); mkdir(cfg.stage12Dir); end

firstGrid = max(1,round(cfg.startIndex));
lastGrid = min(numel(S),cfg.endIndex);
if ~isfinite(lastGrid); lastGrid = numel(S); end
lastGrid = floor(lastGrid);
if firstGrid > lastGrid
    error('Invalid grid range: startIndex=%d, endIndex=%d.',firstGrid,lastGrid);
end

fprintf('Stage 1-2 only: grids %d to %d (total feasible grids %d)\n', ...
    firstGrid,lastGrid,numel(S));
fprintf('Checkpoint chunk size = %d; folder: %s\n',cfg.chunkSize,cfg.stage12Dir);

%% ------------------------ Module 1 + 2 batch ---------------------------
chunkStarts = firstGrid:cfg.chunkSize:lastGrid;
for cc = 1:numel(chunkStarts)
    chunkStart = chunkStarts(cc);
    chunkEnd = min(chunkStart + cfg.chunkSize - 1,lastGrid);
    chunkFile = fullfile(cfg.stage12Dir, ...
        sprintf('stage12_2023_grid_%05d_%05d.mat',chunkStart,chunkEnd));

    if cfg.resume && isCompletedStage12Chunk(chunkFile)
        fprintf('\n[Stage12 chunk %d/%d] grids %d-%d already completed; skipped.\n', ...
            cc,numel(chunkStarts),chunkStart,chunkEnd);
        continue;
    end

    fprintf('\n============================================================\n');
    fprintf('[Stage12 chunk %d/%d] grids %d-%d\n',cc,numel(chunkStarts),chunkStart,chunkEnd);
    fprintf('============================================================\n');

    nThis = chunkEnd-chunkStart+1;
    stage12_chunk = repmat(emptyStage12Record(),nThis,1);
    grid_indices = (chunkStart:chunkEnd).';
    chunkTimer = tic;

    for i = chunkStart:chunkEnd
        gridTimer = tic;
        localIdx = i-chunkStart+1;
        fprintf('\n[%d/%d] grid %d: Module 1-2\n',i,lastGrid,i);

        rec = emptyStage12Record();
        rec.grid_index = i;
        rec.data_year = 2023;
        rec.availability_factor = turbine.availability_factor;
        rec.annual_wind_degradation_rate = turbine.annual_degradation_rate;

        try
            shp = S(i);
            boundary_X = shp.X;
            boundary_Y = shp.Y;
            feasible_region_info = [boundary_X(:),boundary_Y(:)];

            % O(1) in-memory wind-column access.
            gridWind = getCompactGridWind(windDB,i);
            rec.wind_data_index = gridWind.source_index;

            % -------- Module 1: wake-aware RK-NSGA turbine micrositing ----
            t1 = tic;
            micro = optimize_turbine_micrositing(feasible_region_info, ...
                gridWind.wind_speed_140m_mps,gridWind.wind_direction_from_deg, ...
                turbine,turbine.availability_factor, ...
                'PopulationSize',80,'MaxGenerations',150,'AmbientTI',0.0555, ...
                'UseParallel',cfg.useParallelMicrositing,'NumWorkers',cfg.numWorkers, ...
                'ParetoSelection','maxAEP','RandomSeed',100);
            rec.module1_time_s = toc(t1);

            compactPowerMismatch = max(abs( ...
                micro.raw_power_pu_2023(:)-gridWind.power_pu(:)));
            if compactPowerMismatch > 1e-9
                error('Compact/micrositing power mismatch: max |delta p.u.| = %.3g.', ...
                    compactPowerMismatch);
            end

            if micro.capacity_MW <= 0 || isempty(micro.turbine_xy_m)
                error('Micrositing returned no buildable turbines.');
            end

            % Save all Module-1 outputs required by later stages/auditing.
            rec.wind_farm_capacity = micro.capacity_MW;
            rec.turbine_count = micro.turbine_count;
            rec.turbine_xy_m = micro.turbine_xy_m;
            rec.wake_adjusted_power_pu_2023 = double(micro.wake_adjusted_power_pu_2023(:));
            rec.AEP_raw_MWh = micro.AEP_raw_MWh;
            rec.AEP_wake_adjusted_MWh = micro.AEP_wake_adjusted_MWh;
            rec.wake_loss_fraction = micro.wake_loss_from_supplied_power_fraction;
            rec.wake_loss_model_fraction = micro.nsga.selected_wake_full.wake_loss_fraction;
            rec.pareto_capacity_MW = micro.nsga.pareto_capacity_MW;
            rec.pareto_wake_loss_fraction = micro.nsga.pareto_wake_loss_fraction;
            rec.pareto_AEP_MWh = micro.nsga.pareto_AEP_wake_MWh_binned;
            rec.selected_pareto_index = micro.nsga.selected_pareto_index;

            % -------- Module 2: platform co-selection ---------------------
            t2 = tic;
            parNow = platformPar;
            parNow.n_T = micro.turbine_count;
            parNow.distance_port = shp.dist2Port;
            parNow.depth = abs(min(shp.RASTERVALU,0));
            [optimal_platform,min_platform_turbine_cost,platformDetail] = ...
                select_optimal_platform(parNow,life_farm,annual_interest_rate);
            rec.module2_time_s = toc(t2);

            rec.optimal_platform = optimal_platform;
            rec.platform_turbine_cost = min_platform_turbine_cost;
            rec.platform_detail = platformDetail;
            rec.status = 'ok';
            rec.error_message = '';

        catch ME
            rec.status = 'error';
            rec.error_message = getReport(ME,'extended','hyperlinks','off');
            fprintf(2,'Grid %d Stage1-2 failed: %s\n',i,ME.message);

            if cfg.stopOnError
                rec.stage12_time_s = toc(gridTimer);
                stage12_chunk(localIdx) = rec;
                stage12_meta = makeStage12ChunkMeta(chunkStart,i,false,toc(chunkTimer),cfg);
                partialFile = strrep(chunkFile,'.mat','_PARTIAL.mat');
                save(partialFile,'stage12_chunk','grid_indices','stage12_meta','cfg','-v7.3');
                rethrow(ME);
            end
        end

        rec.stage12_time_s = toc(gridTimer);
        stage12_chunk(localIdx) = rec;
        fprintf('>>> Grid %d Stage1-2 total time: %.2f min (%.1f s)\n', ...
            i,rec.stage12_time_s/60,rec.stage12_time_s);
    end

    stage12_meta = makeStage12ChunkMeta(chunkStart,chunkEnd,true,toc(chunkTimer),cfg);
    save(chunkFile,'stage12_chunk','grid_indices','stage12_meta','cfg','-v7.3');
    fprintf('\nStage12 chunk saved: %s\n',chunkFile);
    fprintf('Stage12 chunk elapsed time: %.2f h\n',stage12_meta.elapsed_seconds/3600);
    clear stage12_chunk grid_indices
end

%% ------------------------ Merge Stage12 chunks -------------------------
% One final MAT file is convenient for the Stage3-4 continuation driver.
fprintf('\nMerging completed Stage1-2 chunks into:\n  %s\n',cfg.stage12MasterFile);
[stage12_results,stage12_meta] = mergeStage12Chunks( ...
    cfg.stage12Dir,numel(S),firstGrid,lastGrid,cfg);
cfg_stage12 = cfg; %#ok<NASGU>
stage12_grid_indices = (firstGrid:lastGrid).'; %#ok<NASGU>
save(cfg.stage12MasterFile,'stage12_results','stage12_meta', ...
    'stage12_grid_indices','cfg_stage12','-v7.3');

fprintf('\nStage 1-2 completed and merged.\n');
fprintf('Successful grids: %d\n',sum(strcmp({stage12_results(firstGrid:lastGrid).status},'ok')));
fprintf('Master file: %s\n',cfg.stage12MasterFile);

report = stage12_meta;
report.masterFile = cfg.stage12MasterFile;
end

%% ------------------------ local utilities ------------------------------
function G = getCompactGridWind(W,i)
G = struct();
G.power_pu = W.power_pu(:,i);
G.wind_speed_140m_mps = W.wind_speed_140m_mps(:,i);
G.wind_direction_from_deg = W.wind_direction_from_deg(:,i);
G.num_hours = W.num_hours;
if isfield(W,'source_speed_row_col')
    G.source_index = W.source_speed_row_col(i,:);
else
    G.source_index = i;
end
end

function validateCompactWindDB(W,S,turbine)
required = {'power_pu','wind_speed_140m_mps','wind_direction_from_deg', ...
    'grid_count','num_hours','rated_power_MW','availability_factor'};
for k = 1:numel(required)
    if ~isfield(W,required{k})
        error(['Compact wind DB missing field: %s. Set cfg.rebuildCompactWind=true ', ...
            'and rebuild it from 2023 wind speed.'],required{k});
    end
end
if abs(double(W.rated_power_MW)-double(turbine.rated_power_MW)) > 1e-12
    error('Compact DB rated power differs from the current turbine setting. Rebuild compact DB.');
end
if abs(double(W.availability_factor)-double(turbine.availability_factor)) > 1e-12
    error('Compact DB availability differs from the current turbine setting. Rebuild compact DB.');
end
if W.grid_count ~= numel(S)
    error('Compact DB grid count (%d) differs from shapefile (%d).',W.grid_count,numel(S));
end
if size(W.power_pu,2)~=numel(S) || size(W.wind_speed_140m_mps,2)~=numel(S) || ...
        size(W.wind_direction_from_deg,2)~=numel(S)
    error('Compact DB grid dimension does not match shapefile.');
end
if size(W.power_pu,1)~=W.num_hours || size(W.wind_speed_140m_mps,1)~=W.num_hours || ...
        size(W.wind_direction_from_deg,1)~=W.num_hours
    error('Compact DB hourly dimensions are inconsistent.');
end
end

function tf = isCompletedStage12Chunk(fileName)
tf = false;
if ~isfile(fileName); return; end
try
    A = load(fileName,'stage12_meta');
    tf = isfield(A,'stage12_meta') && isfield(A.stage12_meta,'completed') && ...
        logical(A.stage12_meta.completed);
catch
    tf = false;
end
end

function M = makeStage12ChunkMeta(i1,i2,completed,elapsed,cfg)
M = struct();
M.start_grid = i1;
M.end_grid = i2;
M.completed = logical(completed);
M.elapsed_seconds = elapsed;
M.created_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
M.data_year = 2023;
M.chunk_size = cfg.chunkSize;
M.compact_wind_file = cfg.compactWindFile;
M.modules = 'Module 1 + Module 2 only';
end

function r = emptyStage12Record()
r = struct( ...
    'grid_index',NaN,'status','','data_year',2023, ...
    'availability_factor',NaN,'annual_wind_degradation_rate',NaN, ...
    'wind_data_index',[], ...
    'wind_farm_capacity',NaN,'turbine_count',NaN,'turbine_xy_m',[], ...
    'wake_adjusted_power_pu_2023',[], ...
    'AEP_raw_MWh',NaN,'AEP_wake_adjusted_MWh',NaN, ...
    'wake_loss_fraction',NaN,'wake_loss_model_fraction',NaN, ...
    'pareto_capacity_MW',[],'pareto_wake_loss_fraction',[], ...
    'pareto_AEP_MWh',[],'selected_pareto_index',NaN, ...
    'optimal_platform','','platform_turbine_cost',NaN,'platform_detail',struct(), ...
    'module1_time_s',NaN,'module2_time_s',NaN,'stage12_time_s',NaN, ...
    'error_message','');
end

function [R,M] = mergeStage12Chunks(stage12Dir,nGrid,firstGrid,lastGrid,cfg)
R = repmat(emptyStage12Record(),nGrid,1);
files = dir(fullfile(stage12Dir,'stage12_2023_grid_*.mat'));
if isempty(files)
    error('No Stage1-2 chunk files found in %s.',stage12Dir);
end

loaded = false(nGrid,1);
for k = 1:numel(files)
    fn = fullfile(files(k).folder,files(k).name);
    A = load(fn,'stage12_chunk','grid_indices','stage12_meta');
    if ~isfield(A,'stage12_meta') || ~logical(A.stage12_meta.completed)
        continue;
    end
    idx = A.grid_indices(:);
    if numel(idx) ~= numel(A.stage12_chunk)
        error('Stage12 chunk index/record length mismatch: %s',fn);
    end
    for q = 1:numel(idx)
        ii = idx(q);
        if ii>=1 && ii<=nGrid
            R(ii) = A.stage12_chunk(q);
            loaded(ii) = true;
        end
    end
end

requested = (firstGrid:lastGrid).';
if any(~loaded(requested))
    miss = requested(~loaded(requested));
    error('Stage1-2 merge is incomplete. Missing %d grid records; first missing grid = %d.', ...
        numel(miss),miss(1));
end

M = struct();
M.completed = true;
M.start_grid = firstGrid;
M.end_grid = lastGrid;
M.total_shapefile_grids = nGrid;
M.num_requested_grids = numel(requested);
M.num_success = sum(strcmp({R(requested).status},'ok'));
M.num_error = sum(strcmp({R(requested).status},'error'));
M.created_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
M.data_year = 2023;
M.chunk_size = cfg.chunkSize;
M.stage12_chunk_dir = stage12Dir;
M.compact_wind_file = cfg.compactWindFile;
end
