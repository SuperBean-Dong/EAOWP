function cfg = release_config(overrides)
%RELEASE_CONFIG Portable paths and execution settings for the main model.
%   cfg = release_config();
%   cfg = release_config(struct('dataDir','/path/to/private/data', ...
%                               'outputDir','/path/to/new/outputs'));
%   Input data and generated outputs are deliberately not distributed.
if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('release:BadConfig','overrides must be a scalar structure.');
end
cfg.rootDir = fileparts(mfilename('fullpath'));
cfg.dataDir = fullfile(cfg.rootDir,'data');
cfg.outputDir = fullfile(cfg.rootDir,'outputs');
cfg.startIndex = 1;
cfg.endIndex = Inf;
cfg.gridIDs = [];
cfg.chunkSize = 50;
cfg.resume = true;
cfg.stopOnError = false;
cfg.useParallelMicrositing = false;
cfg.numWorkers = [];
cfg.rebuildCompactWind = false;
cfg.collectionSolver = 'fast';
cfg.ossBlockTimeLimitSeconds = 60;
cfg.ossGridTimeLimitSeconds = 600;
cfg.ossMIPGap = 0.05;
cfg.ossEnablePruning = true;
cfg.solverThreads = 2;
cfg.storageTimeLimitSeconds = 300;
cfg.storageGridTimeLimitSeconds = 1800;
cfg.storageMaxPasses = 10;
cfg.saveStage34Every = 1;
cfg.targetCapacitiesGW = 100:100:1500;

% Resolve directory overrides before calculating dependent filenames.
for name = {'dataDir','outputDir'}
    key = name{1};
    if isfield(overrides,key), cfg.(key) = overrides.(key); end
end
cfg.gisFile = fullfile(cfg.dataDir,'spatial','feasible_final.shp');
cfg.speed140File = fullfile(cfg.dataDir,'weather','south_sea_wind_speed_140m_2023.mat');
cfg.direction140File = fullfile(cfg.dataDir,'weather','south_sea_wind_direction_140m_2023.mat');
cfg.priceFile = fullfile(cfg.dataDir,'economics','electricity_price.mat');
cfg.compactWindFile = fullfile(cfg.outputDir,'cache','wind_2023_feasible_compact.mat');
cfg.stage12File = fullfile(cfg.outputDir,'stage12','stage12_2023_7D_all.mat');
cfg.baselineFile = fullfile(cfg.outputDir,'b0','baseline_results_2023_7D_B0.mat');
cfg.regionalOutputDir = fullfile(cfg.outputDir,'regional');
names = fieldnames(overrides);
for k = 1:numel(names)
    key = names{k};
    if ~isfield(cfg,key)
        error('release:UnknownConfig','Unknown configuration field: %s',key);
    end
    if strcmp(key,'rootDir')
        error('release:RootOverride','rootDir is determined from the code location.');
    end
    cfg.(key) = overrides.(key);
end
if ~ismember(string(cfg.collectionSolver),["fast","solver_init"])
    error('release:BadCollectionSolver','collectionSolver must be fast or solver_init.');
end
end
