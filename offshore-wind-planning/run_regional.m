function reports = run_regional(cfg)
%RUN_REGIONAL Capacity-constrained selection and cross-grid spacing repair.
% Run run_stage1 and run_stage34 first using the same release_config.
% Creates private regional inputs from the newly generated B0 checkpoints.
if nargin<1 || isempty(cfg), cfg=release_config(); end
required={'rootDir','baselineFile','compactWindFile','priceFile', ...
    'regionalOutputDir','targetCapacitiesGW','solverThreads', ...
    'storageTimeLimitSeconds','resume','numWorkers'};
assert(isstruct(cfg) && isscalar(cfg) && all(isfield(cfg,required)), ...
    'run_regional:Config','Use a complete release_config structure.');
validateattributes(cfg.targetCapacitiesGW,{'numeric'},{'vector','finite','positive','nonempty'});
oldPath=path; restorePath=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(cfg.rootDir,'src','regional'),'-begin');
settings=prepare_regional_inputs(cfg);
settings.TargetCapacityMW=double(cfg.targetCapacitiesGW(:).')*1000;
settings.ResultDir=fullfile(cfg.regionalOutputDir,'results');
settings.CacheDir=fullfile(cfg.regionalOutputDir,'cache');
settings.SolverThreads=cfg.solverThreads;
settings.StorageTimeLimitSeconds=cfg.storageTimeLimitSeconds;
settings.Resume=cfg.resume;
settings.RecomputeWorkers=1;
if ~isempty(cfg.numWorkers), settings.RecomputeWorkers=cfg.numWorkers; end
reports=p44_run(settings);
end
