function cfg=p44_config(overrides)
% Confirmed scientific settings; runtime/path settings remain configurable.
if nargin<1; overrides=struct(); end
sourceDir=fileparts(mfilename('fullpath'));
package=fileparts(fileparts(sourceDir));
regionalOutput=fullfile(package,'outputs','regional');
inputRoot=fullfile(regionalOutput,'input_cache');
cfg.PackageDir=package;
cfg.SourceCodeDir=sourceDir;
cfg.ManifestFile=fullfile(inputRoot,'data','input_manifest.json');
cfg.BaselineDir=fullfile(inputRoot,'baseline_snapshot');
cfg.InputFile=fullfile(inputRoot,'data','p44_inputs.mat');
cfg.WindFile=fullfile(package,'outputs','cache','wind_2023_feasible_compact.mat');
cfg.PriceFile=fullfile(inputRoot,'data','p44_price.mat');
cfg.ResultDir=fullfile(regionalOutput,'results');
cfg.CacheDir=fullfile(regionalOutput,'cache');
cfg.TargetCapacityMW=[]; cfg.CapacityToleranceMW=552;
cfg.MinimumSpacingM=7*242; cfg.SpacingNumericalToleranceM=1e-6;
cfg.SolverThreads=2; cfg.SelectionTimeLimitSeconds=60; cfg.SelectionMIPGap=.001;
cfg.StorageTimeLimitSeconds=300; cfg.RecomputeWorkers=1;
cfg.MaxIterations=4000; cfg.Resume=true; cfg.EvaluateFailedPortfolio=false;
cfg.ZeroGridPolicy='skip_and_continue';
cfg.LifeFarmYears=25; cfg.AnnualInterestRate=.03;
cfg.AnnualWindDegradationRate=.0072; cfg.EUR2CNY=8.39;
cfg.StorageParameters=struct('eta_ch',.954,'eta_dis',.954,'self_discharge_per_h',.00054, ...
    'soc_min',.1,'soc_max',1,'duration_h',2,'energy_capex_CNY_per_kWh',1043.82, ...
    'platform_capex_EUR_per_kWh',0,'battery_life_years',15,'om_fraction_per_year',.04,'EUR2CNY',8.39);
cfg.ambient_TI=.0555; cfg.wind_speed_bin_width_mps=1;
cfg.wind_direction_bin_width_deg=15; cfg.rotor_samples=9;
cfg.turbine=struct('rated_power_MW',12,'rotor_diameter_m',242,'availability_factor',.95);
fields=fieldnames(overrides);
allowed={'ManifestFile','SourceCodeDir','BaselineDir','InputFile','WindFile','PriceFile','ResultDir','CacheDir', ...
    'TargetCapacityMW','CapacityToleranceMW','SolverThreads','SelectionTimeLimitSeconds', ...
    'SelectionMIPGap','StorageTimeLimitSeconds','MaxIterations','Resume','RecomputeWorkers', ...
    'EvaluateFailedPortfolio'};
for k=1:numel(fields)
    assert(ismember(fields{k},allowed),'Unknown/frozen setting %s.',fields{k});
    cfg.(fields{k})=overrides.(fields{k});
end
assert(~isempty(cfg.TargetCapacityMW) && all(isfinite(cfg.TargetCapacityMW))&&all(cfg.TargetCapacityMW>0));
assert(isscalar(cfg.CapacityToleranceMW)&&isfinite(cfg.CapacityToleranceMW)&&cfg.CapacityToleranceMW>=0);
assert(cfg.RecomputeWorkers>=1&&cfg.RecomputeWorkers==floor(cfg.RecomputeWorkers));
assert(cfg.SolverThreads>=1&&cfg.SolverThreads==floor(cfg.SolverThreads));
end
