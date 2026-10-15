function report = check_release(cfg)
%CHECK_RELEASE Read-only dependency/input checks. Does not run optimization.
if nargin < 1 || isempty(cfg), cfg = release_config(); end
setup_paths();
report = struct();
report.matlabVersion = version;
report.dependencies.shaperead = ~isempty(which('shaperead'));
report.dependencies.gurobi = ~isempty(which('gurobi'));
report.dependencies.yalmip = ~isempty(which('sdpvar')) && ~isempty(which('optimize'));
report.dependencies.linprog = ~isempty(which('linprog'));
report.dependencies.parallelRequested = logical(cfg.useParallelMicrositing);
report.dependencies.parallelAvailable = ~isempty(which('parpool'));
report.inputs.gis = isfile(cfg.gisFile);
report.inputs.shapefileIndex = isfile(replaceExtension(cfg.gisFile,'.shx'));
report.inputs.shapefileAttributes = isfile(replaceExtension(cfg.gisFile,'.dbf'));
report.inputs.shapefileProjection = isfile(replaceExtension(cfg.gisFile,'.prj'));
report.inputs.windSpeed = isfile(cfg.speed140File);
report.inputs.windDirection = isfile(cfg.direction140File);
report.inputs.compactWind = isfile(cfg.compactWindFile);
report.inputs.electricityPrice = isfile(cfg.priceFile);
report.generated.stage12 = isfile(cfg.stage12File);
report.generated.baseline = isfile(cfg.baselineFile);
gridCheckpoints = dir(fullfile(fileparts(cfg.baselineFile),'checkpoints','grid_*.mat'));
report.generated.stage34Checkpoints = ~isempty(gridCheckpoints);
spatialComplete = report.inputs.gis && report.inputs.shapefileIndex && ...
    report.inputs.shapefileAttributes && report.inputs.shapefileProjection;
usableCompact = report.inputs.compactWind && ~cfg.rebuildCompactWind;
report.ready.stage1 = report.dependencies.shaperead && spatialComplete && ...
    (usableCompact || (report.inputs.windSpeed && report.inputs.windDirection));
report.ready.stage34 = report.dependencies.shaperead && report.dependencies.gurobi && ...
    report.dependencies.yalmip && spatialComplete && report.inputs.electricityPrice && report.generated.stage12;
report.ready.regional = report.dependencies.gurobi && report.inputs.electricityPrice && ...
    report.generated.baseline && report.inputs.compactWind && report.generated.stage34Checkpoints;
fprintf('MATLAB: %s\n', report.matlabVersion);
fprintf('Mapping reader: %d; Gurobi interface: %d; YALMIP: %d\n', ...
    report.dependencies.shaperead,report.dependencies.gurobi,report.dependencies.yalmip);
fprintf('Optional linprog: %d; parallel requested/available: %d/%d\n', ...
    report.dependencies.linprog,report.dependencies.parallelRequested,report.dependencies.parallelAvailable);
fprintf('Inputs present: GIS=%d, wind=%d, compact=%d, price=%d\n', ...
    report.inputs.gis,report.inputs.windSpeed && report.inputs.windDirection, ...
    report.inputs.compactWind,report.inputs.electricityPrice);
fprintf('Generated Stage12=%d; baseline=%d\n',report.generated.stage12,report.generated.baseline);
fprintf('Read-only readiness: Stage1=%d, Stage34=%d, regional=%d\n', ...
    report.ready.stage1,report.ready.stage34,report.ready.regional);
fprintf('Interface detection does not confirm a valid solver license.\n');
fprintf('Missing data in a clean clone is expected; see data/README.md.\n');
end

function file = replaceExtension(file,extension)
[folder,name] = fileparts(file);
file = fullfile(folder,[name extension]);
end
