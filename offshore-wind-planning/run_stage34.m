function report = run_stage34(cfg)
%RUN_STAGE34 Compute collection design and current baseline delivery design.
%   REPORT = RUN_STAGE34() uses RELEASE_CONFIG().
%   Stage-I records must be supplied in cfg.stage12File. No archived baseline,
%   machine assignment, historical result, or bundled research input is used.
%   See docs/STAGE34.md for required inputs and the checkpoint/result schema.
if nargin < 1 || isempty(cfg)
    cfg = release_config();
end
rootDir = fileparts(mfilename('fullpath'));
previousPath = path;
pathGuard = onCleanup(@()path(previousPath));
addpath(fullfile(rootDir,'src','stage34'),'-begin');
report = release_stage34_run(cfg);
end
