function setup_paths()
%SETUP_PATHS Add only the repository's own implementation folders.
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir,'src','stage1','functions'));
addpath(fullfile(rootDir,'src','stage34'));
addpath(fullfile(rootDir,'src','regional'));
% Install/configure YALMIP and the licensed Gurobi MATLAB interface separately.
end
