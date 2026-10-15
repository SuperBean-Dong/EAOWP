function out = bastankhah_gaussian_wake(layout_xy_m, wind_speed_mps, wind_dir_from_deg, turbine, varargin)
%BASTANKHAH_GAUSSIAN_WAKE  Vectorized Gaussian engineering wake model.
%
% 加速版本：
%   1) 相同风向的所有风速状态成组计算；
%   2) downwind/crosswind几何每个风向只计算一次；
%   3) 同一风向下的全部风速状态向量化；
%   4) 支持由外层GA传入预计算的两两几何矩阵。
%
% Wake deficit: Bastankhah & Porte-Agel Gaussian model.
% Wake expansion: k* = 0.3837 I + 0.003678, default I = 0.0555.
% Multiple wakes: RSS of absolute velocity deficits.
%
% INPUTS
% layout_xy_m          N x 2 [xEast,yNorth] [m]
% wind_speed_mps       T x 1 free-stream hub-height wind speed [m/s]
% wind_dir_from_deg    T x 1 meteorological direction [deg FROM north]
% turbine              struct: rotor_diameter_m, rated_power_MW and Ct/Power
%
% OPTIONS
% 'AmbientTI'           scalar or T x 1, default 0.0555
% 'WakeExpansionK'      optional scalar/T x 1; overrides k*(I)
% 'StateWeightsHours'   scalar or T x 1 represented hours, default 1
% 'RotorSamples'        rotor disk sample count, default 1
% 'StoreTurbineStates'  store T x N matrices, default true
% 'PrecomputedGeometry' optional struct:
%                       delta_x_m(i,j)=x_i-x_j
%                       delta_y_m(i,j)=y_i-y_j
% 'DirectionToleranceDeg' grouping tolerance, default 1e-9 deg

p = inputParser;
p.addRequired('layout_xy_m', @(x) isnumeric(x) && size(x,2)==2 && all(isfinite(x(:))));
p.addRequired('wind_speed_mps', @(x) isnumeric(x) && isvector(x));
p.addRequired('wind_dir_from_deg', @(x) isnumeric(x) && isvector(x));
p.addRequired('turbine', @isstruct);
p.addParameter('AmbientTI', 0.0555, @(x) isnumeric(x) && isvector(x));
p.addParameter('WakeExpansionK', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('StateWeightsHours', 1, @(x) isnumeric(x) && isvector(x));
p.addParameter('RotorSamples', 1, @(x) isnumeric(x) && isscalar(x) && x>=1 && isfinite(x));
p.addParameter('StoreTurbineStates', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('PrecomputedGeometry', [], @(x) isempty(x) || isstruct(x));
p.addParameter('DirectionToleranceDeg', 1e-9, @(x) isnumeric(x) && isscalar(x) && x>0);
p.parse(layout_xy_m, wind_speed_mps, wind_dir_from_deg, turbine, varargin{:});
opt = p.Results;

xy = double(layout_xy_m);
U0 = double(wind_speed_mps(:));
WD = mod(double(wind_dir_from_deg(:)),360);
T = numel(U0);
N = size(xy,1);
if numel(WD) ~= T
    error('wind_speed_mps and wind_dir_from_deg must have the same length.');
end
if N < 1
    error('layout_xy_m must contain at least one turbine.');
end
if any(~isfinite(U0)) || any(U0<0) || any(~isfinite(WD))
    error('Wind inputs must be finite; wind speed must be >=0.');
end

if ~isfield(turbine,'rotor_diameter_m') || ~isfield(turbine,'rated_power_MW')
    error('turbine.rotor_diameter_m and turbine.rated_power_MW are required.');
end
D = double(turbine.rotor_diameter_m);
Pr = double(turbine.rated_power_MW);
CtFcn = makeCtEvaluator(turbine);
PowerFcn = makePowerEvaluator(turbine);

TI = expandToLength(opt.AmbientTI,T,'AmbientTI');
weights_h = expandToLength(opt.StateWeightsHours,T,'StateWeightsHours');
if any(~isfinite(TI)) || any(TI<0)
    error('AmbientTI must be finite and nonnegative.');
end
if any(~isfinite(weights_h)) || any(weights_h<0)
    error('StateWeightsHours must be finite and nonnegative.');
end

kstarFitRange = [0.065,0.15];
if isempty(opt.WakeExpansionK)
    kstar = 0.3837.*TI + 0.003678;
    kstarFromTI = true;
    kstarRelationExtrapolated = TI<kstarFitRange(1) | TI>kstarFitRange(2);
else
    kstar = expandToLength(opt.WakeExpansionK,T,'WakeExpansionK');
    kstarFromTI = false;
    kstarRelationExtrapolated = false(T,1);
end
if any(~isfinite(kstar)) || any(kstar<=0)
    error('Wake expansion coefficient k* must be finite and >0.');
end

nRotor = max(1,round(opt.RotorSamples));
[yRotor_m,zRotor_m] = makeRotorSamples(nRotor,D);
storeStates = logical(opt.StoreTurbineStates);

% Pairwise geometry. Outer GA may supply this to avoid rebuilding it.
usedPrecomputedGeometry = false;
if ~isempty(opt.PrecomputedGeometry)
    G = opt.PrecomputedGeometry;
    if ~isfield(G,'delta_x_m') || ~isfield(G,'delta_y_m')
        error('PrecomputedGeometry requires delta_x_m and delta_y_m.');
    end
    deltaX = double(G.delta_x_m);
    deltaY = double(G.delta_y_m);
    if ~isequal(size(deltaX),[N N]) || ~isequal(size(deltaY),[N N])
        error('PrecomputedGeometry matrices must be N x N.');
    end
    usedPrecomputedGeometry = true;
else
    deltaX = xy(:,1)-xy(:,1).';
    deltaY = xy(:,2)-xy(:,2).';
end

farmP_free = zeros(T,1);
farmP_wake = zeros(T,1);
stateLoss = nan(T,1);
if storeStates
    Ueff_all = zeros(T,N);
    Ct_all = zeros(T,N);
    P_all = zeros(T,N);
else
    Ueff_all = [];
    Ct_all = [];
    P_all = [];
end

% Quantize only floating point noise, then group states by identical direction.
qTol = opt.DirectionToleranceDeg;
WDq = mod(round(WD./qTol).*qTol,360);
[uniqueDir,~,dirGroup] = unique(WDq,'sorted');
nUniqueDir = numel(uniqueDir);

for dd = 1:nUniqueDir
    stateIdx = find(dirGroup==dd);
    Uinf = U0(stateIdx);       % nS x 1
    kVec = kstar(stateIdx);    % nS x 1
    nS = numel(stateIdx);

    dirDeg = uniqueDir(dd);
    eFlow = [-sind(dirDeg),-cosd(dirDeg)];
    eCross = [-eFlow(2),eFlow(1)];

    % Geometry is evaluated once for this direction, not once per speed state.
    xDownMat = deltaX.*eFlow(1) + deltaY.*eFlow(2);
    yCrossMat = deltaX.*eCross(1) + deltaY.*eCross(2);
    s = xy*eFlow(:);
    [~,order] = sort(s,'ascend');

    % All wind-speed states of this direction are processed simultaneously.
    Ueff = zeros(nS,N);
    CtEff = zeros(nS,N);
    Pwt = zeros(nS,N);

    for pos = 1:N
        i = order(pos);
        if pos==1
            Ui = Uinf;
        else
            deficitSq = zeros(nS,nRotor);
            for uu = 1:pos-1
                j = order(uu);
                xDown = xDownMat(i,j);
                if xDown<=0
                    continue;
                end

                Ctj = CtEff(:,j);
                Uj = Ueff(:,j);
                active = isfinite(Ctj) & Ctj>0 & Uj>0 & Uinf>0;
                if ~any(active)
                    continue;
                end

                CtSafe = min(max(Ctj,0),0.999999);
                root1mCt = sqrt(max(1-CtSafe,eps));
                beta = (1+root1mCt)./(2.*root1mCt);
                epsilon = 0.2.*sqrt(beta);
                sigmaD = kVec.*(xDown/D) + epsilon;

                sqrtArg = max(1-CtSafe./(8.*sigmaD.^2),eps);
                centerAmp = 1-sqrt(sqrtArg);

                yHub = yCrossMat(i,j);
                yDim = (yHub+yRotor_m)./D;
                zDim = zRotor_m./D;
                rotorR2 = yDim.^2 + zDim.^2;
                deltaFrac = centerAmp .* exp(-0.5.*rotorR2./sigmaD.^2);
                deltaFrac = min(max(deltaFrac,0),0.999999);
                deltaFrac(~active,:) = 0;

                du = deltaFrac.*Uj;
                deficitSq = deficitSq + du.^2;
            end
            Urotor = max(0,Uinf-sqrt(deficitSq));
            Ui = mean(Urotor,2);
        end

        Ui = max(0,Ui);
        Ueff(:,i) = Ui;
        Ct_i = CtFcn(Ui);
        P_i = PowerFcn(Ui);
        Ct_i(~isfinite(Ct_i)) = 0;
        P_i(~isfinite(P_i)) = 0;
        CtEff(:,i) = max(0,Ct_i);
        Pwt(:,i) = min(max(0,P_i),Pr);
    end

    PfreeOne = PowerFcn(Uinf);
    PfreeOne(~isfinite(PfreeOne)) = 0;
    PfreeOne = min(max(0,PfreeOne),Pr);
    farmP_free(stateIdx) = N.*PfreeOne;
    farmP_wake(stateIdx) = sum(Pwt,2);

    validFree = farmP_free(stateIdx)>0;
    tmpLoss = nan(nS,1);
    tmpLoss(validFree) = 1-farmP_wake(stateIdx(validFree))./farmP_free(stateIdx(validFree));
    stateLoss(stateIdx) = tmpLoss;

    if storeStates
        Ueff_all(stateIdx,:) = Ueff;
        Ct_all(stateIdx,:) = CtEff;
        P_all(stateIdx,:) = Pwt;
    end
end

AEP_free = sum(farmP_free.*weights_h);
AEP_wake = sum(farmP_wake.*weights_h);
representedHours = sum(weights_h);
installedCapacity = N*Pr;
if AEP_free>0
    wakeEfficiency = AEP_wake/AEP_free;
    wakeLoss = 1-wakeEfficiency;
else
    wakeEfficiency = NaN;
    wakeLoss = NaN;
end
if representedHours>0 && installedCapacity>0
    CFwake = AEP_wake/(installedCapacity*representedHours);
else
    CFwake = NaN;
end

if N>=2
    dist2 = deltaX.^2+deltaY.^2;
    dist2(1:N+1:end) = inf;
    minPair = sqrt(min(dist2(:)));
else
    minPair = inf;
end

out = struct();
out.model = 'Bastankhah-Porté-Agel Gaussian wake + Niayifar k*(TI) + RSS superposition';
out.wind_direction_convention = 'meteorological FROM north clockwise';
out.num_states = T;
out.num_unique_directions = nUniqueDir;
out.num_turbines = N;
out.rotor_diameter_m = D;
out.rated_power_MW_per_turbine = Pr;
out.installed_capacity_MW = installedCapacity;
out.represented_hours = representedHours;
out.rotor_samples = nRotor;
out.state_weights_hours = weights_h;
out.kstar = kstar;
out.ambient_TI = TI;
out.kstar_from_TI_relation = kstarFromTI;
out.kstar_relation_fit_range_TI = kstarFitRange;
out.kstar_relation_extrapolated = kstarRelationExtrapolated;
out.kstar_baseline_TI = 0.0555;
out.kstar_baseline_value = 0.3837*0.0555+0.003678;
out.farm_power_free_MW = farmP_free;
out.farm_power_wake_MW = farmP_wake;
out.state_wake_loss_fraction = stateLoss;
out.AEP_free_MWh = AEP_free;
out.AEP_wake_MWh = AEP_wake;
out.wake_efficiency = wakeEfficiency;
out.wake_loss_fraction = wakeLoss;
out.capacity_factor_wake = CFwake;
out.min_pair_distance_m = minPair;
out.min_pair_distance_D = minPair/D;
out.used_precomputed_geometry = usedPrecomputedGeometry;
out.vectorized_by_direction = true;
if storeStates
    out.Ueff_mps = Ueff_all;
    out.Ct_effective = Ct_all;
    out.turbine_power_MW = P_all;
end
end

function v = expandToLength(x,n,name)
x = double(x(:));
if isscalar(x)
    v = repmat(x,n,1);
elseif numel(x)==n
    v = x;
else
    error('%s must be scalar or match wind_speed_mps length.',name);
end
end

function f = makeCtEvaluator(turbine)
if isfield(turbine,'CtFcn') && isa(turbine.CtFcn,'function_handle')
    userFcn = turbine.CtFcn;
    f = @(U) evalVectorFunction(userFcn,U,'CtFcn');
elseif isfield(turbine,'Ct_U_mps') && isfield(turbine,'Ct_values')
    x = double(turbine.Ct_U_mps(:)); y = double(turbine.Ct_values(:));
    [x,idx] = sort(x); y = y(idx);
    f = @(U) interp1(x,y,U,'linear',0);
else
    error('Provide turbine.CtFcn or Ct_U_mps/Ct_values.');
end
end

function f = makePowerEvaluator(turbine)
if isfield(turbine,'PowerFcn') && isa(turbine.PowerFcn,'function_handle')
    userFcn = turbine.PowerFcn;
    f = @(U) evalVectorFunction(userFcn,U,'PowerFcn');
elseif isfield(turbine,'Power_U_mps') && isfield(turbine,'Power_MW')
    x = double(turbine.Power_U_mps(:)); y = double(turbine.Power_MW(:));
    [x,idx] = sort(x); y = y(idx);
    f = @(U) interp1(x,y,U,'linear',0);
else
    error('Provide turbine.PowerFcn or Power_U_mps/Power_MW.');
end
end

function y = evalVectorFunction(fcn,x,name)
try
    y = fcn(x);
    if isscalar(y) && ~isscalar(x)
        y = repmat(y,size(x));
    end
    if numel(y)~=numel(x)
        error('size mismatch');
    end
    y = reshape(double(y),size(x));
catch
    try
        y = arrayfun(@(z) double(fcn(z)),x);
    catch ME
        error('%s evaluation failed: %s',name,ME.message);
    end
end
end

function [ys,zs] = makeRotorSamples(n,D)
if n<=1
    ys = 0; zs = 0; return;
end
q = (1:n).';
goldenAngle = pi*(3-sqrt(5));
r = 0.5*D.*sqrt((q-0.5)./n);
a = q.*goldenAngle;
ys = (r.*cos(a)).';
zs = (r.*sin(a)).';
end
