function [best_capacity_value, wind_tubine_point, result] = ga_micrositing_max_capacity(feasible_region_info, turbine_spacing, varargin)
%GA_MICROSITING_MAX_CAPACITY 海上风机微观选址：兼容旧GA，并支持尾流感知NSGA-II。
%
% 该函数尽量保持原主程序接口不变：
%   [best_capacity_value, wind_tubine_point, result] = ...
%       ga_micrositing_max_capacity(feasible_region_info, turbine_spacing, ...)
%
% 当未提供风速/风向时：
%   - 自动进入 legacy-capacity 模式；
%   - 与原代码一致，以“最大风机数量/装机容量”为单目标；
%   - 约束任意两台风机距离 >= turbine_spacing。
%
% 当提供 WindSpeedMps + WindDirFromDeg 时：
%   - 自动进入 wake-nsga2 模式；
%   - 使用 random-key + NSGA-II；
%   - 目标1：最大化装机容量（内部写成最小化 -capacity）；
%   - 目标2：最小化年度尾流损失率；
%   - 尾流损失由 bastankhah_gaussian_wake.m 计算；
%   - C_T(U) 由 GWH242_12MW_Ct_proxy_IEA15.m 提供；
%   - P(U) 由 GWH242_12MW_power_curve.m 提供；可利用率由 AvailabilityFactor 显式传入。
%
% -------------------------------------------------------------------------
% 一、核心空间设定
% -------------------------------------------------------------------------
% turbine_spacing 为硬约束的最小两机距离。当前研究建议：
%   D = 242 m
%   d_min = 7D = 1694 m
%
% “7D-11D”在本函数中不解释为“所有两两距离均<=11D”。11D用于构造一个
% 较低开发密度参考：先估计同一真实可行域在11D最小间距下能够容纳的风机数
% N_11D，再让NSGA-II的目标风机数在 N_11D 与7D高密度上限附近变化。
% 因此个体内部任意两台风机可以相距大于11D；真正的硬约束只有 d_ij >= 7D。
%
% -------------------------------------------------------------------------
% 二、风况输入
% -------------------------------------------------------------------------
% 推荐主程序先由ERA5 140 m的U/V分量计算：
%   [wind_speed, wind_dir_from] = era5_uv_to_wind_speed_direction(u140, v140);
%
% 再传入：
%   'WindSpeedMps', wind_speed, ...
%   'WindDirFromDeg', wind_dir_from, ...
%
% wind_dir_from_deg 采用气象学“FROM”方向：
%   0°=北风，90°=东风，180°=南风，270°=西风。
%
% -------------------------------------------------------------------------
% 三、主要Name-Value参数
% -------------------------------------------------------------------------
% 通用/旧接口：
%   'TurbinePowerMW'     单机额定容量，默认12 MW
%   'PopulationSize'     种群规模，默认80
%   'MaxGenerations'     最大代数，默认150
%   'EliteRatio'         仅legacy模式使用，默认0.05
%   'CrossoverRate'      交叉概率，默认0.85
%   'MutationStart'      初始变异率，默认0.08
%   'MutationEnd'        末期变异率，默认0.02
%   'MutationBoost'      停滞时变异率提升值，默认0.12
%   'StallTrigger'       停滞代数阈值，默认15
%   'NumAngles'          六边形候选库旋转角数量，默认4
%   'PlotFlag'           是否绘图，默认true
%   'RandomSeed'         随机数种子，默认[]
%
% 尾流NSGA-II：
%   'Mode'               'auto'|'legacy'|'wake-nsga2'，默认'auto'
%   'WindSpeedMps'       T×1轮毂高度自由来流风速[m/s]
%   'WindDirFromDeg'     T×1气象学来流方向[deg]
%   'RotorDiameterM'     默认242 m
%   'AvailabilityFactor' 风机可利用率；主流程由外部显式传入（函数默认1.0）
%   'AmbientTI'          默认0.0555
%   'DensityReferenceMaxSpacingD' 默认11；用于估计N_11D
%   'UseBinnedWindStates' 默认true；NSGA内部使用风况分箱加速
%   'WindSpeedBinWidthMps' 默认1 m/s
%   'WindDirectionBinWidthDeg' 默认15 deg（典型8760 h约数百个有效状态）
%   'UseParallel'        默认true；并行评价未缓存的唯一布局
%   'PrecomputePairGeometry' 默认true；候选点两两几何仅计算一次
%   'NSGARotorSamples'   NSGA内部转子采样数，默认1
%   'FinalRotorSamples'  最终代表风况复算转子采样数，默认9
%   'FinalUseRepresentativeStates' 默认true；不逐时执行8760次wake传播
%   'FinalWindSpeedBinWidthMps' 默认0.5 m/s
%   'FinalWindDirectionBinWidthDeg' 默认10 deg
%   'FinalStoreTurbineStates' 仅在关闭代表风况重构时有效
%   'ParetoSelection'    'maxAEP'|'knee'|'maxCapacity'|'minWake'
%                        默认'maxAEP'，用于在保持旧主程序单一布局输出接口时
%                        从Pareto前沿选择一个布局。
%
% 重要说明：完整Pareto前沿仍保存在 result.pareto_* 中。后续如果希望由Stage II
% NPV在多个Pareto布局中最终选择，可直接遍历 result.pareto_layouts，而不需要
% 修改本函数的NSGA求解过程。
%
% -------------------------------------------------------------------------
% 四、输出
% -------------------------------------------------------------------------
% best_capacity_value : 兼容旧主程序的“被选中布局”装机容量[MW]
% wind_tubine_point   : 被选中布局的N×2坐标[m]
% result              : 详细结果；wake-nsga2模式额外包含：
%   .pareto_capacity_MW
%   .pareto_wake_loss_fraction
%   .pareto_AEP_wake_MWh_binned
%   .pareto_layouts
%   .selected_pareto_index
%   .selected_wake_full
%   .wind_power_pu_wake_full  % 可直接替换旧主程序grid_wind_power_pu
%   .wind_power_pu_free_full
%   .N_reference_7D
%   .N_reference_11D
%
% -------------------------------------------------------------------------
% 五、参考模型
% -------------------------------------------------------------------------
% Wake deficit: Bastankhah & Porte-Agel Gaussian model.
% Wake expansion: k* = 0.3837 I + 0.003678, baseline I=0.0555.
% Multi-wake superposition: RSS velocity deficits.
%
% 本函数中的NSGA-II为自包含实现，不依赖gamultiobj；保留原random-key编码、
% 真实非规则/带孔洞可行域和冲突矩阵解码逻辑。
%
% 当前加速版本同时采用：
%   (1) 8760 h -> U-WD代表风况；
%   (2) parfor并行评价唯一布局；
%   (3) 批内布局去重 + 跨代cache；
%   (4) 候选点两两几何预计算；
%   (5) wake函数按风向分组，并向量化同方向下全部风速状态。

% ============================== 参数解析 ================================
p = inputParser;
p.addRequired('feasible_region_info');
p.addRequired('turbine_spacing', @(x) isnumeric(x) && isscalar(x) && x > 0);

% 原函数参数：保持兼容
p.addParameter('TurbinePowerMW', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('PopulationSize', 80, @(x) isnumeric(x) && isscalar(x) && x >= 20);
p.addParameter('MaxGenerations', 150, @(x) isnumeric(x) && isscalar(x) && x >= 5);
p.addParameter('EliteRatio', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 0.5);
p.addParameter('CrossoverRate', 0.85, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('MutationStart', 0.08, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('MutationEnd', 0.02, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('MutationBoost', 0.12, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('StallTrigger', 15, @(x) isnumeric(x) && isscalar(x) && x >= 3);
p.addParameter('NumAngles', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('PlotFlag', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('RandomSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

% 新增尾流/NSGA参数
p.addParameter('Mode', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('WindSpeedMps', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('WindDirFromDeg', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.addParameter('RotorDiameterM', 242, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('AvailabilityFactor', 1.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
p.addParameter('AmbientTI', 0.0555, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
p.addParameter('DensityReferenceMaxSpacingD', 11, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('UseBinnedWindStates', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('WindSpeedBinWidthMps', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('WindDirectionBinWidthDeg', 15, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 180);

% 加速相关参数
p.addParameter('UseParallel', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('AutoStartParallelPool', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('NumWorkers', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
p.addParameter('ParallelMinUniqueLayouts', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
p.addParameter('PrecomputePairGeometry', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('PrecomputeGeometryMaxMB', 512, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);

p.addParameter('NSGARotorSamples', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('FinalRotorSamples', 9, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('FinalUseRepresentativeStates', true, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('FinalWindSpeedBinWidthMps', 0.5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('FinalWindDirectionBinWidthDeg', 10, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 180);
p.addParameter('FinalStoreTurbineStates', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('ParetoSelection', 'maxAEP', @(x) ischar(x) || isstring(x));
p.addParameter('ReferencePackingTrials', 60, @(x) isnumeric(x) && isscalar(x) && x >= 6);
p.addParameter('UpperCountMargin', 0.15, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.parse(feasible_region_info, turbine_spacing, varargin{:});
prm = p.Results;

if ~isempty(prm.RandomSeed)
    rng(prm.RandomSeed);
end

% =========================== 判断运行模式 ===============================
mode = lower(string(prm.Mode));
hasWind = ~isempty(prm.WindSpeedMps) || ~isempty(prm.WindDirFromDeg);
if mode == "auto"
    if hasWind
        mode = "wake-nsga2";
    else
        mode = "legacy";
    end
end
if ~(mode == "legacy" || mode == "wake-nsga2")
    error('Mode must be ''auto'', ''legacy'', or ''wake-nsga2''.');
end

if mode == "wake-nsga2"
    if isempty(prm.WindSpeedMps) || isempty(prm.WindDirFromDeg)
        error('wake-nsga2 mode requires both WindSpeedMps and WindDirFromDeg.');
    end
    [best_capacity_value, wind_tubine_point, result] = ...
        runWakeAwareNSGA(feasible_region_info, turbine_spacing, prm);
else
    [best_capacity_value, wind_tubine_point, result] = ...
        runLegacyCapacityGA(feasible_region_info, turbine_spacing, prm);
end
end

%% ====================== Wake-aware Random-key NSGA-II ===================
function [best_capacity_value, wind_tubine_point, result] = runWakeAwareNSGA(feasible_region_info, turbine_spacing, prm)

% -------------------------- 1. 可行域 -----------------------------------
feasiblePS = buildFeasiblePolyshape(feasible_region_info);
if area(feasiblePS) <= 0
    error('输入的真实可行域面积为0，无法进行微观选址。');
end

D = prm.RotorDiameterM;
Pr = prm.TurbinePowerMW;
dMin = turbine_spacing;
dRefMax = prm.DensityReferenceMaxSpacingD * D;
if dRefMax < dMin
    error('DensityReferenceMaxSpacingD*RotorDiameterM must be >= turbine_spacing.');
end

% -------------------------- 2. 7D候选点库 -------------------------------
candidate_points = generateHexCandidateLibrary(feasiblePS, dMin, prm.NumAngles);
M = size(candidate_points, 1);
if M == 0
    % Isolated repair: only an empty original lattice uses the solver-backed
    % candidate fallback. All nonempty-lattice calculations remain unchanged.
    [candidate_points, candidateFallbackInfo] = ...
        solver_feasible_candidate_fallback(feasiblePS, dMin);
    M = size(candidate_points, 1);
end
conflict = buildConflictMatrix(candidate_points, dMin);

% 预计算候选点两两坐标差；在所有代、所有个体之间复用。
[candidatePairGeom, geometryInfo] = buildCandidatePairGeometry(candidate_points, prm);

% -------------------------- 3. 参考密度边界 -----------------------------
% N_reference_7D / N_reference_11D 都是启发式几何容量参考值，不声称全局最优。
Nref7 = estimateMaxPackingCount(candidate_points, conflict, prm.ReferencePackingTrials);

candRef11 = generateHexCandidateLibrary(feasiblePS, dRefMax, prm.NumAngles);
if isempty(candRef11)
    Nref11 = 1;
else
    conflictRef11 = buildConflictMatrix(candRef11, dRefMax);
    Nref11 = estimateMaxPackingCount(candRef11, conflictRef11, prm.ReferencePackingTrials);
end

Nlower = max(1, min(Nref11, Nref7));
% 给7D启发式容量留出少量上探空间，避免把参考解误当成严格上界。
Nupper = min(M, max(Nlower, ceil(Nref7*(1 + prm.UpperCountMargin))));
Nupper = max(Nupper, min(M, Nref7 + 2));

% -------------------------- 4. 风况 -------------------------------------
Ufull = double(prm.WindSpeedMps(:));
WDfull = mod(double(prm.WindDirFromDeg(:)), 360);
if numel(Ufull) ~= numel(WDfull)
    error('WindSpeedMps and WindDirFromDeg must have the same length.');
end
if any(~isfinite(Ufull)) || any(Ufull < 0) || any(~isfinite(WDfull))
    error('WindSpeedMps/WindDirFromDeg must be finite; wind speed must be >=0.');
end

if logical(prm.UseBinnedWindStates)
    windStates = buildBinnedWindStates(Ufull, WDfull, ...
        prm.WindSpeedBinWidthMps, prm.WindDirectionBinWidthDeg);
    Ueval = windStates.wind_speed_mps;
    WDeval = windStates.wind_dir_from_deg;
    stateHours = windStates.hours;
else
    Ueval = Ufull;
    WDeval = WDfull;
    stateHours = ones(size(Ufull));
    windStates = struct();
    windStates.wind_speed_mps = Ueval;
    windStates.wind_dir_from_deg = WDeval;
    windStates.hours = stateHours;
    windStates.num_states = numel(Ueval);
    windStates.description = 'full input states used inside NSGA';
end

% -------------------------- 5. 风机模型 ---------------------------------
CtData = GWH242_12MW_Ct_proxy_IEA15();
turb = struct();
turb.rotor_diameter_m = D;
turb.rated_power_MW = Pr;
turb.CtFcn = CtData.evaluate;
turb.PowerFcn = @(U) GWH242_12MW_power_curve(U, 'RatedPowerMW', Pr, ...
    'AvailabilityFactor', prm.AvailabilityFactor);

% 并行环境只初始化一次；无并行工具箱时自动回退串行。
parallelInfo = prepareParallelExecution(prm);

% -------------------------- 6. 初始化NSGA-II ----------------------------
popSize = prm.PopulationSize;
maxGen = prm.MaxGenerations;
nGene = M + 1; % 最后一位是目标风机数量/开发密度随机键

population = rand(popSize, nGene);
population(:,1:M) = injectHeuristicSeeds(population(:,1:M), candidate_points);
% 显式铺开密度键，使初代覆盖N_11D至高密度端。
population(:,end) = linspace(0,1,popSize).';
population(:,end) = population(randperm(popSize),end);

% 缓存：不同random-key可能解码成同一布局，且同一布局可跨代重复出现。
cache = containers.Map('KeyType','char','ValueType','any');
evalStatsTotal = emptyEvaluationStats();
[F, met, cache, batchStats] = evaluateWakePopulation(population, candidate_points, conflict, ...
    Nlower, Nupper, turb, Ueval, WDeval, stateHours, prm, cache, ...
    candidatePairGeom, parallelInfo);
evalStatsTotal = addEvaluationStats(evalStatsTotal, batchStats);

max_capacity_history = zeros(maxGen,1);
min_wake_history = nan(maxGen,1);
max_AEP_history = zeros(maxGen,1);
pareto_size_history = zeros(maxGen,1);
mean_capacity_history = zeros(maxGen,1);

stallCounter = 0;
bestFrontAEP = -inf;

% -------------------------- 7. NSGA-II主循环 ----------------------------
for gen = 1:maxGen
    [rank, crowd] = nondominatedRankAndCrowding(F);
    front = find(rank == 1);

    caps = met.count * Pr;
    max_capacity_history(gen) = max(caps);
    mean_capacity_history(gen) = mean(caps);
    if ~isempty(front)
        min_wake_history(gen) = min(met.wake_loss(front));
        max_AEP_history(gen) = max(met.AEP_wake_MWh(front));
        pareto_size_history(gen) = numel(front);
        currentBestFrontAEP = max_AEP_history(gen);
    else
        currentBestFrontAEP = -inf;
    end

    if ~isfinite(bestFrontAEP) || currentBestFrontAEP > bestFrontAEP + max(1e-6, abs(bestFrontAEP)*1e-8)
        bestFrontAEP = currentBestFrontAEP;
        stallCounter = 0;
    else
        stallCounter = stallCounter + 1;
    end

    % 自适应变异率：前高后低，停滞时提升。
    mu = prm.MutationStart + (prm.MutationEnd-prm.MutationStart) * ...
        (gen-1)/max(1,maxGen-1);
    if stallCounter >= prm.StallTrigger
        mu = max(mu, prm.MutationBoost);
    end

    % 生成与父代同规模的后代。
    offspring = zeros(size(population));
    w = 1;
    while w <= popSize
        i1 = nsgaTournamentSelect(rank, crowd);
        i2 = nsgaTournamentSelect(rank, crowd);
        p1 = population(i1,:);
        p2 = population(i2,:);

        if rand < prm.CrossoverRate
            [c1,c2] = uniformCrossover(p1,p2);
        else
            c1 = p1; c2 = p2;
        end
        c1 = randomResetMutation(c1, mu);
        c2 = randomResetMutation(c2, mu);

        offspring(w,:) = c1;
        if w+1 <= popSize
            offspring(w+1,:) = c2;
        end
        w = w + 2;
    end

    [Foff, metOff, cache, batchStats] = evaluateWakePopulation(offspring, candidate_points, conflict, ...
        Nlower, Nupper, turb, Ueval, WDeval, stateHours, prm, cache, ...
        candidatePairGeom, parallelInfo);
    evalStatsTotal = addEvaluationStats(evalStatsTotal, batchStats);

    % 环境选择：父代+子代 -> 按非支配等级和拥挤距离保留popSize个。
    combPop = [population; offspring];
    combF = [F; Foff];
    combMet = concatMetrics(met, metOff);
    [combRank, combCrowd] = nondominatedRankAndCrowding(combF);
    keep = environmentalSelect(combRank, combCrowd, popSize);

    population = combPop(keep,:);
    F = combF(keep,:);
    met = subsetMetrics(combMet, keep);
end

% -------------------------- 8. 提取最终Pareto前沿 ------------------------
[rank, crowd] = nondominatedRankAndCrowding(F); %#ok<ASGLU>
front = find(rank == 1);

paretoKeys = strings(0,1);
paretoChrom = zeros(0,nGene);
paretoCount = zeros(0,1);
paretoWake = zeros(0,1);
paretoAEP = zeros(0,1);
paretoLayouts = {};
paretoMasks = {};

for kk = 1:numel(front)
    ii = front(kk);
    [mask, nT] = decodeChromosome(population(ii,:), conflict, Nlower, Nupper);
    if nT < Nlower
        continue;
    end
    key = layoutMaskKey(mask);
    if any(paretoKeys == string(key))
        continue;
    end
    paretoKeys(end+1,1) = string(key); %#ok<AGROW>
    paretoChrom(end+1,:) = population(ii,:); %#ok<AGROW>
    paretoCount(end+1,1) = nT; %#ok<AGROW>
    paretoWake(end+1,1) = met.wake_loss(ii); %#ok<AGROW>
    paretoAEP(end+1,1) = met.AEP_wake_MWh(ii); %#ok<AGROW>
    paretoMasks{end+1,1} = mask; %#ok<AGROW>
    paretoLayouts{end+1,1} = candidate_points(mask,:); %#ok<AGROW>
end

if isempty(paretoCount)
    error('NSGA-II finished without a valid Pareto layout. Increase population/generations or candidate density.');
end

% 按容量从小到大、同容量下尾流从小到大排序，方便画Pareto曲线。
[~,ord] = sortrows([paretoCount*Pr, paretoWake],[1 2]);
paretoChrom = paretoChrom(ord,:);
paretoCount = paretoCount(ord);
paretoWake = paretoWake(ord);
paretoAEP = paretoAEP(ord);
paretoLayouts = paretoLayouts(ord);
paretoMasks = paretoMasks(ord);

% -------------------------- 9. 兼容旧主函数：选一个Pareto点 -------------
sel = selectParetoSolution(paretoCount*Pr, paretoWake, paretoAEP, prm.ParetoSelection);
wind_tubine_point = paretoLayouts{sel};
best_capacity_value = paretoCount(sel) * Pr;

% ------------------------- 10. 最终8760 h功率曲线 -----------------------
% 默认不再逐小时执行8760次尾流传播。先用更细的代表风况计算尾流效率
% eta_wake(U,theta)，再按每小时所属风况映射回完整时序。逐小时自由来流
% 功率仍由真实U_t计算，因此只近似尾流效率，不把原功率曲线阶梯化。
selectedMask = paretoMasks{sel};
selectedCandidateIdx = find(selectedMask);
selectedGeom = subsetCandidatePairGeometry(candidatePairGeom, selectedCandidateIdx);

if logical(prm.FinalUseRepresentativeStates)
    finalWindStates = buildBinnedWindStates(Ufull, WDfull, ...
        prm.FinalWindSpeedBinWidthMps, prm.FinalWindDirectionBinWidthDeg);

    selectedWakeRepresentative = bastankhah_gaussian_wake( ...
        wind_tubine_point, finalWindStates.wind_speed_mps, ...
        finalWindStates.wind_dir_from_deg, turb, ...
        'AmbientTI', prm.AmbientTI, ...
        'StateWeightsHours', finalWindStates.hours, ...
        'RotorSamples', prm.FinalRotorSamples, ...
        'StoreTurbineStates', false, ...
        'PrecomputedGeometry', selectedGeom);

    selectedWakeFull = reconstructHourlyPowerFromRepresentativeStates( ...
        selectedWakeRepresentative, finalWindStates, Ufull, turb, ...
        best_capacity_value, prm.AmbientTI);
else
    selectedWakeRepresentative = struct();
    finalWindStates = struct();
    selectedWakeFull = bastankhah_gaussian_wake( ...
        wind_tubine_point, Ufull, WDfull, turb, ...
        'AmbientTI', prm.AmbientTI, ...
        'RotorSamples', prm.FinalRotorSamples, ...
        'StoreTurbineStates', logical(prm.FinalStoreTurbineStates), ...
        'PrecomputedGeometry', selectedGeom);
end

if best_capacity_value > 0
    puWakeFull = selectedWakeFull.farm_power_wake_MW ./ best_capacity_value;
    puFreeFull = selectedWakeFull.farm_power_free_MW ./ best_capacity_value;
else
    puWakeFull = zeros(size(Ufull));
    puFreeFull = zeros(size(Ufull));
end

% ----------------------------- 11. 输出 ---------------------------------
result = struct();
result.mode = 'wake-nsga2';
result.algorithm = 'random-key NSGA-II with Bastankhah Gaussian wake';
result.best_turbine_count = paretoCount(sel); % 兼容旧字段名：实际是selected
result.best_capacity_value = best_capacity_value;
result.wind_tubine_point = wind_tubine_point;
result.candidate_points = candidate_points;
if exist('candidateFallbackInfo','var')
    result.candidate_fallback = candidateFallbackInfo;
end
result.feasible_polyshape = feasiblePS;
result.conflict_matrix = conflict;
result.parameters = prm;
result.rotor_diameter_m = D;
result.availability_factor = prm.AvailabilityFactor;
result.minimum_spacing_m = dMin;
result.minimum_spacing_D = dMin / D;
result.density_reference_max_spacing_m = dRefMax;
result.density_reference_max_spacing_D = dRefMax / D;
result.N_reference_7D = Nref7;
result.N_reference_11D = Nref11;
result.N_target_lower = Nlower;
result.N_target_upper = Nupper;

result.wind_states_binned = windStates;
result.wind_full_speed_mps = Ufull;
result.wind_full_dir_from_deg = WDfull;
result.geometry_precompute = geometryInfo;
result.parallel_info = parallelInfo;

result.pareto_chromosomes = paretoChrom;
result.pareto_turbine_count = paretoCount;
result.pareto_capacity_MW = paretoCount * Pr;
result.pareto_wake_loss_fraction = paretoWake;
result.pareto_AEP_wake_MWh_binned = paretoAEP;
result.pareto_layouts = paretoLayouts;
result.pareto_masks = paretoMasks;
result.pareto_selection_method = char(string(prm.ParetoSelection));
result.selected_pareto_index = sel;
result.selected_wake_binned_loss_fraction = paretoWake(sel);
result.selected_AEP_wake_MWh_binned = paretoAEP(sel);
result.selected_wake_full = selectedWakeFull;
result.selected_wake_representative = selectedWakeRepresentative;
result.final_wind_states = finalWindStates;
result.wind_power_pu_wake_full = puWakeFull;
result.wind_power_pu_free_full = puFreeFull;

% 兼容原结果绘图/调试字段：这里含义改为容量历史。
result.best_history = max_capacity_history;
result.mean_history = mean_capacity_history;
result.max_capacity_history_MW = max_capacity_history;
result.mean_capacity_history_MW = mean_capacity_history;
result.min_wake_history = min_wake_history;
result.max_AEP_history_MWh = max_AEP_history;
result.pareto_size_history = pareto_size_history;
result.cache_size = cache.Count;
result.evaluation_stats = evalStatsTotal;
if evalStatsTotal.requested_valid_layouts > 0
    result.cache_hit_ratio = evalStatsTotal.cache_hits / evalStatsTotal.requested_valid_layouts;
    result.unique_wake_evaluation_ratio = evalStatsTotal.unique_wake_evaluations / evalStatsTotal.requested_valid_layouts;
else
    result.cache_hit_ratio = 0;
    result.unique_wake_evaluation_ratio = 0;
end

% ----------------------------- 12. 可视化 -------------------------------
if logical(prm.PlotFlag)
    figure('Name','Wake-aware NSGA-II 微观选址'); hold on; axis equal;
    plot(feasiblePS, 'FaceColor',[0.88 0.92 1.00], 'FaceAlpha',0.35, ...
        'EdgeColor','b','LineWidth',1.1);
    plot(candidate_points(:,1),candidate_points(:,2),'.', ...
        'Color',[0.72 0.72 0.72],'MarkerSize',6);
    plot(wind_tubine_point(:,1),wind_tubine_point(:,2),'ro', ...
        'MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.0);
    title(sprintf('Selected Pareto layout: N=%d, P=%.0f MW, wake=%.2f%%', ...
        paretoCount(sel),best_capacity_value,100*selectedWakeFull.wake_loss_fraction));
    xlabel('x (m)'); ylabel('y (m)'); grid on; box on;

    figure('Name','Capacity-wake Pareto front');
    scatter(result.pareto_capacity_MW,100*result.pareto_wake_loss_fraction,36,'filled'); hold on;
    scatter(best_capacity_value,100*paretoWake(sel),80,'o','LineWidth',1.5);
    xlabel('Installed capacity (MW)'); ylabel('Wake loss (%)'); grid on; box on;
    title('NSGA-II Pareto front');

    figure('Name','NSGA-II convergence'); hold on; box on; grid on;
    yyaxis left;
    plot(max_capacity_history,'LineWidth',1.25);
    ylabel('Max capacity in population (MW)');
    yyaxis right;
    plot(100*min_wake_history,'LineWidth',1.25);
    ylabel('Minimum wake loss on front (%)');
    xlabel('Generation');
    title('NSGA-II convergence indicators');
end
end

%% =========================== Legacy GA =================================
function [best_capacity_value, wind_tubine_point, result] = runLegacyCapacityGA(feasible_region_info, turbine_spacing, prm)
% 原函数逻辑基本保持不变，用于未提供风况数据时的兼容运行。

feasiblePS = buildFeasiblePolyshape(feasible_region_info);
if area(feasiblePS) <= 0
    error('输入的真实可行域面积为0，无法进行微观选址。');
end

candidate_points = generateHexCandidateLibrary(feasiblePS, turbine_spacing, prm.NumAngles);
M = size(candidate_points,1);
if M == 0
    best_capacity_value = 0;
    wind_tubine_point = zeros(0,2);
    result = struct('message','可行域内未生成任何候选点','mode','legacy');
    return;
end
conflict = buildConflictMatrix(candidate_points,turbine_spacing);

popSize = prm.PopulationSize;
maxGen = prm.MaxGenerations;
eliteCnt = max(1,round(prm.EliteRatio*popSize));
population = rand(popSize,M);
population = injectHeuristicSeeds(population,candidate_points);
fitness = zeros(popSize,1);
decoded = cell(popSize,1);

bestHistory = zeros(maxGen,1);
meanHistory = zeros(maxGen,1);
bestCount = -inf;
bestMask = false(1,M);
bestKeys = [];
stallCounter = 0;

for gen = 1:maxGen
    mu = prm.MutationStart + (prm.MutationEnd-prm.MutationStart) * ...
        (gen-1)/max(1,maxGen-1);
    if stallCounter >= prm.StallTrigger
        mu = max(mu,prm.MutationBoost);
    end

    for ii = 1:popSize
        mask = decodeRandomKeys(population(ii,:),conflict);
        decoded{ii} = mask;
        fitness(ii) = sum(mask);
    end

    [genBest,idxBest] = max(fitness);
    if genBest > bestCount
        bestCount = genBest;
        bestMask = decoded{idxBest};
        bestKeys = population(idxBest,:);
        stallCounter = 0;
    else
        stallCounter = stallCounter + 1;
    end
    bestHistory(gen) = bestCount;
    meanHistory(gen) = mean(fitness);

    [~,ord] = sort(fitness,'descend');
    newPop = zeros(size(population));
    newPop(1:eliteCnt,:) = population(ord(1:eliteCnt),:);
    w = eliteCnt+1;
    while w <= popSize
        p1 = legacyTournamentSelect(fitness,3);
        p2 = legacyTournamentSelect(fitness,3);
        if rand < prm.CrossoverRate
            [c1,c2] = uniformCrossover(population(p1,:),population(p2,:));
        else
            c1 = population(p1,:); c2 = population(p2,:);
        end
        c1 = randomResetMutation(c1,mu);
        c2 = randomResetMutation(c2,mu);
        newPop(w,:) = c1;
        if w+1 <= popSize
            newPop(w+1,:) = c2;
        end
        w = w+2;
    end
    population = newPop;
end

wind_tubine_point = candidate_points(bestMask,:);
best_capacity_value = bestCount*prm.TurbinePowerMW;

result = struct();
result.mode = 'legacy';
result.best_turbine_count = bestCount;
result.best_capacity_value = best_capacity_value;
result.wind_tubine_point = wind_tubine_point;
result.candidate_points = candidate_points;
result.feasible_polyshape = feasiblePS;
result.conflict_matrix = conflict;
result.best_history = bestHistory;
result.mean_history = meanHistory;
result.best_keys = bestKeys;
result.parameters = prm;

if logical(prm.PlotFlag)
    figure('Name','Legacy GA 微观选址结果'); hold on; axis equal;
    plot(feasiblePS,'FaceColor',[0.88 0.92 1.00],'FaceAlpha',0.4,'EdgeColor','b','LineWidth',1.2);
    plot(candidate_points(:,1),candidate_points(:,2),'.','Color',[0.6 0.6 0.6],'MarkerSize',8);
    plot(wind_tubine_point(:,1),wind_tubine_point(:,2),'ro','MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.0);
    title(sprintf('最优微观选址：%d台，%.2f MW',bestCount,best_capacity_value));
    xlabel('x'); ylabel('y'); grid on; box on;
end
end

%% ======================== NSGA population evaluation ===================
function [F, met, cache, stats] = evaluateWakePopulation(population, candidate_points, conflict, ...
    Nlower, Nupper, turb, Ueval, WDeval, stateHours, prm, cache, ...
    candidatePairGeom, parallelInfo)
% 三层加速：批量解码 -> cache/批内去重 -> 仅评价未缓存唯一布局。

nPop = size(population,1);
F = zeros(nPop,2);
met.count = zeros(nPop,1);
met.wake_loss = nan(nPop,1);
met.AEP_wake_MWh = zeros(nPop,1);
met.valid = false(nPop,1);

masks = cell(nPop,1);
layoutKeys = cell(nPop,1);
validMask = false(nPop,1);

% 1) 解码全部个体。
for ii = 1:nPop
    [mask,nT] = decodeChromosome(population(ii,:),conflict,Nlower,Nupper);
    masks{ii} = mask;
    met.count(ii) = nT;
    if nT < Nlower
        shortage = Nlower-nT;
        F(ii,1) = 1e6 + shortage*1e3;
        F(ii,2) = 1 + shortage/max(1,Nlower);
        met.wake_loss(ii) = F(ii,2);
    else
        validMask(ii) = true;
        layoutKeys{ii} = layoutMaskKey(mask);
    end
end

validIdx = find(validMask);
stats = emptyEvaluationStats();
stats.requested_layouts = nPop;
stats.requested_valid_layouts = numel(validIdx);
if isempty(validIdx)
    return;
end

% 2) 查询跨代cache；尚未缓存的布局再做批内去重。
missingIdx = zeros(0,1);
for kk = 1:numel(validIdx)
    ii = validIdx(kk);
    if isKey(cache,layoutKeys{ii})
        stats.cache_hits = stats.cache_hits + 1;
    else
        missingIdx(end+1,1) = ii; %#ok<AGROW>
    end
end

if ~isempty(missingIdx)
    missingKeys = layoutKeys(missingIdx);
    [uniqueMissingKeys,ia] = unique(missingKeys,'stable');
    repPopIdx = missingIdx(ia);
    nUnique = numel(repPopIdx);
    stats.duplicate_layouts_in_batch = numel(missingIdx)-nUnique;
    stats.unique_wake_evaluations = nUnique;

    wakeLossNew = zeros(nUnique,1);
    AEPwakeNew = zeros(nUnique,1);
    nTNew = zeros(nUnique,1);

    doParallel = parallelInfo.enabled && nUnique >= prm.ParallelMinUniqueLayouts;
    stats.parallel_batches = double(doParallel);
    stats.serial_batches = double(~doParallel);

    if doParallel
        parfor uu = 1:nUnique
            ii = repPopIdx(uu);
            [wakeLossNew(uu),AEPwakeNew(uu),nTNew(uu)] = evaluateOneLayoutWake( ...
                masks{ii},candidate_points,candidatePairGeom,turb,Ueval,WDeval,stateHours,prm);
        end
    else
        for uu = 1:nUnique
            ii = repPopIdx(uu);
            [wakeLossNew(uu),AEPwakeNew(uu),nTNew(uu)] = evaluateOneLayoutWake( ...
                masks{ii},candidate_points,candidatePairGeom,turb,Ueval,WDeval,stateHours,prm);
        end
    end

    % Map只在主进程中写入，避免parfor共享状态。
    for uu = 1:nUnique
        cache(uniqueMissingKeys{uu}) = [wakeLossNew(uu),AEPwakeNew(uu),nTNew(uu)];
    end
end

% 3) 从cache统一回填种群目标值。
for kk = 1:numel(validIdx)
    ii = validIdx(kk);
    val = cache(layoutKeys{ii});
    wakeLoss = val(1);
    AEPwake = val(2);
    nT = met.count(ii);
    F(ii,1) = -nT*prm.TurbinePowerMW;
    F(ii,2) = wakeLoss;
    met.wake_loss(ii) = wakeLoss;
    met.AEP_wake_MWh(ii) = AEPwake;
    met.valid(ii) = true;
end
end

function [wakeLoss,AEPwake,nT] = evaluateOneLayoutWake(mask,candidate_points,candidatePairGeom, ...
    turb,Ueval,WDeval,stateHours,prm)
idx = find(mask);
nT = numel(idx);
xy = candidate_points(idx,:);
layoutGeom = subsetCandidatePairGeometry(candidatePairGeom,idx);

R = bastankhah_gaussian_wake(xy,Ueval,WDeval,turb, ...
    'AmbientTI',prm.AmbientTI, ...
    'StateWeightsHours',stateHours, ...
    'RotorSamples',prm.NSGARotorSamples, ...
    'StoreTurbineStates',false, ...
    'PrecomputedGeometry',layoutGeom);

wakeLoss = R.wake_loss_fraction;
AEPwake = R.AEP_wake_MWh;
if ~isfinite(wakeLoss); wakeLoss = 1; end
if ~isfinite(AEPwake); AEPwake = 0; end
end

%% ======================== Pareto / NSGA utilities ======================
function [rank,crowd] = nondominatedRankAndCrowding(F)
% 标准NSGA-II非支配排序与拥挤距离；所有目标均按最小化解释。
n = size(F,1);
rank = inf(n,1);
S = cell(n,1);
domCount = zeros(n,1);
fronts = cell(0,1);
first = [];

for p = 1:n
    Sp = [];
    np = 0;
    for q = 1:n
        if p == q
            continue;
        end
        if dominates(F(p,:),F(q,:))
            Sp(end+1) = q; %#ok<AGROW>
        elseif dominates(F(q,:),F(p,:))
            np = np + 1;
        end
    end
    S{p} = Sp;
    domCount(p) = np;
    if np == 0
        rank(p) = 1;
        first(end+1) = p; %#ok<AGROW>
    end
end
fronts{1} = first;

f = 1;
while f <= numel(fronts) && ~isempty(fronts{f})
    next = [];
    for p = fronts{f}
        for q = S{p}
            domCount(q) = domCount(q)-1;
            if domCount(q) == 0
                rank(q) = f+1;
                next(end+1) = q; %#ok<AGROW>
            end
        end
    end
    if ~isempty(next)
        fronts{f+1} = unique(next,'stable'); %#ok<AGROW>
    end
    f = f+1;
end

crowd = zeros(n,1);
for ff = 1:numel(fronts)
    idx = fronts{ff};
    if isempty(idx)
        continue;
    end
    crowd(idx) = crowdingDistance(F(idx,:));
end
end

function tf = dominates(a,b)
tf = all(a <= b) && any(a < b);
end

function d = crowdingDistance(Ffront)
n = size(Ffront,1);
m = size(Ffront,2);
d = zeros(n,1);
if n <= 2
    d(:) = inf;
    return;
end
for jj = 1:m
    [vals,ord] = sort(Ffront(:,jj),'ascend');
    d(ord(1)) = inf;
    d(ord(end)) = inf;
    span = vals(end)-vals(1);
    if span <= 0 || ~isfinite(span)
        continue;
    end
    for kk = 2:n-1
        if ~isinf(d(ord(kk)))
            d(ord(kk)) = d(ord(kk)) + (vals(kk+1)-vals(kk-1))/span;
        end
    end
end
end

function idx = nsgaTournamentSelect(rank,crowd)
n = numel(rank);
a = randi(n); b = randi(n);
if rank(a) < rank(b)
    idx = a;
elseif rank(b) < rank(a)
    idx = b;
elseif crowd(a) > crowd(b)
    idx = a;
elseif crowd(b) > crowd(a)
    idx = b;
else
    if rand < 0.5; idx = a; else; idx = b; end
end
end

function keep = environmentalSelect(rank,crowd,nKeep)
idx = (1:numel(rank)).';
T = [idx,rank,-crowd];
T = sortrows(T,[2 3]);
keep = T(1:nKeep,1);
end

function out = concatMetrics(a,b)
out.count = [a.count;b.count];
out.wake_loss = [a.wake_loss;b.wake_loss];
out.AEP_wake_MWh = [a.AEP_wake_MWh;b.AEP_wake_MWh];
out.valid = [a.valid;b.valid];
end

function out = subsetMetrics(a,idx)
out.count = a.count(idx);
out.wake_loss = a.wake_loss(idx);
out.AEP_wake_MWh = a.AEP_wake_MWh(idx);
out.valid = a.valid(idx);
end

function sel = selectParetoSolution(capacityMW,wakeLoss,AEPwake,method)
method = lower(string(method));
switch method
    case "maxaep"
        [~,sel] = max(AEPwake);
    case "maxcapacity"
        maxC = max(capacityMW);
        cand = find(capacityMW == maxC);
        [~,k] = min(wakeLoss(cand));
        sel = cand(k);
    case "minwake"
        minW = min(wakeLoss);
        cand = find(abs(wakeLoss-minW) <= 1e-12);
        [~,k] = max(capacityMW(cand));
        sel = cand(k);
    case "knee"
        if max(capacityMW) > min(capacityMW)
            capN = (capacityMW-min(capacityMW))/(max(capacityMW)-min(capacityMW));
        else
            capN = ones(size(capacityMW));
        end
        if max(wakeLoss) > min(wakeLoss)
            wakeN = (wakeLoss-min(wakeLoss))/(max(wakeLoss)-min(wakeLoss));
        else
            wakeN = zeros(size(wakeLoss));
        end
        distIdeal = sqrt((1-capN).^2 + wakeN.^2);
        [~,sel] = min(distIdeal);
    otherwise
        error('Unknown ParetoSelection. Use maxAEP, knee, maxCapacity, or minWake.');
end
end

%% ========================= Chromosome decoding ==========================
function [mask,nSelected,targetN] = decodeChromosome(chrom,conflict,Nlower,Nupper)
M = size(conflict,1);
keys = chrom(1:M);
densityKey = min(max(chrom(M+1),0),1);
targetN = round(Nlower + densityKey*(Nupper-Nlower));
targetN = max(Nlower,min(Nupper,targetN));
mask = decodeRandomKeysTarget(keys,conflict,targetN);
nSelected = sum(mask);
end

function mask = decodeRandomKeysTarget(keys,conflict,targetN)
M = numel(keys);
[~,order] = sort(keys,'descend');
blocked = false(1,M);
mask = false(1,M);
count = 0;
for kk = 1:M
    idx = order(kk);
    if ~blocked(idx)
        mask(idx) = true;
        count = count+1;
        blocked = blocked | conflict(idx,:);
        if count >= targetN
            break;
        end
    end
end
end

function key = layoutMaskKey(mask)
idx = find(mask);
if isempty(idx)
    key = 'empty';
else
    key = sprintf('%d_',idx);
end
end

%% ============================= Wind bins ================================
function S = buildBinnedWindStates(U,WD,dU,dDir)
% 将逐小时风况压缩为(U,WD)二维代表状态。
%
% 代表风向采用固定扇区中心，而不是每个U-WD单元各自的圆平均风向。
% 因而同一方向扇区的所有风速状态具有完全相同的方向值，便于wake函数
% 只计算一次该方向的几何投影并复用于所有风速状态。
%
% S.hour_to_state_index 可将代表状态结果映射回原始逐小时时序。

U = double(U(:));
WD = mod(double(WD(:)),360);
if numel(U) ~= numel(WD)
    error('U and WD must have the same length.');
end

speedBin = floor(U./dU) + 1;
nDir = ceil(360/dDir);
dirBin = floor(WD./dDir) + 1;
dirBin(dirBin > nDir) = nDir;

pairID = (speedBin-1)*nDir + dirBin;
[uID,~,grp] = unique(pairID,'sorted');
hours = accumarray(grp,1);
Umean = accumarray(grp,U,[],@mean);

% 固定方向中心，例如15°分箱 -> 7.5,22.5,...,352.5°。
dirIdx = mod(uID-1,nDir) + 1;
WDcenter = mod((double(dirIdx)-0.5).*dDir,360);

% 实际bin内圆平均风向只用于诊断。
sinMean = accumarray(grp,sind(WD),[],@mean);
cosMean = accumarray(grp,cosd(WD),[],@mean);
WDcircularMean = mod(atan2d(sinMean,cosMean),360);

S = struct();
S.wind_speed_mps = Umean;
S.wind_dir_from_deg = WDcenter;
S.wind_dir_circular_mean_deg = WDcircularMean;
S.hours = hours;
S.num_states = numel(uID);
S.num_unique_directions = numel(unique(WDcenter));
S.wind_speed_bin_width_mps = dU;
S.wind_direction_bin_width_deg = dDir;
S.hour_to_state_index = grp;
S.speed_bin_index = floor((uID-1)./nDir) + 1;
S.direction_bin_index = dirIdx;
S.description = ['U-WD representative states; fixed direction-bin centers permit ', ...
    'direction-grouped wake geometry reuse.'];
end

%% ========================= Acceleration helpers =========================
function [G,info] = buildCandidatePairGeometry(candidate_points,prm)
% 预计算候选点两两坐标差；空间复杂度约为2*M^2个double。
M = size(candidate_points,1);
estMB = 2*double(M)*double(M)*8/1024^2;
info = struct();
info.requested = logical(prm.PrecomputePairGeometry);
info.num_candidates = M;
info.estimated_memory_MB = estMB;
info.memory_limit_MB = prm.PrecomputeGeometryMaxMB;
info.enabled = false;
info.reason = '';

if ~logical(prm.PrecomputePairGeometry)
    G = [];
    info.reason = 'disabled by user';
    return;
end
if estMB > prm.PrecomputeGeometryMaxMB
    G = [];
    info.reason = sprintf('estimated %.1f MB exceeds limit %.1f MB',estMB,prm.PrecomputeGeometryMaxMB);
    return;
end

x = double(candidate_points(:,1));
y = double(candidate_points(:,2));
G = struct();
G.delta_x_m = x - x.';
G.delta_y_m = y - y.';
info.enabled = true;
info.reason = 'candidate pairwise delta matrices precomputed';
end

function Gsub = subsetCandidatePairGeometry(G,idx)
if isempty(G)
    Gsub = [];
    return;
end
Gsub = struct();
Gsub.delta_x_m = G.delta_x_m(idx,idx);
Gsub.delta_y_m = G.delta_y_m(idx,idx);
end

function info = prepareParallelExecution(prm)
% 尝试建立/复用并行池；失败时自动回退串行。
info = struct();
info.requested = logical(prm.UseParallel);
info.enabled = false;
info.pool_workers = 0;
info.reason = '';

if ~info.requested
    info.reason = 'disabled by user';
    return;
end

try
    hasPCT = license('test','Distrib_Computing_Toolbox');
catch
    hasPCT = false;
end
if ~hasPCT
    info.reason = 'Parallel Computing Toolbox unavailable; serial fallback';
    return;
end

try
    pool = gcp('nocreate');
    if isempty(pool) && logical(prm.AutoStartParallelPool)
        if isempty(prm.NumWorkers)
            pool = parpool('local');
        else
            pool = parpool('local',round(prm.NumWorkers));
        end
    end
    if isempty(pool)
        info.reason = 'no active parallel pool; serial fallback';
    else
        info.enabled = true;
        info.pool_workers = pool.NumWorkers;
        info.reason = 'parallel pool active';
    end
catch ME
    info.enabled = false;
    info.reason = ['parallel initialization failed; serial fallback: ',ME.message];
end
end

function S = emptyEvaluationStats()
S = struct();
S.requested_layouts = 0;
S.requested_valid_layouts = 0;
S.cache_hits = 0;
S.duplicate_layouts_in_batch = 0;
S.unique_wake_evaluations = 0;
S.parallel_batches = 0;
S.serial_batches = 0;
end

function A = addEvaluationStats(A,B)
fn = fieldnames(A);
for ii = 1:numel(fn)
    A.(fn{ii}) = A.(fn{ii}) + B.(fn{ii});
end
end

function out = reconstructHourlyPowerFromRepresentativeStates(rep,states,Ufull,turb,installedCapacityMW,ambientTI)
% 用代表风况的相对尾流效率重构完整逐小时功率。
% 每小时自由来流功率仍用真实U_t计算，只近似eta_wake(U,theta)。

stateFree = rep.farm_power_free_MW(:);
stateWake = rep.farm_power_wake_MW(:);
eta = ones(size(stateFree));
valid = stateFree > 1e-12;
eta(valid) = stateWake(valid)./stateFree(valid);
eta = min(max(eta,0),1);

map = states.hour_to_state_index(:);
if numel(map) ~= numel(Ufull)
    error('Representative-state hour mapping length does not match full time series.');
end
etaHour = eta(map);

Pone = turb.PowerFcn(Ufull(:));
Pone = double(Pone(:));
Pone(~isfinite(Pone)) = 0;
Pone = min(max(Pone,0),turb.rated_power_MW);
N = rep.num_turbines;
farmFree = N.*Pone;
farmWake = farmFree.*etaHour;

stateLossHour = nan(size(farmFree));
mask = farmFree > 0;
stateLossHour(mask) = 1-etaHour(mask);
AEPfree = sum(farmFree);
AEPwake = sum(farmWake);
if AEPfree > 0
    wakeEff = AEPwake/AEPfree;
    wakeLoss = 1-wakeEff;
else
    wakeEff = NaN;
    wakeLoss = NaN;
end
if installedCapacityMW > 0
    CF = AEPwake/(installedCapacityMW*numel(Ufull));
else
    CF = NaN;
end

out = rep;
out.model = [rep.model,' + representative-state lookup reconstruction'];
out.num_states = numel(Ufull);
out.represented_hours = numel(Ufull);
out.state_weights_hours = ones(numel(Ufull),1);
out.farm_power_free_MW = farmFree;
out.farm_power_wake_MW = farmWake;
out.state_wake_loss_fraction = stateLossHour;
out.AEP_free_MWh = AEPfree;
out.AEP_wake_MWh = AEPwake;
out.wake_efficiency = wakeEff;
out.wake_loss_fraction = wakeLoss;
out.capacity_factor_wake = CF;
out.ambient_TI = repmat(double(ambientTI),numel(Ufull),1);
out.kstar = 0.3837.*out.ambient_TI + 0.003678;
out.reconstruction_method = 'hourly exact free power x representative-state wake efficiency';
out.representative_num_states = states.num_states;
out.representative_num_unique_directions = states.num_unique_directions;
out.representative_wake_efficiency = eta;
out.hour_to_representative_state = map;
out.hourly_wake_efficiency = etaHour;

if isfield(out,'Ueff_mps'); out = rmfield(out,'Ueff_mps'); end
if isfield(out,'Ct_effective'); out = rmfield(out,'Ct_effective'); end
if isfield(out,'turbine_power_MW'); out = rmfield(out,'turbine_power_MW'); end
end

%% ========================= Reference packing ============================
function bestCount = estimateMaxPackingCount(pts,conflict,nTrials)
M = size(pts,1);
if M == 0
    bestCount = 0;
    return;
end
nTrials = max(6,round(nTrials));

seedPop = rand(nTrials,M);
seedPop = injectHeuristicSeeds(seedPop,pts);
bestCount = 0;
for ii = 1:nTrials
    mask = decodeRandomKeys(seedPop(ii,:),conflict);
    bestCount = max(bestCount,sum(mask));
end
end

%% =========================== Geometry helpers ===========================
function feasiblePS = buildFeasiblePolyshape(info)
% 将输入转换为真实可行域polyshape；对多ring采用奇偶包含深度识别孔洞。
if isa(info,'polyshape')
    feasiblePS = simplify(info);
    return;
end
if ischar(info) || isstring(info)
    S = shaperead(info);
    rings = extractRingsFromStruct(S);
    feasiblePS = ringsToFeasiblePolyshape(rings);
    return;
end
if isstruct(info)
    if all(isfield(info,{'X','Y'}))
        rings = extractRingsFromStruct(info);
        feasiblePS = ringsToFeasiblePolyshape(rings);
        return;
    else
        error('struct输入必须包含字段X和Y。');
    end
end
if isnumeric(info) && size(info,2) == 2
    rings = splitXYByNaN(info(:,1),info(:,2));
    feasiblePS = ringsToFeasiblePolyshape(rings);
    return;
end
error('无法识别feasible_region_info的数据类型。');
end

function rings = extractRingsFromStruct(S)
rings = {};
for kk = 1:numel(S)
    ringsK = splitXYByNaN(S(kk).X(:),S(kk).Y(:));
    rings = [rings,ringsK]; %#ok<AGROW>
end
end

function rings = splitXYByNaN(x,y)
nanMask = isnan(x) | isnan(y);
cut = [0;find(nanMask);numel(x)+1];
rings = {};
for ii = 1:numel(cut)-1
    i1 = cut(ii)+1; i2 = cut(ii+1)-1;
    if i2 < i1; continue; end
    xy = [x(i1:i2),y(i1:i2)];
    xy = xy(all(isfinite(xy),2),:);
    if size(xy,1) < 3; continue; end
    if norm(xy(1,:)-xy(end,:)) > 1e-10
        xy = [xy;xy(1,:)]; %#ok<AGROW>
    end
    rings{end+1} = xy; %#ok<AGROW>
end
end

function feasiblePS = ringsToFeasiblePolyshape(rings)
n = numel(rings);
if n == 0
    feasiblePS = polyshape();
    return;
end
areas = zeros(n,1); depth = zeros(n,1);
for ii = 1:n
    xy = rings{ii};
    areas(ii) = abs(polyarea(xy(:,1),xy(:,2)));
end
for ii = 1:n
    pt = rings{ii}(1,:);
    for jj = 1:n
        if ii == jj || areas(jj) <= areas(ii); continue; end
        [in,on] = inpolygon(pt(1),pt(2),rings{jj}(:,1),rings{jj}(:,2));
        if in || on; depth(ii) = depth(ii)+1; end
    end
end
[~,ord] = sort(areas,'descend');
feasiblePS = polyshape();
for kk = 1:numel(ord)
    ii = ord(kk);
    xy = rings{ii};
    pg = polyshape(xy(:,1),xy(:,2),'Simplify',false);
    if mod(depth(ii),2) == 0
        feasiblePS = union(feasiblePS,pg);
    else
        feasiblePS = subtract(feasiblePS,pg);
    end
end
feasiblePS = simplify(feasiblePS);
end

function cand = generateHexCandidateLibrary(regionPS,d,numAngles)
% 多角度+多平移六边形/三角密铺候选库。
[xmin,xmax,ymin,ymax] = boundingBoxOfPolyshape(regionPS);
ctr = [(xmin+xmax)/2,(ymin+ymax)/2];
dx = d; dy = sqrt(3)/2*d;
L = hypot(xmax-xmin,ymax-ymin) + 2*d;
angles = linspace(0,60,numAngles+1); angles(end)=[];
offsets = [0,0;0.5*dx,0;0,0.5*dy;0.5*dx,0.5*dy];

cand = zeros(0,2);
yrows = -L:dy:L;
xbase = -L:dx:L;
for aa = 1:numel(angles)
    th = deg2rad(angles(aa));
    R = [cos(th),-sin(th);sin(th),cos(th)];
    for oo = 1:size(offsets,1)
        ox = offsets(oo,1); oy = offsets(oo,2);
        pts = zeros(0,2);
        for rr = 1:numel(yrows)
            y = yrows(rr)+oy;
            xshift = ox + mod(rr,2)*dx/2;
            xrow = xbase+xshift;
            pts = [pts;[xrow(:),y*ones(numel(xrow),1)]]; %#ok<AGROW>
        end
        pts = (R*pts.').';
        pts = pts+ctr;
        inside = isinterior(regionPS,pts(:,1),pts(:,2));
        cand = [cand;pts(inside,:)]; %#ok<AGROW>
    end
end
if isempty(cand); return; end
tol = d*1e-3;
cand = round(cand/tol)*tol;
cand = unique(cand,'rows');
end

function conflict = buildConflictMatrix(pts,dmin)
n = size(pts,1);
if n == 0
    conflict = false(0,0);
    return;
end
% 避免依赖pdist2：直接使用坐标差构造欧氏距离平方。
dx = pts(:,1)-pts(:,1).';
dy = pts(:,2)-pts(:,2).';
d2 = dx.^2+dy.^2;
conflict = d2 < (dmin-1e-9)^2;
conflict(1:n+1:end) = true;
end

function [xmin,xmax,ymin,ymax] = boundingBoxOfPolyshape(ps)
[bx,by] = boundary(ps);
good = isfinite(bx) & isfinite(by);
bx = bx(good); by = by(good);
xmin = min(bx); xmax = max(bx);
ymin = min(by); ymax = max(by);
end

%% ============================ GA helpers ================================
function pop = injectHeuristicSeeds(pop,pts)
[nPop,nGene] = size(pop);
if isempty(pts) || nPop < 1 || nGene < 2
    return;
end
ctr = mean(pts,1);
dist2ctr = vecnorm(pts-ctr,2,2);
[~,o1] = sortrows(pts,[1 2]);
[~,o2] = sortrows(pts,[-1 2]);
[~,o3] = sortrows(pts,[2 1]);
[~,o4] = sortrows(pts,[-2 1]);
[~,o5] = sort(dist2ctr,'ascend');
[~,o6] = sort(dist2ctr,'descend');
orders = {o1,o2,o3,o4,o5,o6};
for kk = 1:min(numel(orders),nPop)
    pop(kk,:) = orderToKeys(orders{kk},nGene);
end
end

function keys = orderToKeys(order,n)
keys = zeros(1,n);
keys(order) = linspace(1,0,n);
end

function mask = decodeRandomKeys(keys,conflict)
M = numel(keys);
[~,order] = sort(keys,'descend');
blocked = false(1,M);
mask = false(1,M);
for kk = 1:M
    idx = order(kk);
    if ~blocked(idx)
        mask(idx) = true;
        blocked = blocked | conflict(idx,:);
    end
end
end

function idx = legacyTournamentSelect(fitness,tsize)
n = numel(fitness);
cand = randi(n,[tsize,1]);
[~,loc] = max(fitness(cand));
idx = cand(loc);
end

function [c1,c2] = uniformCrossover(p1,p2)
mask = rand(size(p1)) < 0.5;
c1 = p1; c2 = p2;
c1(mask) = p2(mask);
c2(mask) = p1(mask);
end

function child = randomResetMutation(child,mu)
mut = rand(size(child)) < mu;
child(mut) = rand(1,nnz(mut));
end
