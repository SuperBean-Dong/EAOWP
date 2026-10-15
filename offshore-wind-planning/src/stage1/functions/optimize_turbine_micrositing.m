function out = optimize_turbine_micrositing(feasible_region_info, wind_speed_140m_mps, ...
    wind_dir_from_deg, turbine, availability_factor, varargin)
%OPTIMIZE_TURBINE_MICROSITING Per-farm wake-aware RK-NSGA micrositing.
%
% Power is recalculated directly from the supplied 140-m wind speed using
% GWH242_12MW_power_curve.m. The turbine availability factor is an explicit
% function input and is therefore not hidden inside the power-curve model.
%
% Core NSGA-II objectives:
%   1) maximize installed capacity;
%   2) minimize annual wake-loss fraction.
% The final single layout is selected from the Pareto front by maximum AEP.
%
% INPUTS
%   feasible_region_info : actual nonconvex/perforated feasible polygon
%   wind_speed_140m_mps  : T x 1 hub-height wind speed [m/s]
%   wind_dir_from_deg    : T x 1 meteorological FROM direction [deg]
%   turbine              : struct with rated_power_MW, rotor_diameter_m,
%                          minimum_spacing_D
%   availability_factor  : turbine availability in [0,1]
%
% OUTPUT
%   Includes selected layout, full Pareto result, recalculated no-wake 2023
%   power p.u. and wake-adjusted 2023 power p.u. for Stage II.

p = inputParser;
p.addParameter('PopulationSize',80,@(x)isnumeric(x)&&isscalar(x)&&x>=20);
p.addParameter('MaxGenerations',150,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
p.addParameter('AmbientTI',0.0555,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('RandomSeed',1,@(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('UseParallel',true,@(x)islogical(x)||(isnumeric(x)&&isscalar(x)));
p.addParameter('NumWorkers',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=1));
p.addParameter('ParetoSelection','maxAEP',@(x)ischar(x)||isstring(x));
p.parse(varargin{:});
opt = p.Results;

validateattributes(availability_factor,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1},mfilename,'availability_factor');

U = double(wind_speed_140m_mps(:));
WD = mod(double(wind_dir_from_deg(:)),360);
if numel(U) ~= numel(WD)
    error('wind speed and wind direction must have the same length.');
end
if any(~isfinite(U)) || any(U<0) || any(~isfinite(WD))
    error('Wind speed/direction contains invalid values.');
end

Pr = turbine.rated_power_MW;
D = turbine.rotor_diameter_m;
dMin = turbine.minimum_spacing_D * D;

% Exact no-wake 2023 power reconstructed from wind speed using the supplied
% project curve. pRaw is actual available p.u. output, i.e. availability is
% already included.
[PoneRaw_MW,~,pRaw] = GWH242_12MW_power_curve(U, ...
    'RatedPowerMW',Pr,'AvailabilityFactor',availability_factor);
pRaw = double(pRaw(:));

[capacityMW, xy, R] = ga_micrositing_max_capacity( ...
    feasible_region_info, dMin, ...
    'Mode','wake-nsga2', ...
    'TurbinePowerMW',Pr, ...
    'RotorDiameterM',D, ...
    'WindSpeedMps',U, ...
    'WindDirFromDeg',WD, ...
    'AvailabilityFactor',availability_factor, ...
    'AmbientTI',opt.AmbientTI, ...
    'PopulationSize',opt.PopulationSize, ...
    'MaxGenerations',opt.MaxGenerations, ...
    'UseBinnedWindStates',true, ...
    'WindSpeedBinWidthMps',2.0, ...
    'WindDirectionBinWidthDeg',45, ...
    'UseParallel',logical(opt.UseParallel), ...
    'NumWorkers',opt.NumWorkers, ...
    'PrecomputePairGeometry',true, ...
    'NSGARotorSamples',1, ...
    'FinalUseRepresentativeStates',true, ...
    'FinalWindSpeedBinWidthMps',1, ...
    'FinalWindDirectionBinWidthDeg',15, ...
    'FinalRotorSamples',9, ...
    'ParetoSelection',opt.ParetoSelection, ...
    'PlotFlag',false, ...
    'RandomSeed',opt.RandomSeed);

% ga_micrositing_max_capacity reconstructs the full hourly free/wake power
% using the same GWH242 curve and the same availability input. Therefore the
% Stage-II hourly p.u. series can be taken directly from that result.
pWake = double(R.wind_power_pu_wake_full(:));
pFreeCheck = double(R.wind_power_pu_free_full(:));
if numel(pWake) ~= numel(pRaw) || numel(pFreeCheck) ~= numel(pRaw)
    error('Hourly power length returned by RK-NSGA does not match 2023 wind series.');
end

% Numerical consistency check. This should be machine-precision small,
% because both paths use the same power-curve function and availability.
maxPowerMismatch = max(abs(pFreeCheck-pRaw));
if maxPowerMismatch > 1e-9
    error(['No-wake power reconstruction mismatch (max |delta p.u.| = %.3g). ', ...
        'Check power-curve/availability settings.'],maxPowerMismatch);
end

pWake = min(max(pWake,0),1);

out = struct();
out.capacity_MW = capacityMW;
out.turbine_count = size(xy,1);
out.turbine_xy_m = xy;
out.availability_factor = availability_factor;
out.raw_power_pu_2023 = pRaw;
out.raw_single_turbine_power_MW_2023 = PoneRaw_MW;
out.wake_adjusted_power_pu_2023 = pWake;
out.AEP_raw_MWh = capacityMW * sum(pRaw);
out.AEP_wake_adjusted_MWh = capacityMW * sum(pWake);
out.wake_loss_from_supplied_power_fraction = ...
    1 - out.AEP_wake_adjusted_MWh/max(out.AEP_raw_MWh,eps);
out.power_reconstruction_max_abs_pu_error = maxPowerMismatch;
out.nsga = R;
end
