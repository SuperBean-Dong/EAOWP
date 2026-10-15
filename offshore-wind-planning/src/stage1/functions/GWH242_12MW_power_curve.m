function [P_MW, p_pu, p_available_pu] = GWH242_12MW_power_curve(U_mps, varargin)
%GWH242_12MW_POWER_CURVE  GWH 242-12 MW segmented wind-speed/power curve.
%
% This implementation follows the supplied project function exactly for the
% wind-speed/power relation, but turbine availability is now an EXTERNAL
% function input rather than an internal fixed constant.
%
% Default aerodynamic power curve:
%   cut-in  = 2.5 m/s
%   rated   = 10.0 m/s
%   cut-out = 25.0 m/s
%
% Ideal power coefficient before availability:
%   p(U) = 0,                                      U < 2.5
%   p(U) = (U^3-2.5^3)/(10^3-2.5^3),             2.5 <= U <= 10
%   p(U) = 1,                                      10 < U <= 25
%   p(U) = 0,                                      U > 25
%
% Available single-turbine power:
%   P_MW = RatedPowerMW * p(U) * AvailabilityFactor
%
% OUTPUTS
%   P_MW            : available turbine output [MW]
%   p_pu            : ideal aerodynamic p.u. coefficient before availability
%   p_available_pu  : available p.u. output = p_pu*AvailabilityFactor
%
% EXAMPLE
%   P = GWH242_12MW_power_curve(U, ...
%       'RatedPowerMW',12,'AvailabilityFactor',0.95);
%
% AvailabilityFactor is intentionally provided by the caller so that the
% assumed turbine availability is explicit and auditable in the main model.

p = inputParser;
p.addRequired('U_mps', @isnumeric);
p.addParameter('RatedPowerMW', 12, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('AvailabilityFactor', 1.0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 1);
p.addParameter('CutInMps', 2.5, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('RatedWindMps', 10.0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('CutOutMps', 25.0, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.parse(U_mps, varargin{:});
prm = p.Results;

U = double(U_mps);
p_pu = zeros(size(U));

ci = prm.CutInMps;
vr = prm.RatedWindMps;
co = prm.CutOutMps;
if ~(ci < vr && vr < co)
    error('Require CutInMps < RatedWindMps < CutOutMps.');
end

ramp = isfinite(U) & U >= ci & U <= vr;
rated = isfinite(U) & U > vr & U <= co;

p_pu(ramp) = (U(ramp).^3 - ci^3) ./ (vr^3 - ci^3);
p_pu(rated) = 1;
p_pu = min(max(p_pu,0),1);

p_available_pu = p_pu .* prm.AvailabilityFactor;
P_MW = prm.RatedPowerMW .* p_available_pu;
end
