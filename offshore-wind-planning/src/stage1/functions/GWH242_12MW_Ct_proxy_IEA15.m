function out = GWH242_12MW_Ct_proxy_IEA15()
%GWH242_12MW_CT_PROXY_IEA15  Ct(U) proxy for wake calculations.
%
% GWH242 public Ct data were not supplied with the project. Following the
% established model used by the micrositing code, this helper uses the
% publicly available IEA/NREL 15-MW reference-turbine operational thrust
% coefficient curve as a nondimensional proxy for the GWH242-12 MW rotor.
%
% The IEA reference curve has a 240-m rotor; Ct is nondimensional. The first
% point is extended from 3.0 to the GWH project cut-in speed 2.5 m/s so that
% the proxy is defined over the same active wind-speed range as the supplied
% GWH power curve. Outside [2.5,25] m/s, Ct is set to zero.
%
% Output fields:
%   .U_mps      wind-speed knots [m/s]
%   .Ct_values  thrust-coefficient knots [-]
%   .evaluate   vectorized function handle Ct = evaluate(U)

U = [ ...
     2.5, 3.0, 3.5, 4.0, 4.5, 4.8, 5.0, 5.2, 6.0, 6.2, ...
     6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.5, 8.0, 8.5, ...
     9.0, 9.5,10.0,10.2,10.5,10.6,10.7,10.8,10.9,11.0, ...
    11.2,11.5,11.8,12.0,13.0,14.0,15.0,17.5,20.0,22.5,25.0]';

Ct = [ ...
    0.820,0.820,0.801,0.808,0.822,0.823,0.823,0.831,0.835,0.834, ...
    0.832,0.829,0.826,0.820,0.816,0.810,0.807,0.805,0.805,0.804, ...
    0.804,0.804,0.803,0.802,0.802,0.769,0.699,0.672,0.635,0.607, ...
    0.549,0.501,0.461,0.426,0.321,0.251,0.201,0.126,0.085,0.061,0.046]';

out = struct();
out.U_mps = U;
out.Ct_values = Ct;
out.evaluate = @(u) evalCt(u,U,Ct);
out.proxy_turbine = 'IEA/NREL 15-MW reference turbine (240-m rotor)';
out.target_turbine = 'GWH242-12 MW (Ct proxy only)';
out.source_note = ['Operational Ct knots based on the public IEA/NREL 15-MW ', ...
    'reference-turbine power/thrust curve; 3.0-m/s Ct extended to 2.5 m/s.'];
end

function y = evalCt(u,U,Ct)
uq = double(u);
y = interp1(U,Ct,uq,'linear',0);
y(~isfinite(y)) = 0;
y = min(max(y,0),0.999999);
end
