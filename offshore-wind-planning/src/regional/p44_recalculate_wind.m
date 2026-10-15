function out = p44_recalculate_wind(r,keepMask,wind,cfg)
%P44_RECALCULATE_WIND Re-evaluate an unchanged-coordinate turbine subset.
% Matches the production final evaluation: 1 m/s, 15-degree representative
% states, 9 rotor samples, and exact hourly free-stream power. Only turbines
% within this original grid participate in wake evaluation. No NSGA is run.
% wind requires wind_speed_140m_mps and wind_direction_from_deg vectors.
if nargin<4; cfg=struct(); end
xy=double(r.turbine_xy_m);
keep=p44_keep_mask(keepMask,size(xy,1));
U=double(wind.wind_speed_140m_mps(:));
WD=mod(double(wind.wind_direction_from_deg(:)),360);
assert(numel(U)==numel(WD) && ~isempty(U) && all(isfinite(U)) ...
    && all(U>=0) && all(isfinite(WD)), ...
    'p44:WindInput','Invalid or unequal-length wind inputs.');
Pr=12; D=242;
availability=0.95;
if isfield(r,'availability_factor') && ~isempty(r.availability_factor)
    availability=double(r.availability_factor);
end
if isfield(cfg,'turbine')
    Pr=p44_option(cfg.turbine,'rated_power_MW',Pr);
    D=p44_option(cfg.turbine,'rotor_diameter_m',D);
    availability=p44_option(cfg.turbine,'availability_factor',availability);
end
TI=p44_option(cfg,'ambient_TI',0.0555);
dU=p44_option(cfg,'wind_speed_bin_width_mps',1);
dDir=p44_option(cfg,'wind_direction_bin_width_deg',15);
samples=p44_option(cfg,'rotor_samples',9);
assert(Pr>0 && D>0 && availability>=0 && availability<=1 ...
    && dU>0 && dDir>0 && dDir<=180 && samples>=1, ...
    'p44:WindParameters','Invalid wake model parameters.');
n=sum(keep); capacity=n*Pr;
out=struct('capacity_MW',capacity,'wind_farm_capacity',capacity, ...
    'turbine_count',n,'turbine_xy_m',xy(keep,:), ...
    'kept_original_indices',find(keep), ...
    'Pwind_MW',zeros(numel(U),1), ...
    'wake_adjusted_power_pu_2023',zeros(numel(U),1), ...
    'AEP_raw_MWh',0,'AEP_wake_adjusted_MWh',0, ...
    'wake_loss_fraction',0,'actual_wind_peak_MW',0, ...
    'representative_num_states',0,'ambient_TI',TI, ...
    'wind_speed_bin_width_mps',dU,'wind_direction_bin_width_deg',dDir, ...
    'rotor_samples',samples,'availability_factor',availability, ...
    'wake_scope','within_original_grid', ...
    'method','original final representative-state evaluation; coordinates fixed');
if n==0; return; end
C=GWH242_12MW_Ct_proxy_IEA15();
turb=struct('rotor_diameter_m',D,'rated_power_MW',Pr, ...
    'CtFcn',C.evaluate,'PowerFcn',@(v)GWH242_12MW_power_curve(v, ...
    'RatedPowerMW',Pr,'AvailabilityFactor',availability));
speedBin=floor(U./dU)+1;
nDir=ceil(360/dDir);
dirBin=min(floor(WD./dDir)+1,nDir);
pairID=(speedBin-1)*nDir+dirBin;
[uID,~,grp]=unique(pairID,'sorted');
hours=accumarray(grp,1);
Umean=accumarray(grp,U,[],@mean);
dirIdx=mod(uID-1,nDir)+1;
WDcenter=mod((double(dirIdx)-0.5).*dDir,360);
rep=bastankhah_gaussian_wake(xy(keep,:),Umean,WDcenter,turb, ...
    'AmbientTI',TI,'StateWeightsHours',hours,'RotorSamples',samples, ...
    'StoreTurbineStates',false);
eta=ones(size(rep.farm_power_free_MW(:)));
valid=rep.farm_power_free_MW(:)>1e-12;
eta(valid)=rep.farm_power_wake_MW(valid)./rep.farm_power_free_MW(valid);
eta=min(max(eta,0),1);
Pone=double(turb.PowerFcn(U));
Pone(~isfinite(Pone))=0;
Pone=min(max(Pone(:),0),Pr);
farmFree=n.*Pone;
out.Pwind_MW=farmFree.*eta(grp);
out.wake_adjusted_power_pu_2023=min(max(out.Pwind_MW/capacity,0),1);
out.AEP_raw_MWh=sum(farmFree);
out.AEP_wake_adjusted_MWh=sum(out.Pwind_MW);
out.wake_loss_fraction=1-out.AEP_wake_adjusted_MWh/max(out.AEP_raw_MWh,eps);
out.actual_wind_peak_MW=max(out.Pwind_MW);
out.representative_num_states=numel(uID);
assert(all(isfinite(out.Pwind_MW)) && all(out.Pwind_MW>=0) ...
    && all(out.Pwind_MW<=capacity+1e-8), ...
    'p44:WindOutput','Invalid reconstructed wind output.');
end
