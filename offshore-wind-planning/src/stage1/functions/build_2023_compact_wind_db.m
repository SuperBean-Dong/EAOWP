function wind2023 = build_2023_compact_wind_db(shpFile,speedFile,directionFile,outFile, ...
    ratedPowerMW,availabilityFactor,varargin)
%BUILD_2023_COMPACT_WIND_DB Build compact 2023 wind data for feasible grids.
%
% IMPORTANT REVISION
%   Wind power is NOT read from south_sea_wind_power_pu_2023.mat.
%   It is recalculated directly from the 2023 140-m wind-speed series using
%   GWH242_12MW_power_curve.m.
%
%   The resulting stored field power_pu is the AVAILABLE no-wake farm power
%   coefficient, i.e. single-turbine available output divided by rated power:
%
%       power_pu = P_MW / RatedPowerMW
%                = p_curve(U) * AvailabilityFactor.
%
%   Thus the same availability factor is consistently inherited by Stage I
%   AEP calculation and Stage II hourly wind-power input.
%
% The output is organized as [hour x feasible-grid-index] and saved as
% -v7.3. Normal optimization runs load this compact database fully into RAM.
%
% Example:
%   build_2023_compact_wind_db('feasible_final.shp',speedFile,directionFile, ...
%       'wind_2023_feasible_compact_from_speed.mat',12,0.95);

validateattributes(ratedPowerMW,{'numeric'}, ...
    {'scalar','real','finite','>',0},mfilename,'ratedPowerMW');
validateattributes(availabilityFactor,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<=',1},mfilename,'availabilityFactor');

p = inputParser;
p.addParameter('CoordinateToleranceDeg',0.02, ...
    @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Verbose',true, ...
    @(x)islogical(x)||(isnumeric(x)&&isscalar(x)));
p.parse(varargin{:});
opt = p.Results;

assert(isfile(shpFile),'Cannot find shapefile: %s',shpFile);
assert(isfile(speedFile),'Cannot find 2023 140-m wind-speed file: %s',speedFile);
assert(isfile(directionFile),'Cannot find 2023 140-m wind-direction file: %s',directionFile);

S = shaperead(shpFile);
nGrid = numel(S);
gridLon = nan(nGrid,1);
gridLat = nan(nGrid,1);
for i = 1:nGrid
    gridLon(i) = double(S(i).X1);
    gridLat(i) = double(S(i).Y1);
end

if logical(opt.Verbose)
    fprintf('Building compact 2023 wind DB for %d feasible grids.\n',nGrid);
    fprintf('Wind power will be recalculated from 140-m wind speed.\n');
    fprintf('Rated power = %.3f MW; availability factor = %.4f.\n', ...
        ratedPowerMW,availabilityFactor);
end

% --------------------------- 1) 140-m speed -----------------------------
Ssrc = loadGeoCube(speedFile, ...
    {'wind_speed_140m_southchina_sea','wind_speed_140m', ...
     'wind_speed_southchina_sea','wind_speed'}, ...
    '140-m wind speed',[],[]);
[wind_speed_140m_mps,speedIndex] = extractAllFeasibleSeries( ...
    Ssrc.data,Ssrc.lat,Ssrc.lon,gridLon,gridLat, ...
    opt.CoordinateToleranceDeg,'140-m wind speed');
numHours = size(wind_speed_140m_mps,1);
speedVarName = Ssrc.data_name;
speedLatName = Ssrc.lat_name;
speedLonName = Ssrc.lon_name;
speedCoordLat = Ssrc.lat;
speedCoordLon = Ssrc.lon;
clear Ssrc

% ------------------------- 2) 140-m direction ---------------------------
Dsrc = loadGeoCube(directionFile, ...
    {'wind_direction_140m_southchina_sea','wind_direction_140m', ...
     'wind_direction_southchina_sea','wind_direction'}, ...
    '140-m wind direction',speedCoordLat,speedCoordLon);
[wind_direction_from_deg,directionIndex] = extractAllFeasibleSeries( ...
    Dsrc.data,Dsrc.lat,Dsrc.lon,gridLon,gridLat, ...
    opt.CoordinateToleranceDeg,'140-m wind direction');
if size(wind_direction_from_deg,1) ~= numHours
    error('2023 wind-direction series has %d hours but speed data has %d.', ...
        size(wind_direction_from_deg,1),numHours);
end
directionVarName = Dsrc.data_name;
directionLatName = Dsrc.lat_name;
directionLonName = Dsrc.lon_name;
clear Dsrc speedCoordLat speedCoordLon

% -------------------- 3) Recalculate no-wake power ---------------------
wind_speed_140m_mps = double(wind_speed_140m_mps);
wind_direction_from_deg = mod(double(wind_direction_from_deg),360);

if any(~isfinite(wind_speed_140m_mps(:))) || any(wind_speed_140m_mps(:)<0)
    error('Compact 140-m wind-speed data contains invalid values.');
end
if any(~isfinite(wind_direction_from_deg(:)))
    error('Compact 140-m wind-direction data contains invalid values.');
end

% GWH242_12MW_power_curve follows the supplied project power curve:
% cut-in 2.5 m/s, rated 10 m/s, cut-out 25 m/s. Availability is supplied
% explicitly here rather than being hard-coded inside the function.
P_one_MW = GWH242_12MW_power_curve(wind_speed_140m_mps, ...
    'RatedPowerMW',ratedPowerMW, ...
    'AvailabilityFactor',availabilityFactor);
power_pu = P_one_MW ./ ratedPowerMW;
clear P_one_MW
power_pu = min(max(double(power_pu),0),1);

% ----------------------------- output -----------------------------------
wind2023 = struct();
wind2023.data_year = 2023;
wind2023.num_hours = numHours;
wind2023.grid_count = nGrid;
wind2023.grid_index = (1:nGrid).';
wind2023.grid_lon = gridLon;
wind2023.grid_lat = gridLat;
wind2023.power_pu = power_pu;                              % T x nGrid
wind2023.wind_speed_140m_mps = wind_speed_140m_mps;        % T x nGrid
wind2023.wind_direction_from_deg = wind_direction_from_deg;% T x nGrid
wind2023.source_speed_row_col = speedIndex;
wind2023.source_direction_row_col = directionIndex;
wind2023.source_files = struct('speed140',speedFile,'direction140',directionFile);
wind2023.source_variables = struct( ...
    'speed140',speedVarName,'speed_lat',speedLatName,'speed_lon',speedLonName, ...
    'direction140',directionVarName,'direction_lat',directionLatName, ...
    'direction_lon',directionLonName);
wind2023.rated_power_MW = ratedPowerMW;
wind2023.availability_factor = availabilityFactor;
wind2023.power_curve = struct( ...
    'function','GWH242_12MW_power_curve', ...
    'cut_in_mps',2.5,'rated_wind_mps',10.0,'cut_out_mps',25.0, ...
    'description',['p=0 below cut-in; cubic normalized ramp from 2.5 to 10 m/s; ', ...
                   'p=1 above rated through 25 m/s; available power multiplied by availability.']);
wind2023.power_provenance = [ ...
    'Recalculated from the supplied 2023 140-m wind-speed series using ', ...
    'GWH242_12MW_power_curve; power_pu includes AvailabilityFactor.'];
wind2023.direction_convention = ...
    'meteorological FROM north clockwise: 0=N, 90=E, 180=S, 270=W';
wind2023.coordinate_tolerance_deg = opt.CoordinateToleranceDeg;
wind2023.shapefile = shpFile;
wind2023.created_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

outDir = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
save(outFile,'wind2023','-v7.3');

if logical(opt.Verbose)
    info = dir(outFile);
    fprintf('Compact 2023 wind DB saved: %s\n',outFile);
    fprintf('  dimensions: %d h x %d feasible grids x 3 data fields\n',numHours,nGrid);
    fprintf('  power source: recalculated from 140-m wind speed\n');
    fprintf('  file size: %.2f GB\n',info.bytes/1024^3);
end
end

function src = loadGeoCube(fileName,preferredDataNames,role,fallbackLat,fallbackLon)
W = whos('-file',fileName);
names = {W.name};

dataName = findPreferredName(names,preferredDataNames);
if isempty(dataName)
    isNum = arrayfun(@(x)ismember(x.class,{'double','single'}),W);
    isCube = arrayfun(@(x)numel(x.size)>=3 && prod(x.size)>1,W);
    cand = find(isNum & isCube);
    if isempty(cand)
        error('No numeric 3-D data variable found for %s in %s.',role,fileName);
    end
    [~,k] = max([W(cand).bytes]);
    dataName = W(cand(k)).name;
end

latName = findPreferredName(names,{'lat','latitude','LAT','Latitude','grid_lat','lat_grid'});
lonName = findPreferredName(names,{'lon','longitude','LON','Longitude','grid_lon','lon_grid'});

vars = {dataName};
if ~isempty(latName); vars{end+1}=latName; end %#ok<AGROW>
if ~isempty(lonName); vars{end+1}=lonName; end %#ok<AGROW>
A = load(fileName,vars{:});
data = A.(dataName);
if ndims(data) < 3
    error('%s variable %s in %s is not a 3-D [space x space x time] array.', ...
        role,dataName,fileName);
end

if ~isempty(latName) && ~isempty(lonName)
    lat = double(A.(latName));
    lon = double(A.(lonName));
else
    if isempty(fallbackLat) || isempty(fallbackLon)
        error(['Cannot find latitude/longitude variables in %s. ', ...
            'Expected names such as lat/lon or latitude/longitude.'],fileName);
    end
    lat = double(fallbackLat);
    lon = double(fallbackLon);
    latName = '<reused from wind-speed file>';
    lonName = '<reused from wind-speed file>';
end

src = struct('data',data,'lat',lat,'lon',lon, ...
    'data_name',dataName,'lat_name',latName,'lon_name',lonName);
end

function name = findPreferredName(names,preferred)
name = '';
for k = 1:numel(preferred)
    j = find(strcmpi(names,preferred{k}),1);
    if ~isempty(j)
        name = names{j};
        return;
    end
end
end

function [X,idxRC] = extractAllFeasibleSeries(cube,lat,lon,gridLon,gridLat,tol,role)
sz = size(cube);
T = sz(3);
nGrid = numel(gridLon);
X = zeros(T,nGrid,'double');
idxRC = zeros(nGrid,2);

for i = 1:nGrid
    [r,c] = locateGridIndex(lat,lon,sz(1:2),gridLon(i),gridLat(i),tol);
    if isempty(r) || isempty(c)
        error('Cannot locate feasible grid %d (lon=%.6f, lat=%.6f) in %s source grid.', ...
            i,gridLon(i),gridLat(i),role);
    end
    v = squeeze(cube(r,c,:));
    if numel(v) ~= T
        error('Unexpected time-series length for feasible grid %d in %s.',i,role);
    end
    X(:,i) = double(v(:));
    idxRC(i,:) = [r,c];
end
end

function [r,c] = locateGridIndex(lat,lon,spatialSize,targetLon,targetLat,tol)
r = []; c = [];

% Case A: 2-D coordinate grids. Preserve the indexing convention used in
% the original optimization code: row is identified from longitude and
% column is identified from latitude.
if ~isvector(lat) && ~isvector(lon)
    if isequal(size(lat),spatialSize) && isequal(size(lon),spatialSize)
        [~,cLat] = find(abs(lat-targetLat)<=tol,1,'first');
        [rLon,~] = find(abs(lon-targetLon)<=tol,1,'first');
        if ~isempty(rLon) && ~isempty(cLat)
            r = rLon; c = cLat;
            return;
        end
    elseif isequal(fliplr(size(lat)),spatialSize) && ...
            isequal(fliplr(size(lon)),spatialSize)
        latT = lat.'; lonT = lon.';
        [~,cLat] = find(abs(latT-targetLat)<=tol,1,'first');
        [rLon,~] = find(abs(lonT-targetLon)<=tol,1,'first');
        if ~isempty(rLon) && ~isempty(cLat)
            r = rLon; c = cLat;
            return;
        end
    end
end

% Case B: coordinate vectors. Infer which spatial dimension is lon/lat.
if isvector(lat) && isvector(lon)
    latv = lat(:); lonv = lon(:);
    ilat = find(abs(latv-targetLat)<=tol,1,'first');
    ilon = find(abs(lonv-targetLon)<=tol,1,'first');
    if isempty(ilat) || isempty(ilon)
        return;
    end
    if numel(lonv)==spatialSize(1) && numel(latv)==spatialSize(2)
        r = ilon; c = ilat;
        return;
    elseif numel(latv)==spatialSize(1) && numel(lonv)==spatialSize(2)
        r = ilat; c = ilon;
        return;
    end
end

% Case C: same-size 2-D coordinates but unusual shape. Use nearest joint
% coordinate only as a structural fallback, without interpolation.
if numel(lat)==prod(spatialSize) && numel(lon)==prod(spatialSize)
    lat2 = reshape(lat,spatialSize);
    lon2 = reshape(lon,spatialSize);
    d2 = (lat2-targetLat).^2 + (lon2-targetLon).^2;
    [dmin,lin] = min(d2(:));
    if isfinite(dmin) && sqrt(dmin) <= sqrt(2)*tol
        [r,c] = ind2sub(spatialSize,lin);
    end
end
end
