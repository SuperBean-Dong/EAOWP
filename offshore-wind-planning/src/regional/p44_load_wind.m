function wind = p44_load_wind(gridID,cfg)
% Read two hourly series only; do not load the regional wind database.
persistent cachedFile ids
if isempty(cachedFile) || ~strcmp(cachedFile,cfg.WindFile)
    ids=double(h5read(cfg.WindFile,'/wind2023/grid_index'));
    ids=ids(:); cachedFile=cfg.WindFile;
end
column=find(ids==gridID); assert(isscalar(column),'Missing/duplicate wind grid %d.',gridID);
wind.wind_speed_140m_mps=readColumn(cfg.WindFile,'/wind2023/wind_speed_140m_mps',column);
wind.wind_direction_from_deg=readColumn(cfg.WindFile,'/wind2023/wind_direction_from_deg',column);
wind.grid_index=gridID;
end
function v=readColumn(file,path,column)
info=h5info(file,path); shape=info.Dataspace.Size;
if shape(1)==8760
    v=h5read(file,path,[1 column],[8760 1]);
elseif shape(2)==8760
    v=h5read(file,path,[column 1],[1 8760]);
else
    error('p44:WindShape','Unexpected wind dimensions.');
end
v=double(v(:)); assert(numel(v)==8760 && all(isfinite(v)));
end
