function [keepGlobal,keepByGrid,info]=p44_spatial_repair(D,adjacency,keepGlobal,newGridIndices)
% Existing turbines are never removed. New grids are processed by B0 NPV/MW.
newGridIndices=newGridIndices(:);
score=cellfun(@(r)double(r.NPV_best)/double(r.wind_farm_capacity),D.records(newGridIndices));
[~,order]=sortrows([-score,double(D.grid_ids(newGridIndices))],[1 2]);
ordered=newGridIndices(order); removed=zeros(numel(ordered),1);
for k=1:numel(ordered)
    g=ordered(k);
    pts=double(D.point_start(g)):double(D.point_start(g+1))-1;
    assert(~any(keepGlobal(pts)),'Attempted to spatially repair an already-retained grid.');
    conflicting=adjacency(pts,:)*double(keepGlobal)>0;
    keepGlobal(pts)=~conflicting;
    removed(k)=sum(conflicting);
end
keepByGrid=cell(numel(D.records),1);
for g=1:numel(D.records)
    pts=double(D.point_start(g)):double(D.point_start(g+1))-1;
    keepByGrid{g}=keepGlobal(pts);
end
info=struct('processed_grid_indices',ordered,'processed_grid_ids',D.grid_ids(ordered), ...
    'deleted_counts',removed,'deleted_total',sum(removed), ...
    'zero_survivor_grid_indices',ordered(cellfun(@(v)~any(v),keepByGrid(ordered))), ...
    'rule','baseline NPV/MW descending, grid ID ascending; previously retained turbines protected');
end
