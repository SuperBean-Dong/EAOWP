function pairs = p44_cross_grid_pairs(xy,owner,radius)
%P44_CROSS_GRID_PAIRS Cross-owner pairs within a closed Euclidean radius.
% Base-MATLAB spatial bins replace scipy.spatial.cKDTree.query_pairs only
% during input preparation. IDs are one-based; same-owner pairs are omitted.
validateattributes(xy,{'numeric'},{'2d','real','finite'});
assert(size(xy,2)==2,'p44:PairCoordinates','Coordinates must be N by 2.');
owner=owner(:); n=size(xy,1);
assert(numel(owner)==n && all(isfinite(owner)),'p44:PairOwners','Invalid owners.');
validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
pairs=zeros(0,2);
if n<2, return; end
xy=double(xy);
bins=floor((xy-min(xy,[],1))/radius);
[~,order]=sortrows(bins,[1 2]); b=bins(order,:);
first=[1;find(any(diff(b,1,1),2))+1]; last=[first(2:end)-1;n];
occupied=b(first,:)+1;
lookup=sparse(occupied(:,1),occupied(:,2),(1:numel(first)).', ...
    max(occupied(:,1)),max(occupied(:,2)));
offsets=[0 0;0 1;1 -1;1 0;1 1];
chunks=cell(numel(first),1);
for g=1:numel(first)
    left=order(first(g):last(g)); local=zeros(0,2);
    for k=1:size(offsets,1)
        q=occupied(g,:)+offsets(k,:);
        if any(q<1) || q(1)>size(lookup,1) || q(2)>size(lookup,2), continue; end
        h=full(lookup(q(1),q(2))); if h==0, continue; end
        right=order(first(h):last(h));
        [ii,jj]=ndgrid(left,right); ii=ii(:); jj=jj(:);
        keep=owner(ii)~=owner(jj);
        if h==g, keep=keep & ii<jj; end
        ii=ii(keep); jj=jj(keep);
        delta=xy(ii,:)-xy(jj,:);
        keep=sum(delta.^2,2)<=radius^2;
        local=[local;sort([ii(keep),jj(keep)],2)]; %#ok<AGROW>
    end
    chunks{g}=local;
end
pairs=sortrows(vertcat(chunks{:}),[1 2]);
end
