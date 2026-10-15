function out = p44_repair_collection(r,keepMask,cfg)
%P44_REPAIR_COLLECTION Prune/rescale an original radial collection tree.
% Original node 1 is OSS; node t+1 is original turbine t. A deleted turbine
% location remains a junction if retained generation still traverses it.
% Original route lengths, parent nodes and OSS coordinates are unchanged.
if nargin<3; cfg=struct(); end
n=size(r.turbine_xy_m,1);
keep=p44_keep_mask(keepMask,n);
E=r.collection_edge_detail;
assert(isstruct(E) && numel(E)==n, ...
    'p44:CollectionTree','Expected one original cable edge per turbine.');
Pr=12;
if isfield(cfg,'turbine'); Pr=p44_option(cfg.turbine,'rated_power_MW',Pr); end
capacities=double(p44_option(cfg,'cable_capacity_MW',[40 50 63]));
costs=double(p44_option(cfg,'cable_unit_cost_CNY_per_m',[423 474 554]*8.39));
capacities=capacities(:); costs=costs(:);
assert(Pr>0 && numel(capacities)==numel(costs) ...
    && all(isfinite(capacities)) && all(capacities>0) ...
    && all(isfinite(costs)) && all(costs>=0), ...
    'p44:CableParameters','Invalid cable capacities or unit costs.');
from=double([E.from_node]).'; parentTo=double([E.to_node]).';
len=double([E.length_m]).';
assert(isequal(sort(from),(2:n+1).') && all(parentTo>=1) ...
    && all(parentTo<=n+1) && all(parentTo==round(parentTo)) ...
    && all(from~=parentTo) && all(isfinite(len)) && all(len>=0), ...
    'p44:CollectionTree','Original topology has invalid nodes or route lengths.');
parent=zeros(n+1,1); edgeOfNode=zeros(n+1,1);
parent(from)=parentTo; edgeOfNode(from)=(1:n).';
% Validate every original node, including branches that will be pruned.
for t=2:n+1
    node=t; visited=false(n+1,1);
    while node~=1
        assert(~visited(node) && parent(node)>=1, ...
            'p44:CollectionTree','Original topology contains a cycle or detached node.');
        visited(node)=true; node=parent(node);
    end
end
flow=zeros(n,1);
for t=find(keep).'
    node=t+1;
    while node~=1
        j=edgeOfNode(node); flow(j)=flow(j)+Pr; node=parent(node);
    end
end
active=flow>0;
junction=find(~keep & active(edgeOfNode(2:n+1)));
edges=E(active);
activeOriginal=find(active);
for j=1:numel(edges)
    old=activeOriginal(j);
    edges(j).flow_MW=flow(old);
    if all(keep)
        % An unchanged project retains exactly its frozen baseline design.
        assert(edges(j).cable_capacity_MW+1e-6>=flow(old), ...
            'p44:CableCapacity','Original selected cable cannot carry subtree generation.');
    else
        feasible=find(capacities+1e-8>=flow(old));
        assert(~isempty(feasible),'p44:CableCapacity', ...
            'Retained generation exceeds every available cable capacity.');
        [~,at]=min(costs(feasible)); k=feasible(at);
        edges(j).cable_type=k;
        edges(j).cable_capacity_MW=capacities(k);
        edges(j).unit_cost_CNY_per_m=costs(k);
    end
end
internal=sum(double([edges.length_m]).*double([edges.unit_cost_CNY_per_m]));
if all(keep) && isfield(r,'internal_cable_cost')
    assert(abs(internal-double(r.internal_cable_cost)) ...
        <=max(1,1e-7*abs(double(r.internal_cable_cost))), ...
        'p44:BaselineCollectionCost','Original edges and saved cable cost disagree.');
    internal=double(r.internal_cable_cost);
end
out=struct('internal_cable_cost_CNY',internal,'edges',edges, ...
    'junction_original_indices',junction(:),'junction_count',numel(junction), ...
    'kept_original_indices',find(keep), ...
    'pruned_original_edge_indices',find(~active), ...
    'retained_original_edge_indices',activeOriginal, ...
    'route_policy','original OSS, parent nodes and cable route lengths fixed', ...
    'junction_policy','deleted ancestor locations retain connection facilities');
if isfield(r,'beta'); out.oss_xy_m=double(r.beta); end
end
