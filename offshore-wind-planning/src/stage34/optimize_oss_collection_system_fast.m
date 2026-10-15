function out = optimize_oss_collection_system_fast(turbine_capacity_MW,cable_capacity_MW, ...
    cable_cost_per_meter,export_cable_cost_per_meter,boundary_X,boundary_Y, ...
    daoan_x,daoan_y,wind_tubine_point,varargin)
%OPTIMIZE_OSS_COLLECTION_SYSTEM_FAST Budgeted, safely pruned original model.
% All original arcs, cable choices and the perforated OSS feasible region are
% retained. Each call clears YALMIP and uses ordinary optimize calls. No
% compiled model, warm start or parallel worker is reused. A feasible star
% layout guarantees a usable incumbent even when the time budget expires.
%
% Name-value options (seconds): BlockTimeLimitSeconds=60,
% GridTimeLimitSeconds=600, MIPGap=.05, SolverThreads=2, EnablePruning=true,
% StopWhenGlobalGapReached=true, Verbose=0.
% GridTimeLimitSeconds is a wall-clock admission budget: YALMIP compilation
% cannot be interrupted; each admitted Gurobi call has its own TimeLimit.

yalmip('clear');
gridTimer = tic;
ip = inputParser;
ip.addParameter('BlockTimeLimitSeconds',60,@positiveBudget);
ip.addParameter('GridTimeLimitSeconds',600,@positiveBudget);
ip.addParameter('MIPGap',.05,@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0&&v<1);
ip.addParameter('SolverThreads',2,@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1&&v==floor(v));
ip.addParameter('EnablePruning',true,@logicalScalar);
ip.addParameter('StopWhenGlobalGapReached',true,@logicalScalar);
ip.addParameter('Verbose',0,@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=0);
ip.parse(varargin{:});
cfg = ip.Results;

wt_xy = double(wind_tubine_point);
cable_capacity_MW = double(cable_capacity_MW(:).');
cable_cost_per_meter = double(cable_cost_per_meter(:).');
assert(size(wt_xy,2)==2 && ~isempty(wt_xy) && all(isfinite(wt_xy(:))), ...
    'Wind turbine coordinates must be a finite, nonempty n-by-2 matrix.');
assert(isscalar(turbine_capacity_MW)&&isfinite(turbine_capacity_MW)&&turbine_capacity_MW>0, ...
    'Turbine capacity must be positive.');
assert(~isempty(cable_capacity_MW)&&numel(cable_capacity_MW)==numel(cable_cost_per_meter) ...
    &&all(isfinite(cable_capacity_MW))&&all(cable_capacity_MW>=0) ...
    &&all(isfinite(cable_cost_per_meter))&&all(cable_cost_per_meter>=0), ...
    'Cable capacities and nonnegative costs must have matching finite entries.');
assert(isscalar(export_cable_cost_per_meter)&&isfinite(export_cable_cost_per_meter) ...
    &&export_cable_cost_per_meter>=0&&isscalar(daoan_x)&&isscalar(daoan_y) ...
    &&all(isfinite([daoan_x,daoan_y])), 'Invalid export cable cost or shore point.');
export_point = double([daoan_x;daoan_y]);
num_wt = size(wt_xy,1);
N = num_wt+1;
K = numel(cable_capacity_MW);
G = [zeros(1,2);wt_xy];

[rings,~] = splitBoundaryToRings(boundary_X,boundary_Y);
[feasible_ps,ring_depth] = buildFeasiblePolyshapeFromRings(rings);
TR = triangulation(feasible_ps);
tri_pts = TR.Points;
tri_conn = TR.ConnectivityList;
num_tri = size(tri_conn,1);
assert(num_tri>0,'The feasible region has no nondegenerate triangles.');
tri_poly_list = cell(num_tri,1);
tri_box = zeros(num_tri,4);
for tt=1:num_tri
    V = tri_pts(tri_conn(tt,:),:);
    tri_poly_list{tt} = polyshape(V(:,1),V(:,2),'Simplify',false);
    tri_box(tt,:) = [min(V(:,1)),max(V(:,1)),min(V(:,2)),max(V(:,2))];
end
feas_box = [min(tri_pts(:,1)),max(tri_pts(:,1)),min(tri_pts(:,2)),max(tri_pts(:,2))];
dist_matrix = zeros(N,N);
for j=2:N
    dist_matrix(2:end,j) = vecnorm(wt_xy-wt_xy(j-1,:),2,2);
end
nearest_other = inf(num_wt,1);
if num_wt>1
    D = dist_matrix(2:end,2:end);
    D(1:num_wt+1:end) = inf;
    nearest_other = min(D,[],1).';
end

% Precompute a rigorous lower bound for every nonempty rectangle before any
% expensive solve. Distance to its bounding rectangle underestimates distance
% to its true polygon intersection, including disconnected pieces and holes.
num_cols_div=7; num_rows_div=7; num_subregions=49;
blank = emptyBlock();
results = repmat(blank,num_subregions,1);
inter_poly_list = cell(num_subregions,1);
nonempty = false(num_subregions,1);
block_lb = nan(num_subregions,1);
min_internal_unit_cost = min(cable_cost_per_meter);
cell_size_x = (feas_box(2)-feas_box(1))/num_cols_div;
cell_size_y = (feas_box(4)-feas_box(3))/num_rows_div;
for r=1:num_rows_div
    for c=1:num_cols_div
        cid = (r-1)*num_cols_div+c;
        results(cid).cell_id=cid;
        box = [feas_box(1)+(c-1)*cell_size_x,feas_box(1)+c*cell_size_x, ...
            feas_box(3)+(r-1)*cell_size_y,feas_box(3)+r*cell_size_y];
        results(cid).cell_box = box;
        rect = polyshape([box(1),box(2),box(2),box(1)], ...
            [box(3),box(3),box(4),box(4)]);
        inter_poly_list{cid} = intersect(feasible_ps,rect);
        results(cid).area_m2 = area(inter_poly_list{cid});
        if results(cid).area_m2<=1e-8
            results(cid).status='empty_region';
            continue;
        end
        nonempty(cid)=true;
        d_root = distanceToBox(wt_xy,box);
        d_shore = distanceToBox(export_point.',box);
        raw_lb = export_cable_cost_per_meter*d_shore ...
            +min_internal_unit_cost*sum(min(d_root,nearest_other));
        block_lb(cid)=conservativeCostBound(raw_lb);
        results(cid).lower_bound_CNY=block_lb(cid);
        results(cid).status='not_solved';
    end
end
assert(any(nonempty),'The feasible region has no positive-area grid blocks.');

% A star is feasible for the original unrestricted collection model: one
% sufficiently large cable from each turbine to the OSS. Pick the best of
% triangle centroids using that same topology, rather than fixing the OSS
% for subsequent MISOCP solves.
eligible = find(cable_capacity_MW>=turbine_capacity_MW);
assert(~isempty(eligible),'No collection cable can carry one turbine.');
[~,cheapest] = min(cable_cost_per_meter(eligible));
star_type=eligible(cheapest);
centroids=(tri_pts(tri_conn(:,1),:)+tri_pts(tri_conn(:,2),:) ...
    +tri_pts(tri_conn(:,3),:))/3;
seedCosts=zeros(num_tri,1);
for tt=1:num_tri
    seedCosts(tt)=cable_cost_per_meter(star_type)*sum(vecnorm(wt_xy-centroids(tt,:),2,2)) ...
        +export_cable_cost_per_meter*norm(centroids(tt,:).'-export_point);
end
[~,seedTriangle]=min(seedCosts);
seedBeta=centroids(seedTriangle,:).';
seedX=zeros(N,N); seedY=zeros(N,N,K);
seedX(1,2:end)=1; seedY(1,2:end,star_type)=1;
seed = reconstructFeasibleCandidate(seedX,seedY,seedBeta,wt_xy,feasible_ps, ...
    feas_box,turbine_capacity_MW,cable_capacity_MW,cable_cost_per_meter, ...
    export_cable_cost_per_meter,export_point,dist_matrix);
assert(seed.valid,'Internal error: the fallback star failed physical validation.');
bestSol=seed.solution;
upper_bound=bestSol.cost;
initial_feasible_cost=upper_bound;
selected_block=pointBlock(seedBeta,feas_box,cell_size_x,cell_size_y);
assert(nonempty(selected_block),'The seed centroid was assigned to an empty block.');
results(selected_block)=recordCandidate(results(selected_block),bestSol,'fallback_star');
bestSource='fallback_star';
num_solver_calls=0; total_opt_time=0; budget_exhausted=false;
stop_reason='all_blocks_processed';
if toc(gridTimer)>=cfg.GridTimeLimitSeconds
    budget_exhausted=true;
    stop_reason='grid_time_budget_exhausted';
else

% Original full model, including all directed arcs and all cable types.
% Building an ordinary YALMIP model here retains the original solve workflow.
x=binvar(N,N,'full');
P=sdpvar(N,N,'full');
y=binvar(N,N,K,'full');
beta=sdpvar(2,1); t=sdpvar(N,1);
w=sdpvar(N,K,'full'); texport=sdpvar(1,1);
z_tri=binvar(num_tri,1,'full');
lambda_tri=sdpvar(num_tri,3,'full');
M=0;
for j=1:num_wt
    M=max(M,max(vecnorm(tri_pts-wt_xy(j,:),2,2)));
end
M=M+1;
obj=export_cable_cost_per_meter*texport;
for j=2:N
    for k=1:K
        obj=obj+cable_cost_per_meter(k)*w(j,k);
    end
end
for i=2:N
    for j=2:N
        if i~=j
            for k=1:K
                obj=obj+cable_cost_per_meter(k)*dist_matrix(i,j)*y(i,j,k);
            end
        end
    end
end
C=[sum(z_tri)==1,lambda_tri(:)>=0];
beta_expr=[0;0];
for tt=1:num_tri
    V=tri_pts(tri_conn(tt,:),:);
    C=[C,sum(lambda_tri(tt,:))==z_tri(tt)]; %#ok<AGROW>
    beta_expr=beta_expr+V'*lambda_tri(tt,:)';
end
C=[C,beta==beta_expr,t(1)==0,t>=0,w(1,:)==0];
for j=2:N
    C=[C,cone(beta-G(j,:)',t(j))]; %#ok<AGROW>
    for k=1:K
        C=[C,w(j,k)>=t(j)-(1-y(1,j,k))*M,w(j,k)<=t(j), ...
            w(j,k)<=M*y(1,j,k),w(j,k)>=0]; %#ok<AGROW>
    end
end
C=[C,cone(beta-export_point,texport),texport>=0];
for i=1:N
    for j=1:N
        if i==j
            C=[C,x(i,i)==0,P(i,i)==0]; %#ok<AGROW>
            for k=1:K
                C=[C,y(i,i,k)==0]; %#ok<AGROW>
            end
        else
            C=[C,x(i,j)==sum(y(i,j,:))]; %#ok<AGROW>
            cap_ij=0;
            for k=1:K
                cap_ij=cap_ij+cable_capacity_MW(k)*y(i,j,k);
            end
            C=[C,P(i,j)>=turbine_capacity_MW*x(i,j),P(i,j)<=cap_ij,P(i,j)>=0]; %#ok<AGROW>
        end
    end
end
C=[C,sum(x(:,1))==0];
for j=2:N
    C=[C,sum(x(:,j))==1,sum(P(j,:))+turbine_capacity_MW==sum(P(:,j))]; %#ok<AGROW>
end
C=[C,sum(P(1,2:end))==num_wt*turbine_capacity_MW];
ops=sdpsettings('solver','gurobi','verbose',cfg.Verbose);
% Pass only intended parameters: old YALMIP releases may carry Gurobi
% defaults that a newer installed solver no longer accepts.
ops.gurobi=struct('MIPGap',cfg.MIPGap,'Threads',cfg.SolverThreads,'TuneTimeLimit',0);
ops.savesolveroutput=1;

[~,order_local]=sort(block_lb(nonempty),'ascend');
active_cells=find(nonempty);
order=active_cells(order_local);
for kk=1:numel(order)
    cid=order(kk);
    if cfg.EnablePruning && block_lb(cid)>=upper_bound+costSafety(upper_bound)
        results(cid).status='pruned_by_safe_bound';
        results(cid).pruned=true;
        continue;
    end
    gap_now=relativeGlobalGap(upper_bound,block_lb(nonempty));
    if cfg.StopWhenGlobalGapReached && gap_now<=cfg.MIPGap
        stop_reason='global_gap_target_reached';
        break;
    end
    remaining=cfg.GridTimeLimitSeconds-toc(gridTimer);
    if remaining<=0
        budget_exhausted=true;
        stop_reason='grid_time_budget_exhausted';
        break;
    end
    box=results(cid).cell_box;
    inter_poly=inter_poly_list{cid};
    active_tri=false(num_tri,1);
    overlaps=tri_box(:,2)>=box(1)&tri_box(:,1)<=box(2) ...
        &tri_box(:,4)>=box(3)&tri_box(:,3)<=box(4);
    for tt=find(overlaps).'
        active_tri(tt)=area(intersect(tri_poly_list{tt},inter_poly))>1e-8;
    end
    if ~any(active_tri)
        results(cid).status='geometry_filter_unresolved';
        continue;
    end
    C_temp=[C,beta(1)>=box(1),beta(1)<=box(2),beta(2)>=box(3),beta(2)<=box(4), ...
        z_tri(~active_tri)==0,lambda_tri(~active_tri,:)==0];
    [ix,iy]=boundary(inter_poly);
    vertices=[ix(:),iy(:)]; vertices=vertices(all(isfinite(vertices),2),:);
    for j=2:N
        Mj=max(vecnorm(vertices-wt_xy(j-1,:),2,2))+1;
        for k=1:K
            C_temp=[C_temp,w(j,k)>=t(j)-(1-y(1,j,k))*Mj,w(j,k)<=Mj*y(1,j,k)]; %#ok<AGROW>
        end
    end
    remaining=cfg.GridTimeLimitSeconds-toc(gridTimer);
    if remaining<=0
        budget_exhausted=true;
        stop_reason='grid_time_budget_exhausted';
        break;
    end
    ops.gurobi.TimeLimit=min(cfg.BlockTimeLimitSeconds,remaining);
    solveTimer=tic;
    num_solver_calls=num_solver_calls+1;
    results(cid).solved=true;
    try
        result=optimize(C_temp,obj,ops);
        results(cid).elapsed_seconds=toc(solveTimer);
        results(cid).solver_problem=result.problem;
        results(cid).solver_info=result.info;
        results(cid).status='no_valid_incumbent';
        raw=struct();
        if isfield(result,'solveroutput')&&isfield(result.solveroutput,'result')
            raw=result.solveroutput.result;
        end
        if isfield(raw,'status'); results(cid).solver_status=raw.status; end
        if isfield(raw,'runtime'); results(cid).gurobi_runtime_seconds=raw.runtime; end
        candidate=struct('valid',false,'reason','No finite solver incumbent.');
        if isfield(raw,'x')&&isfield(raw,'objval')&&isfinite(raw.objval)
            candidate=reconstructFeasibleCandidate(value(x),value(y),value(beta), ...
                wt_xy,feasible_ps,box,turbine_capacity_MW,cable_capacity_MW, ...
                cable_cost_per_meter,export_cable_cost_per_meter,export_point,dist_matrix);
        end
        if candidate.valid
            results(cid)=recordCandidate(results(cid),candidate.solution,'validated_solver_topology');
            if strcmp(results(cid).solver_status,'OPTIMAL')
                results(cid).status='solver_gap_target_met';
            else
                results(cid).status='feasible_solver_limit';
            end
            if candidate.solution.cost<upper_bound
                upper_bound=candidate.solution.cost;
                bestSol=candidate.solution;
                selected_block=cid;
                bestSource='validated_solver_topology';
            end
        else
            results(cid).acceptance_note=candidate.reason;
        end
        % Bounds from numeric/infeasible/error statuses are not used. A
        % finite lower bound at normal solve limits remains useful even if
        % the solver found no feasible solution in this block.
        trusted_status={'OPTIMAL','TIME_LIMIT','NODE_LIMIT','ITERATION_LIMIT', ...
            'WORK_LIMIT','SOLUTION_LIMIT','INTERRUPTED','USER_OBJ_LIMIT'};
        if isfield(raw,'objbound')&&isfinite(raw.objbound) ...
                &&any(strcmp(results(cid).solver_status,trusted_status))
            solver_lb=conservativeCostBound(raw.objbound);
            if isfinite(results(cid).cost)
                if solver_lb>results(cid).cost+costSafety(results(cid).cost)
                    results(cid).acceptance_note=[results(cid).acceptance_note, ...
                        ' Solver bound exceeded the independently verified block incumbent; ignored.'];
                    solver_lb=NaN;
                else
                    solver_lb=min(solver_lb,results(cid).cost);
                end
            end
            if isfinite(solver_lb)
                results(cid).solver_lower_bound_CNY=solver_lb;
                block_lb(cid)=max(block_lb(cid),solver_lb);
            end
        end
        results(cid).lower_bound_CNY=block_lb(cid);
    catch ME
        results(cid).elapsed_seconds=toc(solveTimer);
        results(cid).status='solver_exception';
        results(cid).solver_info=ME.message;
        results(cid).acceptance_note=ME.message;
    end
    total_opt_time=total_opt_time+results(cid).elapsed_seconds;
end
end % The time budget may already have expired after geometry and fallback.

% Unvisited blocks keep their geometric lower bound. They are never omitted
% from the global gap, including on time-limit or solver-error exits.
final_gap=relativeGlobalGap(upper_bound,block_lb(nonempty));
if final_gap<=cfg.MIPGap
    optimality_status='gap_target_met';
elseif budget_exhausted
    optimality_status='feasible_time_limit';
else
    optimality_status='feasible_gap_not_met';
end
for cid=find(nonempty).'
    results(cid).lower_bound_CNY=block_lb(cid);
    results(cid).maymin=block_lb(cid);
    if strcmp(results(cid).status,'not_solved')
        if budget_exhausted; results(cid).status='not_solved_time_budget';
        elseif strcmp(stop_reason,'global_gap_target_reached'); results(cid).status='not_solved_global_gap_met';
        end
    end
end
out=struct();
out.oss_xy_m=bestSol.beta;
out.internal_cable_cost_CNY=bestSol.internal_cable_cost;
out.search_objective_cost_CNY=bestSol.cost;
out.oss_to_shore_distance_km=bestSol.oss_to_shore_distance_m/1000;
out.final_gap=final_gap;
out.total_optimization_time_s=total_opt_time;
out.selected_block_index=selected_block;
out.edge_detail=buildEdgeDetail(bestSol.x,bestSol.y,bestSol.P_MW, ...
    bestSol.t_oss_to_node_m,dist_matrix,cable_capacity_MW,cable_cost_per_meter);
out.solver_blocks=results;
out.feasible_polyshape=feasible_ps;
out.ring_depth=ring_depth;
out.triangulation=TR;
out.has_feasible_solution=true;
out.optimality_status=optimality_status;
out.budget_exhausted=budget_exhausted;
out.global_lower_bound_CNY=min(block_lb(nonempty));
out.diagnostics=struct('block_records',results,'initial_feasible_cost_CNY',initial_feasible_cost, ...
    'num_nonempty_blocks',sum(nonempty),'num_solver_calls',num_solver_calls, ...
    'num_pruned_blocks',sum([results.pruned]), ...
    'num_unsolved_blocks',sum(nonempty&~[results.solved].'&~[results.pruned].'), ...
    'wall_time_s',toc(gridTimer),'requested_gap',cfg.MIPGap,'stop_reason',stop_reason, ...
    'incumbent_source',bestSource,'settings',cfg, ...
    'grid_budget_semantics','Wall-clock admission budget; YALMIP compilation is not interruptible.');
end

function ok=positiveBudget(v)
ok=isnumeric(v)&&isscalar(v)&&~isnan(v)&&v>0;
end

function ok=logicalScalar(v)
ok=isscalar(v)&&(islogical(v)||(isnumeric(v)&&isfinite(v)&&ismember(v,[0,1])));
end

function B=emptyBlock()
B=struct('cell_id',NaN,'beta',[NaN;NaN],'cost',NaN,'internal_cable_cost',NaN,'maymin',NaN, ...
    'cell_box',nan(1,4),'status','not_solved','x',[],'y',[],'P_MW',[], ...
    't_oss_to_node_m',[],'oss_to_shore_distance_m',NaN,'area_m2',0, ...
    'lower_bound_CNY',NaN,'solver_lower_bound_CNY',NaN,'incumbent_cost_CNY',NaN, ...
    'solver_status','','solver_problem',NaN,'solver_info','','solved',false,'pruned',false, ...
    'has_feasible_incumbent',false,'elapsed_seconds',0,'gurobi_runtime_seconds',NaN, ...
    'acceptance_note','');
end

function B=recordCandidate(B,S,note)
if isfinite(B.cost)&&B.cost<=S.cost
    B.acceptance_note=[B.acceptance_note,'; ',note,' (earlier block incumbent retained)'];
    return;
end
names={'beta','cost','internal_cable_cost','x','y','P_MW', ...
    't_oss_to_node_m','oss_to_shore_distance_m'};
for k=1:numel(names); B.(names{k})=S.(names{k}); end
B.incumbent_cost_CNY=S.cost; B.has_feasible_incumbent=true; B.acceptance_note=note;
end

function d=distanceToBox(points,box)
dx=max(max(box(1)-points(:,1),0),points(:,1)-box(2));
dy=max(max(box(3)-points(:,2),0),points(:,2)-box(4));
d=hypot(dx,dy);
end

function t=costSafety(c)
t=max(1e-6,1e-9*abs(c));
end

function lb=conservativeCostBound(v)
lb=max(0,v-costSafety(v));
end

function gap=relativeGlobalGap(U,lower_bounds)
assert(all(isfinite(lower_bounds)),'Every nonempty block needs a finite safe lower bound.');
if U==0; gap=0; else; gap=max(0,(U-min(lower_bounds))/abs(U)); end
end

function cid=pointBlock(b,box,dx,dy)
c=min(7,max(1,floor((b(1)-box(1))/dx)+1));
r=min(7,max(1,floor((b(2)-box(3))/dy)+1));
cid=(r-1)*7+c;
end

function R=reconstructFeasibleCandidate(xv,yv,b,wt_xy,ps,box,q,capacities,costs, ...
    export_cost,shore,dist_matrix)
% Verify rounded topology, reconstruct exact subtree power flows and actual
% lengths, and evaluate a feasible original-model candidate independently of
% the solver's auxiliary distance/flow values. This also validates limit exits.
R=struct('valid',false,'reason','','solution',struct());
N=size(wt_xy,1)+1; K=numel(capacities);
if ~isequal(size(xv),[N,N])||numel(yv)~=N*N*K||numel(b)~=2 ...
        ||~all(isfinite(xv(:)))||~all(isfinite(yv(:)))||~all(isfinite(b(:)))
    R.reason='Missing/nonfinite primal topology or OSS coordinate.'; return;
end
b=b(:); yv=reshape(yv,N,N,K);
int_tol=1e-5;
if max(abs(xv(:)-round(xv(:))))>int_tol||max(abs(yv(:)-round(yv(:))))>int_tol ...
        ||any(xv(:)<-int_tol)||any(xv(:)>1+int_tol) ...
        ||any(yv(:)<-int_tol)||any(yv(:)>1+int_tol)
    R.reason='Topology/cable variables are not integral binary values.'; return;
end
xv=round(xv); yv=round(yv);
if any(any(xv~=sum(yv,3)))||any(diag(xv))||any(xv(:,1)) ...
        ||any(sum(xv(:,2:end),1)~=1)
    R.reason='Topology violates root, cable choice or unique-parent constraints.'; return;
end
% Do not project an infeasible OSS or silently accept a point inside a hole.
% isinterior(polyshape,...) includes boundary points in MATLAB.
if ~isinterior(ps,b(1),b(2))||b(1)<box(1)||b(1)>box(2)||b(2)<box(3)||b(2)>box(4)
    R.reason='OSS point is outside the true feasible region or its block.'; return;
end
parent=zeros(N,1);
for j=2:N; parent(j)=find(xv(:,j),1); end
flow=zeros(N,N);
for turbine=2:N
    j=turbine; visited=false(N,1);
    while j~=1
        if visited(j)||parent(j)<1
            R.reason='Topology contains a cycle or does not reach the OSS.'; return;
        end
        visited(j)=true;
        flow(parent(j),j)=flow(parent(j),j)+q;
        j=parent(j);
    end
end
length_to_root=[0;vecnorm(wt_xy-b.',2,2)];
internal=0;
for j=2:N
    i=parent(j); k=find(yv(i,j,:),1);
    if flow(i,j)>capacities(k)+1e-9
        R.reason='Reconstructed subtree power exceeds selected cable capacity.'; return;
    end
    if i==1; len=length_to_root(j); else; len=dist_matrix(i,j); end
    internal=internal+costs(k)*len;
end
shore_distance=norm(b-shore);
S=struct('beta',b,'cost',internal+export_cost*shore_distance, ...
    'internal_cable_cost',internal,'x',xv,'y',yv,'P_MW',flow, ...
    't_oss_to_node_m',length_to_root,'oss_to_shore_distance_m',shore_distance);
R.valid=true; R.reason='Independently validated physical topology and exact flows.'; R.solution=S;
end

function E=buildEdgeDetail(xv,yv,Pv,tv,dist_matrix,capacities,costs)
E=struct('from_node',{},'to_node',{},'cable_type',{},'flow_MW',{}, ...
    'length_m',{},'cable_capacity_MW',{},'unit_cost_CNY_per_m',{});
N=size(xv,1);
for j=2:N
    i=find(xv(:,j),1); k=find(yv(i,j,:),1);
    if i==1; len=tv(j); else; len=dist_matrix(i,j); end
    E(end+1)=struct('from_node',j,'to_node',i,'cable_type',k,'flow_MW',Pv(i,j), ...
        'length_m',len,'cable_capacity_MW',capacities(k),'unit_cost_CNY_per_m',costs(k)); %#ok<AGROW>
end
end

function [rings,areas]=splitBoundaryToRings(X,Y)
X=X(:); Y=Y(:);
assert(numel(X)==numel(Y),'Boundary X/Y lengths differ.');
split_idx=[0;find(isnan(X)|isnan(Y));numel(X)+1];
rings={}; areas=[];
for s=1:numel(split_idx)-1
    i1=split_idx(s)+1; i2=split_idx(s+1)-1;
    if i2<i1; continue; end
    xy=[X(i1:i2),Y(i1:i2)]; xy=xy(all(isfinite(xy),2),:);
    xy=removeConsecutiveDuplicatePoints(xy);
    if size(xy,1)<3; continue; end
    if norm(xy(1,:)-xy(end,:))>1e-10; xy=[xy;xy(1,:)]; end %#ok<AGROW>
    ar=abs(polyarea(xy(:,1),xy(:,2)));
    if ar<=1e-8; continue; end
    rings{end+1}=xy; areas(end+1,1)=ar; %#ok<AGROW>
end
end

function xy=removeConsecutiveDuplicatePoints(xy)
if isempty(xy); return; end
keep=true(size(xy,1),1);
for i=2:size(xy,1)
    if norm(xy(i,:)-xy(i-1,:))<=1e-12; keep(i)=false; end
end
xy=xy(keep,:);
if size(xy,1)>=2&&norm(xy(1,:)-xy(end,:))<=1e-12; xy=xy(1:end-1,:); end
end

function [feasible_ps,depth]=buildFeasiblePolyshapeFromRings(rings)
n=numel(rings); ring_ps=cell(n,1); ring_area=zeros(n,1);
for i=1:n
    xy=rings{i}; ring_ps{i}=polyshape(xy(:,1),xy(:,2),'Simplify',false);
    ring_area(i)=area(ring_ps{i});
end
depth=zeros(n,1);
for i=1:n
    for j=1:n
        if i==j||ring_area(j)<=ring_area(i)+1e-8; continue; end
        inter_ij=intersect(ring_ps{i},ring_ps{j});
        if abs(area(inter_ij)-ring_area(i))<=max(1e-6,1e-8*ring_area(i))
            depth(i)=depth(i)+1;
        end
    end
end
[~,order]=sort(ring_area,'descend'); feasible_ps=polyshape();
for kk=1:numel(order)
    i=order(kk);
    if mod(depth(i),2)==0; feasible_ps=union(feasible_ps,ring_ps{i});
    else; feasible_ps=subtract(feasible_ps,ring_ps{i}); end
end
feasible_ps=simplify(feasible_ps);
end
