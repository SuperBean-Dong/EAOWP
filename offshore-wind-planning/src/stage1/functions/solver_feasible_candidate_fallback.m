function [candidatePoints, info] = solver_feasible_candidate_fallback(regionPS, dMin)
%SOLVER_FEASIBLE_CANDIDATE_FALLBACK Repair an empty original lattice only.
% Triangulation uses the real polyshape, including disconnected components
% and holes. Longest-edge subdivision bounds each triangle diameter by
% dMin/2. Every candidate is the solution of a linear program maximizing
% its distance from the triangle edges, then checked against regionPS.
% These are candidate locations, not turbine layouts: the unchanged GA,
% hard-spacing conflict matrix and wake model still select/evaluate layouts.
%
% No random numbers, polygon convex hulls, rounded coordinates, constant
% capacity substitutes or zero-result substitutes are used here.

if ~isa(regionPS,'polyshape') || area(regionPS) <= 0
    error('EmptyCandidateRepair:InvalidRegion', ...
        'The fallback requires a positive-area, real feasible polyshape.');
end
if ~isscalar(dMin) || ~isfinite(dMin) || dMin <= 0
    error('EmptyCandidateRepair:InvalidSpacing','dMin must be positive.');
end

solver = chooseLPSolver();
mesh = triangulation(regionPS);
conn = mesh.ConnectivityList;
vertices = mesh.Points;
if isempty(conn)
    error('EmptyCandidateRepair:NoTriangles', ...
        'The positive-area feasible region did not produce any triangles.');
end

% Each subdivision stays within its original triangle, so it cannot bridge
% a concavity, disconnected components, or a polygon hole.
targetEdge = dMin/2;
pending = cell(size(conn,1),1);
for ii = 1:size(conn,1)
    pending{ii} = vertices(conn(ii,:),:);
end
head = 1;
candidatePoints = zeros(0,2);
radiusM = zeros(0,1);
refinedCount = 0;
degenerateCount = 0;
maxObservedEdge = 0;
maxLPViolation = 0;
% Fail explicitly if malformed geometry would produce an excessive library;
% never silently truncate candidate coverage or fabricate an optimization.
maxRefinedTriangles = 50000;
while head <= numel(pending)
    tri = pending{head};
    pending{head} = [];
    head = head+1;
    edgeIdx = [1 2;2 3;3 1];
    edgeVec = tri(edgeIdx(:,2),:) - tri(edgeIdx(:,1),:);
    edgeLength = sqrt(sum(edgeVec.^2,2));
    [longest, edge] = max(edgeLength);
    twiceArea = cross2d(tri(2,:)-tri(1,:),tri(3,:)-tri(1,:));
    if longest == 0 || abs(twiceArea) <= 64*eps(longest^2)
        degenerateCount = degenerateCount+1;
        continue;
    end
    if longest > targetEdge*(1+1e-12)
        ia = edgeIdx(edge,1);
        ib = edgeIdx(edge,2);
        ic = setdiff(1:3,[ia ib]);
        midpoint = (tri(ia,:)+tri(ib,:))/2;
        pending{end+1,1} = [tri(ia,:);midpoint;tri(ic,:)]; %#ok<AGROW>
        pending{end+1,1} = [midpoint;tri(ib,:);tri(ic,:)]; %#ok<AGROW>
        if numel(pending)-head+1+refinedCount > maxRefinedTriangles
            error('EmptyCandidateRepair:LibraryTooLarge', ...
                'Triangle refinement exceeded the explicit safety limit.');
        end
        continue;
    end

    % Translate and scale metre coordinates before solving the LP, avoiding
    % loss of precision for small patches in large projected coordinates.
    origin = mean(tri,1);
    localTri = (tri-origin)/longest;
    if twiceArea < 0
        localTri = localTri([1 3 2],:);
    end
    next = [2 3 1];
    edges = localTri(next,:) - localTri;
    lengths = sqrt(sum(edges.^2,2));
    normals = [edges(:,2),-edges(:,1)]./lengths;
    b = sum(normals.*localTri,2);
    A = [normals,ones(3,1)];

    % LP: max r, subject to n_i*[x;y]+r <= b_i and r >= 0.
    % A positive r certifies an actual interior feasible point.
    [z,status] = solveInteriorLP(A,b,solver);
    violation = max([A*z-b;-z(3);0]);
    if ~strcmp(status,'OPTIMAL') || any(~isfinite(z)) || ...
            z(3) <= 0 || violation > 1e-7
        error('EmptyCandidateRepair:LPNotFeasible', ...
            'Triangle interior LP failed: status=%s, violation=%g.',status,violation);
    end
    xy = origin + longest*z(1:2).';
    if ~all(isfinite(xy)) || ~isinterior(regionPS,xy(1),xy(2))
        error('EmptyCandidateRepair:PointOutsideRegion', ...
            'Solver point failed validation against the true feasible region.');
    end
    candidatePoints(end+1,:) = xy; %#ok<AGROW>
    radiusM(end+1,1) = longest*z(3); %#ok<AGROW>
    refinedCount = refinedCount+1;
    maxObservedEdge = max(maxObservedEdge,longest);
    maxLPViolation = max(maxLPViolation,violation);
end

if isempty(candidatePoints)
    error('EmptyCandidateRepair:NoFeasibleCandidates', ...
        'No solver-certified interior candidates were generated.');
end
[candidatePoints,keep] = unique(candidatePoints,'rows','stable');
radiusM = radiusM(keep);
if ~all(isinterior(regionPS,candidatePoints(:,1),candidatePoints(:,2)))
    error('EmptyCandidateRepair:ValidationFailed', ...
        'The final candidate library contains a point outside the true region.');
end

info = struct();
info.trigger = 'original multi-angle hexagonal candidate library was empty';
info.method = 'true-polyshape triangulation + longest-edge subdivision + interior LP';
info.solver = solver.name;
info.lp_objective = 'maximize the minimum distance to all triangle edges';
info.source_triangle_count = size(conn,1);
info.refined_triangle_count = refinedCount;
info.degenerate_triangle_count = degenerateCount;
info.num_candidates = size(candidatePoints,1);
info.target_max_triangle_edge_m = targetEdge;
info.observed_max_triangle_edge_m = maxObservedEdge;
info.lp_normalized_max_constraint_violation = maxLPViolation;
info.candidate_triangle_interior_radius_m = radiusM;
info.all_candidates_in_true_region = true;
info.num_regions = regionPS.NumRegions;
info.num_holes = regionPS.NumHoles;
info.turbine_minimum_spacing_m = dMin;
info.layout_selection = 'unchanged random-key NSGA-II and wake evaluation';
end

function solver = chooseLPSolver()
solver = struct();
if exist('gurobi','file') == 3 || exist('gurobi','file') == 2
    solver.name = 'Gurobi MATLAB LP';
    solver.backend = 'gurobi';
    solver.options = struct('OutputFlag',0,'Threads',1, ...
        'FeasibilityTol',1e-9,'OptimalityTol',1e-9);
    return;
end
% Gurobi's examples folder may contain a linprog.m script. Accept only the
% Optimization Toolbox implementation, never a shadowing example script.
linprogPath = which('linprog');
optimRoot = fullfile(matlabroot,'toolbox','optim');
if license('test','Optimization_Toolbox') && ...
        strncmpi(linprogPath,optimRoot,numel(optimRoot))
    solver.name = 'MathWorks linprog LP';
    solver.backend = 'linprog';
    solver.options = optimoptions('linprog','Display','none', ...
        'ConstraintTolerance',1e-9,'OptimalityTolerance',1e-9);
    return;
end
error('EmptyCandidateRepair:NoLPSolver', ...
    'No Gurobi MATLAB or unshadowed Optimization Toolbox linprog is available.');
end

function [z,status] = solveInteriorLP(A,b,solver)
f = [0;0;-1];
lb = [-inf;-inf;0];
if strcmp(solver.backend,'gurobi')
    model = struct();
    model.A = sparse(A);
    model.obj = f;
    model.rhs = b;
    model.sense = repmat('<',3,1);
    model.lb = lb;
    model.modelsense = 'min';
    solution = gurobi(model,solver.options);
    status = solution.status;
    if isfield(solution,'x')
        z = solution.x;
    else
        z = nan(3,1);
    end
else
    [z,~,exitflag] = linprog(f,A,b,[],[],lb,[],solver.options);
    if exitflag == 1
        status = 'OPTIMAL';
    else
        status = sprintf('linprog exitflag %d',exitflag);
    end
    if isempty(z)
        z = nan(3,1);
    end
end
end

function value = cross2d(a,b)
value = a(1)*b(2)-a(2)*b(1);
end
