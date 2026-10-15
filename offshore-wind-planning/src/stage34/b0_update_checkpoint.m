function checkpoint = b0_update_checkpoint(checkpoint,out,epsilonIndices)
% Merge only the supplied epsilon results. Nonoptimal LP incumbents are diagnostics.
epsilonIndices = double(epsilonIndices(:).');
assert(numel(epsilonIndices)==numel(out.epsilon_set_MW));
assert(all(abs(checkpoint.epsilon_set_MW(epsilonIndices)- ...
    double(out.epsilon_set_MW(:).'))<=1e-9));
summary = b0_compact_candidates(out.by_epsilon);
statuses = out.storage_status_by_epsilon(:).';
for j=1:numel(epsilonIndices)
    % TIME_LIMIT incumbents do not reproduce the original optimal LP search.
    if ~(statuses(j).is_optimal && statuses(j).accepted_feasible)
        summary.feasible(1:2,j)=false;
        names=fieldnames(summary);
        for q=1:numel(names)
            name=names{q}; value=summary.(name);
            if isnumeric(value) && size(value,1)==4 && size(value,2)==numel(epsilonIndices)
                summary.(name)(1:2,j)=NaN;
            end
        end
        summary.storage_energy_MWh(1,j)=NaN;
    end
end
names = fieldnames(checkpoint.candidate_summary);
for q=1:numel(names)
    name = names{q};
    if strcmp(name,'row_labels'); continue; end
    assert(isfield(summary,name),'Candidate summary schema mismatch: %s',name);
    checkpoint.candidate_summary.(name)(:,epsilonIndices)=summary.(name);
end
if ~isfield(checkpoint,'status_history_by_epsilon')
    checkpoint.status_history_by_epsilon=cell(size(checkpoint.epsilon_set_MW));
end
for j=1:numel(epsilonIndices)
    k=epsilonIndices(j);
    previous=checkpoint.status_by_epsilon(k);
    if previous.attempted || previous.certified_pruned
        history=checkpoint.status_history_by_epsilon{k};
        if isempty(history); history=previous; else; history(end+1)=previous; end %#ok<AGROW>
        checkpoint.status_history_by_epsilon{k}=history;
    end
end
checkpoint.status_by_epsilon(epsilonIndices)=statuses;
if ~isfield(checkpoint,'attempt_count_by_epsilon')
    checkpoint.attempt_count_by_epsilon=zeros(size(checkpoint.epsilon_set_MW));
end
checkpoint.attempt_count_by_epsilon(epsilonIndices)= ...
    checkpoint.attempt_count_by_epsilon(epsilonIndices)+double([statuses.attempted]);

S = checkpoint.candidate_summary;
values = S.stage2_npv_CNY;
values(~S.feasible | ~isfinite(values)) = -Inf;
[maxNPV,linearIndex] = max(values(:));
assert(isfinite(maxNPV),'No finite feasible candidate after checkpoint update.');
[bestRow,bestColumn] = ind2sub(size(values),linearIndex);
cases = {'storage_optimized','storage_optimized','no_storage','no_storage'};
modes = {'HVAC','VSC-HVDC','HVAC','VSC-HVDC'};
epsilon = checkpoint.epsilon_set_MW(bestColumn);
matched = false;
for j=1:numel(out.candidates)
    candidate = out.candidates(j);
    if abs(candidate.epsilon_MW-epsilon)<=1e-9 && ...
            strcmp(candidate.storage_case,cases{bestRow}) && ...
            strcmp(candidate.export_mode,modes{bestRow}) && candidate.feasible
        assert(abs(candidate.stage2_npv_CNY-maxNPV)<=max(1,1e-10*abs(maxNPV)));
        checkpoint.selected_candidate = candidate;
        matched = true; break;
    end
end
if ~matched
    candidate = checkpoint.selected_candidate;
    assert(abs(candidate.epsilon_MW-epsilon)<=1e-9 && ...
        strcmp(candidate.storage_case,cases{bestRow}) && ...
        strcmp(candidate.export_mode,modes{bestRow}) && ...
        abs(candidate.stage2_npv_CNY-maxNPV)<=max(1,1e-10*abs(maxNPV)), ...
        'Checkpoint best-candidate dispatch is unavailable.');
end
checkpoint.search_complete = all(([checkpoint.status_by_epsilon.is_optimal]& ...
    [checkpoint.status_by_epsilon.accepted_feasible]) | ...
    [checkpoint.status_by_epsilon.certified_pruned]);
checkpoint.candidate_coverage_complete = all([checkpoint.status_by_epsilon.is_optimal]& ...
    [checkpoint.status_by_epsilon.accepted_feasible]);
checkpoint.updated_utc = char(datetime('now','TimeZone','UTC', ...
    'Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
end
