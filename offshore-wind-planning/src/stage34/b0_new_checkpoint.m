function checkpoint = b0_new_checkpoint(D,column,cfg,noStorageOut)
% New checkpoint. Every explicit no-storage branch is evaluated analytically.
assert(all(~[noStorageOut.storage_status_by_epsilon.attempted]), ...
    'Initialization must not spend an LP solve.');
checkpoint = struct();
checkpoint.schema_version = 'B0_STAGE4_PLATFORM_ZERO_V1';
checkpoint.model_signature = cfg.ModelSignature;
checkpoint.dataset_id = char(D.dataset_id);
checkpoint.grid_index = double(D.grid_ids(column));
checkpoint.input_column_index = double(column);
checkpoint.capacity_MW = double(D.capacity_MW(column));
checkpoint.distance_to_shore_km = double(D.distance_to_shore_km(column));
checkpoint.platform_turbine_cost_CNY = double(D.platform_turbine_cost_CNY(column));
checkpoint.internal_cable_cost_CNY = double(D.internal_cable_cost_CNY(column));
checkpoint.module3_gap = double(D.module3_gap(column));
checkpoint.storage_parameters = cfg.StorageParameters;
checkpoint.epsilon_set_MW = double(noStorageOut.epsilon_set_MW(:).');
checkpoint.status_by_epsilon = noStorageOut.storage_status_by_epsilon(:).';
checkpoint.candidate_summary = b0_compact_candidates(noStorageOut.by_epsilon);
checkpoint.selected_candidate = noStorageOut.best;
checkpoint.pass_count = 0;
checkpoint.wall_time_seconds = 0;
checkpoint.attempt_count_by_epsilon = zeros(size(checkpoint.epsilon_set_MW));
checkpoint.status_history_by_epsilon = cell(size(checkpoint.epsilon_set_MW));
checkpoint.search_complete = false;
checkpoint.candidate_coverage_complete = false;
checkpoint.seed_provenance = struct();
checkpoint.last_error_message = '';
checkpoint.updated_utc = utcTimestamp();
end

function value = utcTimestamp()
value = char(datetime('now','TimeZone','UTC', ...
    'Format',"yyyy-MM-dd'T'HH:mm:ss'Z'"));
end
