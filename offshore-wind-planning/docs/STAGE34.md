# Farm electrical design and delivered-shore baseline

Run `run_stage1(cfg)` first, then `run_stage34(cfg)`. The second command computes
new Stage-3 collection designs from the saved Stage-I layouts and then computes
the current Stage-4 baseline at an effective storage purchase cost of
**1043.82 CNY/kWh**. It never loads a previously archived B0 or a historical
economic scenario. No results, weather files, GIS files, or other research
inputs are distributed with this source package.

```matlab
cfg = release_config();
run_stage1(cfg);
report = run_stage34(cfg);
```

For a partial run, set the same original GIS grid indices in the configuration:

```matlab
cfg = release_config(struct('gridIDs', [1 2 3]));
run_stage34(cfg);
```

Only completed Stage-I records can be passed to Stage 3. `startIndex` and
`endIndex` also filter the original GIS indices. `gridIDs=[]` means all indices
available in the Stage-I master within that range; there is no fixed grid list,
five-computer assignment, or special rerun subset.

## External requirements

Provide MATLAB, Mapping Toolbox for `shaperead`, YALMIP, and a separately
installed and licensed Gurobi MATLAB interface. These dependencies are not
vendored. The unchanged storage implementation calls Gurobi directly; collection
design uses YALMIP and Gurobi. Normal execution does not require Python.

The SHA-256 helper uses MATLAB's JVM when available. Its archived no-JVM fallback
is Windows PowerShell; a non-Windows MATLAB session should retain the JVM.

## Inputs

All paths come from `release_config`. Put your inputs outside the tracked source
tree, or under the ignored local data directory.

| Configuration | Required contents |
| --- | --- |
| `gisFile` | Original ordered feasible-grid SHP and its required sidecars. `X`, `Y` contain projected boundary rings in metres, including NaN-separated holes/components. **`X_UTM`, `Y_UTM` are the onshore landing-point coordinates**, not the farm centroid. Use a common projected coordinate system for boundaries, turbines, and landing points. |
| `stage12File` | MAT produced by `run_stage1`, containing `stage12_results`, `stage12_grid_indices`, and `stage12_meta`. Original one-based SHP indices are retained. |
| `priceFile` | MAT containing at least 8760 first-year hourly selling prices in CNY/MWh. The first 8760 values are used and must be finite; zero and negative prices are not rejected by the wrapper. Recognised variable names are `price_CNY_per_MWh`, `Price_yuan_per_MWh`, `power_price`, and `price`. An otherwise unambiguous numeric variable with at least 8760 values is also accepted. |

Each requested Stage-I record must have `status='ok'`, `data_year=2023`,
`grid_index`, `wind_data_index`, `wind_farm_capacity`, `turbine_count`,
`turbine_xy_m`, `wake_adjusted_power_pu_2023`, `availability_factor`,
`annual_wind_degradation_rate`, `optimal_platform`, `platform_turbine_cost`, and
`platform_detail`. The wind profile has 8760 entries and already incorporates
the 95% availability factor and layout-dependent wake losses. It must not be
multiplied by availability a second time. `platform_turbine_cost` is a lifecycle
present cost in CNY, rather than just initial support CAPEX.

The baseline uses 12-MW turbines, the archived 7D Stage-I layouts, 25 project
years, 3% discounting, 0.72% annual generation-performance decline, and
8.39 CNY/EUR. The collection cables retain capacities `[40,50,63]` MW and
prices `[423,474,554]` EUR/m. The shore-distance objective proxy remains
815 EUR/m. The exact scientific parameter structure is stored in
`baseline_meta.scientific_parameters`.

## Collection solver and finite budgets

`cfg.collectionSolver='fast'` is the original main production strategy. It
constructs a physically feasible star incumbent, retains the complete joint
collection model, and uses safe block bounds/pruning. Explicitly select
`'solver_init'` to use the archived supplementary strategy that accepts only
solver-generated incumbents and raises an error if none is found. Both
scientific function files are included unchanged; the driver does not silently
switch strategies.

The historical final B0 combined an earlier main run and supplementary runs.
A fresh portable run follows the strategy you choose for every requested grid.
Its finite-budget solutions can differ from the historical merged result and
are not promised to reproduce it byte for byte.

| Runtime field | Default |
| --- | --- |
| `collectionSolver` | `'fast'` |
| `ossBlockTimeLimitSeconds` | 60 |
| `ossGridTimeLimitSeconds` | 600 |
| `ossMIPGap` | 0.05 |
| `ossEnablePruning` | `true` |
| `solverThreads` | 2 |
| `storageTimeLimitSeconds` | 300 per LP |
| `storageGridTimeLimitSeconds` | 1800 per pass |
| `storageMaxPasses` | 10 |
| `saveStage34Every` | 1, controls the small summary CSV interval |

The storage candidate set is 50:10:installed-capacity MW, supplemented by the
exact actual wind-power peak and exact installed capacity, with nonpositive and
duplicate entries removed. When installed capacity is below 50 MW, the two
exact endpoints remain. Both storage-optimised and explicit no-storage branches
are evaluated for each export candidate. Their HVAC/VSC-HVDC cable alternatives
are compared using delivered-shore NPV. The original wind-side LP followed by
export/loss valuation is preserved; this does not create a new feedback
optimisation of dispatch against every electrical design.

The driver initialises analytical no-storage branches, solves only missing LP
candidates, saves after every attempted candidate, and retries missing
candidates for up to `storageMaxPasses` per invocation. It uses no imported cache
seeds, pruning certificates, or old financial results. Incomplete LP searches
are explicitly marked `pending_stage4`, not accepted as completed baselines.

## Outputs and resume

`cfg.baselineFile` defaults to
`outputs/b0/baseline_results_2023_7D_B0.mat`. Its variables are:

- `results`: a structure array indexed by the original SHP grid index;
- `baseline_meta`: scientific settings, selected solver, input SHA-256 values,
  requested and completed grid counts, and `all_complete`;
- `price_CNY_per_MWh`: the exact 8760-value price vector used for the run;
- `cfg_stage34`: the run configuration.

Each grid also has `outputs/b0/checkpoints/grid_00001.mat`, etc., containing
three variables: `record`, `checkpoint`, and `run_signature`. `record` is the
same grid record used in the master. `checkpoint` retains per-epsilon storage
statuses, numerical candidate tables, attempt history, the selected dispatch,
and `candidate_coverage_complete`. Its legacy internal schema identifier is
retained to keep the archived checkpoint helper functions unchanged.

`baseline_results_summary.csv` is a small generated index. The large master is
installed at the end of the invocation (or before a configured stop-on-error),
while candidate checkpoints are saved throughout. This avoids rewriting a
multi-gigabyte master after every grid. If a run is interrupted before the master
is installed, rerun the same configuration to reconstruct it from checkpoints.

Resume checks the content hashes of Stage-I, prices, SHP and available sidecars,
the copied source code, scientific parameters, collection solver, and key
budgets. A changed identity raises an error rather than reusing another run's
checkpoint. Use a different `baselineFile` directory for another strategy or
budget. `resume=false` intentionally starts fresh in the chosen output location.

`status='ok'` means a feasible collection incumbent exists and every required
storage LP candidate has completed. It does **not** prove that every collection
solve met the target: consult `internal_cable_gap`, `module3_optimality_status`,
`module3_diagnostics`, and `optimization_complete`. Metadata `all_complete`
requires an `ok` record for every GIS grid, so a partial run is not mistaken for
a complete regional baseline.

## Result fields used by regional planning

The records retain the Stage-I geometry, original `wind_data_index`, AEP,
platform details, and wake-adjusted profile. Collection outputs include `beta`,
`internal_cable_cost`, `internal_collection_objective_cost`,
`internal_cable_gap`, `collection_edge_detail`, and `distance_to_shore_km`.
The shore-distance proxy is not charged again as a physical collection cost.

Economic outputs include `NPV_best`, `total_cost_best`,
`discounted_onshore_energy_kWh`, and `LCOE_min`. The legacy `LCOE_min` field is
the delivered-shore LCOE of the NPV-selected configuration; no separate
LCOE-minimisation run is implied.

`record.selected_candidate` is the unchanged core's complete selected candidate:

- `epsilon_MW`, `storage_energy_MWh`, `storage_case`, `export_mode`;
- `transmission_detail` (selected cable and circuit configuration);
- `farm_delivery_MW`, `onshore_power_kW`;
- `storage_detail` with `charge_MW`, `discharge_MW`, `energy_MWh`,
  `curtailment_MW`, and physical/solver status;
- `annual_delivered_energy_kWh`, `annual_revenue_CNY`, `pv_revenue_CNY`;
- `export_system_cost_CNY` and `battery_cost_CNY`.

The regional adapter can read these per-grid checkpoints and the saved price
vector without loading historical CSVs, historical B0 results, or all hourly
master records at once.
