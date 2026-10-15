# Stage I inputs and outputs

`run_stage1` runs the paper's baseline turbine micrositing and platform
selection modules. It uses the 2023 hourly input convention, a 12 MW turbine,
a 242 m rotor, a minimum spacing of 7D, and availability of 0.95.
GIS screening, meteorological preprocessing, research inputs, and computed
results are not included in this repository.

## Entry point

From the repository root:

```matlab
cfg = release_config();
report = run_stage1(cfg);
```

Configure the input locations through `release_config(overrides)`:

```matlab
o = struct();
o.dataDir = '/path/to/your/inputs';
o.outputDir = '/path/to/new/outputs';
o.gisFile = '/path/to/your/feasible_final.shp';
o.speed140File = '/path/to/your/speed_140m_2023.mat';
o.direction140File = '/path/to/your/direction_140m_2023.mat';
cfg = release_config(o);
report = run_stage1(cfg);
```

The default requests all shapefile records. For a subset, set `startIndex`
and `endIndex` in the overrides. These are one-based **shapefile record
positions**, not an `OBJECTID`, geographic grid identifier, or chunk number.
Every chunk in the requested range is processed; no machine-specific slice
is embedded in this entry point.

The eight scientific functions under `src/stage1/functions` are copied
without numerical or algorithm changes from the archived baseline package.
The released driver changes paths, the function entry point, and batch
orchestration. In particular, the saved empty-candidate-library fallback
remains available. The current source defaults to serial micrositing;
`useParallelMicrositing` and `numWorkers` control optional parallel execution.

## Feasible-domain shapefile

`cfg.gisFile` identifies a complete, user-supplied shapefile dataset. Retain
its required `.shp`, `.shx`, and `.dbf` components and its CRS description
(`.prj`). Record ordering must agree with all compact-wind columns.

The model accesses these fields:

| Field | Meaning and required units |
| --- | --- |
| `X`, `Y` | Feasible polygon boundary coordinates in projected metres. Preserve NaN-separated rings, internal holes, and disconnected components. The original case uses WGS 84 / UTM zone 50N, EPSG:32650. |
| `X1`, `Y1` | Grid wind-data lookup coordinates: longitude and latitude in degrees, respectively. These are distinct from projected polygon coordinates. |
| `dist2Port` | Distance to port in kilometres, used by the platform lifecycle-cost model. |
| `RASTERVALU` | Signed seabed elevation in metres: negative values below sea level. The driver computes depth as `abs(min(RASTERVALU,0))`; a positive depth magnitude must not be placed directly in this field. |

`X1` and `Y1` are used when the compact wind database must be built.
Stage II additionally requires the landing-point/site fields documented in
its own input contract; retain those fields when preparing a shared GIS
dataset. Stage I does not reproject the shapefile or infer a landing point.

## Compact wind database: preferred prepared input

If `cfg.compactWindFile` exists and `rebuildCompactWind` is false, the raw
speed/direction cubes are not opened. The MAT file must contain a scalar
structure named **`wind2023`** with these fields:

| Field | Shape or value |
| --- | --- |
| `wind_speed_140m_mps` | `T x N` numeric array of hub-height wind speed in m/s. |
| `wind_direction_from_deg` | `T x N` numeric array of meteorological FROM direction in degrees: north=0, east=90, south=180, west=270. |
| `power_pu` | `T x N` available no-wake single-turbine power divided by 12 MW, using the included power curve and availability 0.95. Availability must not be applied twice. |
| `grid_count` | `N`, equal to the number of shapefile records. |
| `num_hours` | `T`; the paper's non-leap 2023 input uses 8760 hourly observations. |
| `rated_power_MW` | `12`. |
| `availability_factor` | `0.95`. |

Column `i` must describe shapefile record `i`, for every array. Array sizes
alone cannot establish geographic alignment. Ensure the same hourly order
and coordinate mapping before running the model.

The builder also records `grid_index`, `grid_lon`, `grid_lat`, provenance,
and optional `source_speed_row_col` / `source_direction_row_col` fields.
`source_speed_row_col`, if present, is copied to each result's
`wind_data_index`. Without it the driver uses the record position.

The driver recomputes no-wake power inside micrositing and checks it against
`power_pu` to numerical tolerance. Do not supply an aerodynamic curve without
the specified availability factor or a wake-adjusted series as `power_pu`.

By default the compact cache is generated below the configured output
directory. To use a separately prepared compact MAT, set
`compactWindFile` explicitly. To rebuild it, set `rebuildCompactWind=true`.

## Optional already-processed wind cubes

Only if the compact file is missing, or rebuilding is requested,
`build_2023_compact_wind_db` reads `speed140File` and `direction140File` and
saves `wind2023` at `compactWindFile`. These inputs are already processed
140 m wind fields. The repository does not download ERA5 data, perform
ERA5/ERA5-Land blending, extrapolate wind speed to hub height, or generate
the GIS exclusions.

Each data array is spatial-dimension 1 by spatial-dimension 2 by hour.
Supported variable names, in preference order, are:

- Speed: `wind_speed_140m_southchina_sea`, `wind_speed_140m`,
  `wind_speed_southchina_sea`, `wind_speed`.
- Direction: `wind_direction_140m_southchina_sea`, `wind_direction_140m`,
  `wind_direction_southchina_sea`, `wind_direction`.

The existing loader can fall back to the largest numeric three-dimensional
variable when a preferred name is absent. Use a preferred name to avoid
ambiguity. Latitude/longitude can be coordinate grids or compatible vectors
named `lat`/`lon`, `latitude`/`longitude`, `LAT`/`LON`,
`Latitude`/`Longitude`, `grid_lat`/`grid_lon`, or `lat_grid`/`lon_grid`.
If direction coordinates are omitted, the builder uses the speed-file
coordinates. The source array dimensions must still match those coordinates.

The retained source convention first treats the first spatial axis as
longitude and the second as latitude when compatible. It also accepts a
latitude-first vector grid when its axis lengths establish that ordering.
For equal-length spatial axes, supply unambiguous coordinate grids that
follow the source convention rather than relying on vector lengths.
The default coordinate lookup tolerance is 0.02 degrees; the builder matches
grid locations without spatial interpolation.

## Checkpoints and handoff to Stage II

The master file is written to `cfg.stage12File`. Chunk files are written to
a sibling `chunks` directory. The master contains:

- `stage12_results`: one record for every shapefile position. A subset run
  leaves unrequested positions as empty records.
- `stage12_grid_indices`: the one-based positions requested for this run.
- `stage12_meta`: requested range, record counts, successful/failed counts,
  timestamps, and input/output provenance.
- `cfg_stage12`: the runtime configuration.

Successful records have `status='ok'` and retain the selected turbine
coordinates, capacity, turbine count, hourly wake-adjusted power,
availability, annual degradation, platform choice and cost breakdown,
Pareto candidate summaries, and diagnostic timings. In particular,
`wake_adjusted_power_pu_2023` is retained for the Stage II hourly calculations.

Chunk metadata `completed=true` means the chunk loop finished. It does not
mean every grid solved successfully: inspect `stage12_meta.num_error` and
per-grid `status` / `error_message`. `report` returns lightweight master
metadata and `masterFile`; it does not return the hourly database in memory.

`resume=true` retains the original checkpoint behavior: already completed
chunks are skipped. It does not validate a scientific input signature and
does not automatically retry individual error records inside a completed
chunk. Use a new output directory when inputs, scientific parameters, or
the grid ordering change; inspect error records before proceeding to Stage II.

## MATLAB dependencies

- The Stage I shapefile reader requires Mapping Toolbox.
- The RK-NSGA-II implementation is self-contained and does not call
  `gamultiobj` or require Global Optimization Toolbox.
- Parallel Computing Toolbox is optional. The source implementation falls
  back to serial evaluation when parallel execution is unavailable.
- Empty original candidate libraries use the retained LP fallback. This
  requires either a configured Gurobi MATLAB interface or the unshadowed
  MathWorks `linprog` from Optimization Toolbox. No solver distribution is
  bundled in this repository.

The archived implementation was checked with MATLAB R2022a. Publishing this
source package does not imply that the full regional optimization has been
rerun in a new MATLAB or solver environment.
