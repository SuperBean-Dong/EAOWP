# Private inputs — files are not distributed

Supply these files locally, or override their paths with `release_config(struct(...))`:

```text
data/
├── spatial/
│   ├── feasible_final.shp
│   ├── feasible_final.shx
│   ├── feasible_final.dbf
│   └── feasible_final.prj
├── weather/
│   ├── south_sea_wind_speed_140m_2023.mat
│   └── south_sea_wind_direction_140m_2023.mat
└── economics/
    └── electricity_price.mat
```

These are expected filenames, **not included files**. No dataset is silently downloaded. Meteorological preprocessing and marine-exclusion GIS construction are outside this main-model release.

## Spatial input

The shapefile contains the inherited feasible farm polygons after marine exclusions. `X`/`Y` boundaries are in projected metre coordinates, with nonconvex boundaries and interior rings retained. Attributes used by the model include geographic reference coordinates `X1`/`Y1`, `dist2Port` in km, negative bathymetric elevation `RASTERVALU` in m, and the landing/shore information required by Stage 3/4. The exact fields and alternatives are documented in `docs/STAGE1_INPUTS.md` and `docs/STAGE34.md`.

Do not reorder grids independently of weather columns or saved records. The manuscript used 5287 feasible grids; code availability does not provide access to that private spatial dataset.

## Meteorological input

The planning input is the 2023 ERA5-Land/ERA5 hybrid wind-speed and meteorological wind-direction series extrapolated to 140 m. See the Stage 1 input documentation for the required MAT structures and dimensions. An existing compatible compact wind file may be supplied through `compactWindFile`; otherwise Stage 1 constructs a compact file under `outputs/cache/` from the two weather MAT inputs.

The compact wind file is runtime data, not source code, and is not distributed.

## Electricity prices

Supply at least 8760 finite hourly prices in CNY/MWh, in the MAT format required by `docs/STAGE34.md`. Stage 3/4 uses the first 8760 values, retaining the original production policy; do not substitute different units or a different hour alignment. The regional model additionally requires those prices to be strictly positive, as required by its original core. The Stage 3/4 documentation specifies the accepted variable names.

## Intermediate inputs

You do **not** need pre-existing research results to start the pipeline:

1. Stage 1 writes a new `stage12_results` master file.
2. Stage 3/4 reads that master and writes new baseline `results`.
3. The regional adapter derives its own input cache from that baseline and the configured price file.

These are generated artifacts in `outputDir`, not additional datasets bundled with the repository. Missing inputs should stop execution with an actionable error; the program must not replace missing scientific data with synthetic defaults.
