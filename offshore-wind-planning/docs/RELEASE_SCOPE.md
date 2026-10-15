# Scope and provenance of the code-only release

## Included

- Archived 7D Stage 1–2 micrositing, wake and platform-selection functions, including the validated empty-candidate-library fallback.
- Main production collection-system solver and the optional solver-incumbent initialization variant.
- Current B0 export/storage LP and dispatch normalization/physical-verification helpers, with effective storage purchase cost 1043.82 CNY/kWh.
- Section 4.4 regional selection/repair/revaluation model.
- New portable entry points, relative-path configuration, dependency/input checks and input schemas.

## Excluded

All observations, model inputs, computed results, checkpoints, archived baseline MAT files, weather MAT/NC/GRIB, GIS layers and sidecars, tabular result exports, figures, manuscripts, reviewer correspondence, inventories, local audit reports, Python runtimes and third-party solver/toolbox distributions are excluded.

The scope also excludes sensitivity/ablation experiments, 9D/11D drivers, high-storage-cost experiments, meteorological processing, plotting and historical one-grid/batch recovery scripts. The archived final-results directory is not part of this repository.

## Packaging changes

The scientific solver bodies are retained from the selected production sources. Public entry points replace machine-specific production orchestration:

- Full requested intervals replace hard-coded rerun batch slices.
- Data and output paths come from one repository-root configuration.
- Stage 3/4 consumes a newly generated Stage 1–2 master.
- Regional preparation consumes a newly generated B0 baseline and generates a local runtime cache.
- Stage interfaces retain global grid identifiers and the fields needed downstream.

Such packaging is not a new experiment or a claim of stronger optimality. Model limitations, finite budgets and the distinction between a clean rerun and historical merged results remain documented.

The release folder includes no `LICENSE` file because the authors have not selected a public reuse license. No repository URL, DOI or author account has been fabricated.
