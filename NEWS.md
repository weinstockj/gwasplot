# gwasplot (development version)

- No notable changes yet.

# gwasplot 0.2.0

## Added

- Added fixed-effects meta-analysis via `meta_analyze_fe()` for `GWASFormatter`, `data.frame`, and `tibble` inputs.
- Added allele harmonization for meta-analysis, including support for simple `REF`/`ALT` swaps by flipping the second study's `BETA`.
- Added focused meta-analysis regression tests, including `GWASFormatter` coverage and extreme-tail p-value checks.
- Added package-load and runtime checks for the DuckDB `stochastic` community extension used by DuckDB-backed meta-analysis.

## Changed

- Changed `GWASFormatter` materialization to use unique DuckDB-backed table names instead of a shared `summary_stats` table.
- Changed DuckDB-backed meta-analysis to require the `stochastic` extension and compute p-values with `dist_normal_cdf_complement()`.
- Changed the minimum `duckdb` dependency to `>= 1.3.2`.
- Changed `find_nearest_gene()` to use overlap-aware nearest-gene logic, treating variants inside genes as distance `0` and restricting candidates to protein-coding genes with non-null `gene_name`.

## Fixed

- Fixed cross-object table collisions when multiple `GWASFormatter` objects are created in the same working directory.
- Fixed DuckDB meta-analysis p-value behavior in the extreme tails by removing the old SQL approximation path and standardizing on `stochastic`.

## Documentation

- Expanded README coverage for meta-analysis workflows and standard GWAS output behavior.
