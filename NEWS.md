# gwasplot (development version)

## Changed

- `manhattan()` gene labels are now lifted above their points and repelled along
  the y-axis only, so each name sits directly above its variant instead of
  overlapping the points. The lift is controlled by the new `label_nudge_y`
  argument (default: 4% of the plotted y-range; set to 0 to disable).
- `qqplot()` now renders with a fixed square panel (`aspect.ratio`, default 1)
  instead of `coord_equal()`. A strong GWAS with very large observed
  `-log10(p)` no longer stretches the plot into a tall, thin rectangle. Pass
  `aspect.ratio` to adjust. The y = x reference line is still drawn but no
  longer sits at a literal 45 degrees.

## Fixed

- Fixed `manhattan()` labeling variants that do not clear genome-wide
  significance. Labels are now gated to variants passing the significance
  threshold (default 5e-8), so `label_top_n` labels up to N *significant* hits
  rather than always labeling exactly N points down to `lower_logp_threshold`.

## Added

- Added `exclude_aberrant_pvalue_loci()` to drop suspicious "lonely peak"
  loci: a very strong lead variant (`PVALUE < lead_pvalue_threshold`, default
  `5e-10`) with no nearby supporting variant (none reaching
  `support_pvalue_threshold`, default `5e-5`, within `window_kb`, default
  `500`). The number of variants dropped is logged. Supported for
  `GWASFormatter`, `data.frame`, and `tibble` inputs.
- Added `keep = "union"` to `meta_analyze_fe()`. The default `"overlap"` keeps
  only variants present in both studies; `"union"` keeps variants present in
  either study, passing single-study variants through with that study's estimate
  and `N_studies = 1` (`matched_by` of `"x_only"` / `"y_only"`). Supported for
  `GWASFormatter` (DuckDB full join), `data.frame`, and `tibble` inputs.
- Added `write_summary_statistics()` to persist a (filtered) summary-statistic
  object to Parquet or CSV. For a `GWASFormatter` the write streams directly
  from DuckDB without collecting into R, capturing filters applied via
  `filter_variants()` / `exclude_difficult_regions()`; `data.frame` and
  `tibble` inputs are supported via a temporary in-memory DuckDB connection.
- Added `label_pvalue_threshold` to `manhattan()` to control the labeling
  cutoff independently of the plotted points; defaults to the genome-wide
  significance line. Genes listed in `highlight_genes` are always labeled
  regardless of this threshold.

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
