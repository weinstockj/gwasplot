# Repository Guidelines

## Project Structure & Module Organization
`gwasplot` is an R package organized with standard package conventions.
- `R/`: package source code (`manhattan.R`, `locuszoom.R`, `annotate.R`, etc.).
- `src/`: C++ helpers via Rcpp (`RcppExports.cpp`, `qqplot.cpp`).
- `man/`: generated Rd documentation (do not hand-edit).
- `data/`: packaged `.rda` reference datasets.
- `data-raw/`: scripts used to build data artifacts.
- `.github/workflows/pkgdown.yaml`: CI workflow for site/docs publishing.
- `renv/` and `renv.lock`: reproducible dependency environment.

## Build, Test, and Development Commands
Run from the repository root:
- `R -q -e "renv::restore()"`: install locked dependencies.
- `R -q -e "devtools::load_all()"`: load package for interactive development.
- `R -q -e "devtools::document()"`: regenerate `NAMESPACE` and `man/*.Rd` from roxygen comments.
- `R -q -e "devtools::test()"`: run the `testthat` suite.
- `R -q -e "devtools::check()"` or `R CMD check .`: run package checks.
- `R -q -e "pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)"`: build docs site (matches CI).
- `air format R/`: format R sources (config in `air.toml`).

See `CLAUDE.md` for the high-level architecture (the DuckDB-backed `GWASFormatter`
object, the S3 method pattern, and the in-place filtering pipeline).

## Coding Style & Naming Conventions
- Format R sources with `air` (`air.toml`: 80-col, 2-space indent). Keep
  formatting-only changes in their own commit, separate from behavior changes.
- Prefer `snake_case` for function and object names (for example, `reformat_summary_statistics`).
- Keep exported APIs documented with roxygen2 (`#'`) and regenerate docs with `devtools::document()`.
- Follow existing S3 method patterns (for `GWASFormatter`, `data.frame`, `tbl_df`) when adding generic behavior.

## Testing Guidelines
Tests live under `tests/testthat/` (testthat 3e), named `test-<feature>.R`. Before opening a PR:
- Run `devtools::test()` and `devtools::check()`; resolve all warnings/errors.
- Add `testthat` tests under `tests/testthat/` for non-trivial logic.
- Manually validate changed plotting/annotation paths with a small GWAS sample.

## Commit & Pull Request Guidelines
- Keep commit messages short, imperative, and specific (examples from history: `fix manhattan bug`, `improve API for qqplot function`).
- Separate refactors from behavior changes when possible.
- PRs should include:
  - a concise summary of user-visible changes,
  - linked issue(s) when relevant,
  - representative plots/screenshots for visualization changes,
  - confirmation that `devtools::document()` and `devtools::check()` were run.
