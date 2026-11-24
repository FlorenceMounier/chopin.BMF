
<!-- README.md is generated from README.Rmd. Edit README.Rmd, not README.md. -->

# chopin.BMF

<!-- badges: start -->

[![R-CMD-check](https://github.com/FlorenceMounier/chopin.BMF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/FlorenceMounier/chopin.BMF/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`chopin.BMF` is a reproducible R research package for assessing
**field-based predator-to-diet contaminant concentration ratios**
(operational biomagnification factors, BMFs) in juvenile common sole
(*Solea solea*) and benthic prey from the Seine estuary, France.

The analysis covers polychlorinated biphenyls (PCBs),
hexabromocyclododecane (HBCDD) isomers, and per- and polyfluoroalkyl
substances (PFAS). It combines juvenile sole chemistry, benthic prey
chemistry, stomach-content information, and literature-derived
individual prey masses to reconstruct diet-based concentration
denominators and investigate methodological sensitivity.

These BMFs are **field-based, concentration-normalised ratios**. They
should not be interpreted as direct measurements of contaminant uptake,
trophic transfer efficiency, or the physical contaminant concentration
of a reconstructed meal.

The repository contains the analyses, supplementary sensitivity
assessments, and source files for the main manuscript and its
Supplementary Information (SI).

## Repository contents

| Path | Purpose |
|:---|:---|
| `R/`, `data/`, `man/` | Package code, generated package datasets, and documentation (where present). |
| `data-raw/00_CHOPIN_config.R` | Central analysis settings, including sole sampling season. |
| `data-raw/00_RUN_ALL_FROM_SCRATCH.R` | Main entry point for rebuilding the analytical outputs. |
| `data-raw/00_CHOPIN_preflight.R` | Checks the availability of source data, required software, and input workbook. |
| `data-raw/01_CHOPIN_contamination_gamma.R` | Prepares biotic contaminant data, compound decisions, and censoring scenarios. |
| `data-raw/02_CHOPIN_diet.R` | Reconstructs diet proportions using stomach contents and reference prey masses. |
| `data-raw/03_CHOPIN_BMF_gamma.R` | Computes field-based BMF estimates and methodological comparisons. |
| `data-raw/04_CHOPIN_gamma_vs_lognormal.R` | Compares gamma and lognormal approaches for censored concentrations. |
| `data-raw/05_CHOPIN_sample_inventory_supplementary.R` | Produces sampling-inventory material for the SI. |
| `data-raw/06_CHOPIN_BMF_mass_quantiles.R` | Propagates Q25, median, and Q75 prey-mass assumptions to deterministic BMF sensitivity ranges. |
| `data-raw/07_CHOPIN_leave_one_prey_out.R` | Tests the influence of individual chemically analysed prey species. |
| `data-raw/08_CHOPIN_prey_season_sensitivity.R` | Examines the sensitivity to the seasonal composition of prey-chemistry samples. |
| `data-raw/09_CHOPIN_HBCDD_isomer_bounds.R` | Computes censoring-bound analyses for HBCDD isomers. |
| `data-raw/10_CHOPIN_postrun_provenance.R` | Checks the provenance and consistency of regenerated outputs. |
| `data-raw/11_CHOPIN_checked_submission_snapshot.R` | Creates a dated submission snapshot after verification (optional). |
| `manuscript/Article.Rmd`, `manuscript/Article_SI.Rmd` | Main article and SI source documents. |
| `manuscript/tables/`, `manuscript/figures/` | Generated tables and figures used by the manuscripts. |
| `results/` | Intermediate results, scenario comparisons, sensitivity analyses, and audit outputs. |
| `inst/extdata/Gogina_data-2021-0002-01.xlsx` | Literature-derived prey-mass source workbook required for the analysis. |

The source scripts in `data-raw/` are needed to reproduce the study;
they are **not** included in a standard built-package installation
because `data-raw/` is listed in `.Rbuildignore`.

## Requirements and source data

- R **\>= 4.1** and `dplyr` **\>= 1.1.0**.
- The source-data package
  [`chopin.data`](https://github.com/FlorenceMounier/chopin.data),
  containing `soles_contam`, `benthos_contam`, `soles_stomach`,
  `soles_G0_pools_biometrics`, `contam_info`, and `species_classif`.
- The reference workbook `inst/extdata/Gogina_data-2021-0002-01.xlsx`
  supplied with the source repository.
- R packages used by the preflight and analytical pipeline: `tidyverse`,
  `Btrait`, `readxl`, `janitor`, `survival`, `usethis`, `ggrepel`,
  `scales`, `tidytext`, and `lubridate`. Optional or output-specific
  dependencies include `flextable`, `officer`, `officedown`, `knitr`,
  `rmarkdown`, and `writexl`.

`chopin.data` supplies **source observations**. Rebuilding from scratch
means starting without **derived outputs** in `data/`, `results/`, and
the generated manuscript figure/table folders; it does **not** mean
rebuilding without the required source data.

The literature-derived individual prey-mass reference workbook
(`inst/extdata/Gogina_data-2021-0002-01.xlsx`) is version-controlled
through its SHA-256 checksum, recorded in
`data-raw/source_checksums.csv`. The preflight script verifies the
workbook’s integrity before running the analysis pipeline. This ensures
that reconstructed dietary biomass contributions are calculated using
the documented source file.

### Install dependencies

``` r
install.packages(c(
  "tidyverse", "readxl", "janitor", "survival", "usethis",
  "ggrepel", "scales", "tidytext", "lubridate",
  "flextable", "officer", "officedown", "knitr", "rmarkdown",
  "writexl", "pak"
))

# Install Btrait from the appropriate available source, if not installed.
# Install the source-data package (pin to an exact commit for publication):
pak::pak("FlorenceMounier/chopin.data@4c2c81ff26ddd5dce74873613a780f2c692baf99")
```

The unpinned installation above is intended for development only. **For
archival reproduction, replace it with the exact commit or release of
`chopin.data` used for the published analyses.** The specific version
has not yet been pinned in this README.

A standard installation of the research package itself may be useful for
inspecting package objects:

``` r
pak::pak("FlorenceMounier/chopin.BMF")
```

However, to rerun the full analysis, **clone or download the complete
source repository**, rather than relying on this installation alone.

## Full analytical rebuild

### 1. Obtain and open the source repository

Clone it with Git (or download the repository ZIP from GitHub), then
open `chopin.BMF.Rproj` in RStudio:

``` bash
git clone https://github.com/FlorenceMounier/chopin.BMF.git
cd chopin.BMF
```

Run all commands below from the **project root**: the directory
containing `DESCRIPTION` and `data-raw/`.

### 2. Check the analysis settings

The central configuration file is `data-raw/00_CHOPIN_config.R`. In the
archived pipeline, the main setting is:

``` r
SOLE_SEASON <- "autumn"
SAVE_TO_PACKAGE <- TRUE
```

`SOLE_SEASON = "autumn"` selects the age-0 sole chemistry used for the
reference analyses; `SAVE_TO_PACKAGE = TRUE` allows the scripts to
regenerate relevant objects in `data/`.

### 3. Run the pipeline

``` r
source("data-raw/00_RUN_ALL_FROM_SCRATCH.R")
```

This launcher executes the preflight and scripts **01–09 in order**,
each in a separate R environment. It then runs
`10_CHOPIN_postrun_provenance.R`. It does **not** create a submission
snapshot by default.

The main intermediate results include:

- `results/scenarios.rds`: censored-contaminant scenarios;
- `results/prey.rds`: chemically characterised prey data;
- `results/prop_diet.rds`: reconstructed dietary weights;
- `results/bmf_compare.rds`: reference and alternative field-based BMF
  estimates;
- `results/mass_quantiles/`: deterministic prey-mass sensitivity
  results;
- `results/reviewer_checks/`: provenance checks.

Scripts also create tables and figures under `manuscript/tables/` and
`manuscript/figures/`.

**Additional zero-compatible sensitivity analysis:** where included in a
revision of the repository,
`data-raw/03b_CHOPIN_zero_safe_sensitivity.R` must be run after script
03. It is **not listed in the current `00_RUN_ALL_FROM_SCRATCH.R`
launcher**; add it to that launcher or run it separately before
rendering SI sections that depend on its outputs.

``` r
# Only if this script is included in your working project:
source("data-raw/03b_CHOPIN_zero_safe_sensitivity.R")
```

### 4. Render the manuscript and Supplementary Information

The manuscript sources are `manuscript/Article.Rmd` and
`manuscript/Article_SI.Rmd`. They produce Word documents using the
project templates and bibliography.

Open each Rmd in RStudio and use **Knit**, or render programmatically
once paths to generated tables and figures are consistent with the
document’s working directory:

``` r
rmarkdown::render(
  "manuscript/Article.Rmd",
  knit_root_dir = normalizePath("manuscript")
)

rmarkdown::render(
  "manuscript/Article_SI.Rmd",
  knit_root_dir = normalizePath("manuscript")
)
```

The manuscript is **not** rendered by `00_RUN_ALL_FROM_SCRATCH.R`: run
the analytical scripts first, then render both documents. If a table or
figure cannot be found during knitting, inspect its path relative to
`manuscript/` and verify that the producing script was executed.

## Methodological scope and interpretation

The workflow reconstructs dietary chemical denominators using measured
concentrations of chemically characterised prey and dietary weights
estimated from stomach contents and literature-derived individual
masses. Chemical concentrations for PCB and HBCDD are lipid-normalised,
while PFAS follow their respective analytical normalisation. Diet
weighting does not turn these normalised concentrations into a physical
reconstructed meal.

Sensitivity assessments examine alternative treatment of censored
values, diet weighting, prey-mass Q25/median/Q75 assumptions, omission
of chemically characterised prey species, prey sampling season, and
HBCDD-specific censoring bounds.

- The **Q25/median/Q75 ranges are deterministic sensitivity ranges**,
  not confidence intervals or credible intervals; the minimum and
  maximum are calculated across the tested scenarios and need not
  correspond to Q25 and Q75 respectively.
- The **leave-one-prey-out** analysis tests the influence of species
  represented in the chemical dataset, **not** uncertainty about the
  identity of prey actually consumed.
- HBCDD censoring bounds are analytical bounds under specified
  assumptions, **not** statistical confidence intervals.
- Ratios above or below unity should not, by themselves, be interpreted
  as proof of trophic transfer efficiency, causal biomagnification, or
  in-vivo metabolism.

Refer to `manuscript/Article.Rmd` and `manuscript/Article_SI.Rmd` for
the detailed methods, selection rules, results, and limitations.

## Validation and provenance

The pipeline includes source and result checks. To assess package
consistency and record the computational environment:

``` r
# Run after the analysis from the package root
source("data-raw/10_CHOPIN_postrun_provenance.R")

# Package checks (additional development dependencies may be needed)
# install.packages("devtools")
devtools::check()

# Capture R / package versions
writeLines(capture.output(sessionInfo()),
           "results/sessionInfo_rebuild.txt")
```

For a publication-grade reproducibility record, document the **exact
`chopin.BMF` Git commit**, the **exact `chopin.data` commit/version**,
the workbook’s **SHA-256 fingerprint**, and the R/package versions. At
present, this README does not claim those revisions are frozen. A
successful provenance check should not be confused with an independently
verified clean-environment rebuild.

Once the analyses and manuscripts have been reviewed, the optional
snapshot can be created:

``` r
source("data-raw/11_CHOPIN_checked_submission_snapshot.R")
```

Alternatively, change the final call in
`data-raw/00_RUN_ALL_FROM_SCRATCH.R` to
`run_all(include_snapshot = TRUE)` for a validated submission run.
Snapshot creation is not enabled by default.

## Citation and availability

The manuscript and associated data publication are in preparation. Add
the **final article citation, data-paper citation, permanent archive
DOI, and exact software release/tag** once they are available. Do not
cite the repository’s current `main` branch as a fixed archival version.

**Repository:** <https://github.com/FlorenceMounier/chopin.BMF>\
**Source-data repository:**
<https://github.com/FlorenceMounier/chopin.data>

## License

MIT License; see `LICENSE`.
