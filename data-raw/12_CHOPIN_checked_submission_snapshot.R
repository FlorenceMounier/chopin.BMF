
# =============================================================================
# CHOPIN | 12 - Create a validated submission snapshot
# =============================================================================
#
# PURPOSE
# -------
# Create a timestamped, non-overwritable snapshot of the CHOPIN analyses,
# manuscript sources, supporting data and reproducibility information.
#
# The snapshot is intended to preserve the files used for a particular
# manuscript submission or reviewer response.
#
# IMPORTANT
# ---------
# Run this script ONLY AFTER:
#   1. Executing the complete CHOPIN production pipeline;
#   2. Running all reviewer-oriented sensitivity analyses;
#   3. Checking that the final results and tables are up to date;
#   4. Finalising the manuscript and Supplementary Information.
#
# This script does NOT recompute the BMF analyses.
# It checks key results, then freezes selected existing files.
#
# OUTPUT
# ------
# An external directory:
#
#   ../CHOPIN_submission_snapshots/
#       CHOPIN_validated_YYYYMMDD_HHMMSS/
#
# containing selected project files, a checksum manifest and sessionInfo.
#
# The Git history (.git/) and RStudio settings (.Rproj.user/) are excluded
# because the snapshot copies only explicitly selected files.
#
# REQUIREMENT
# -----------
# Run from the root directory of the chopin.BMF project.
# =============================================================================


# =============================================================================
# STEP 1 - Load packages and run provenance checks
# =============================================================================

# Load tidyverse for data manipulation and CSV export.
suppressPackageStartupMessages(
  library(tidyverse)
)

# Re-run the post-production provenance checks.
#
# Script 10 checks consistency between regenerated analytical results
# and the original source data, notably the biological pools.
#
# If these checks fail, execution should stop before any snapshot is created.
source(
  "data-raw/10_CHOPIN_postrun_provenance.R",
  local = TRUE
)


# =============================================================================
# STEP 2 - Check the presence of required analytical results
# =============================================================================

# Three RDS files are required to validate the reference BMF calculations:
#
# 1. bmf_compare.rds:
#    BMF estimates across dietary and censoring scenarios.
#
# 2. scenarios.rds:
#    Stored analytical scenarios.
#
# 3. reference_check.rds:
#    Independent comparison between gamma-based calculations
#    and the reference estimates produced by script 03.

RDS <- c(
  "results/bmf_compare.rds",
  "results/scenarios.rds",
  "results/gamma_vs_lognormal/reference_check.rds"
)

# Identify any missing files.
missing_rds <- RDS[!file.exists(RDS)]

# Stop if at least one required file is absent.
if (length(missing_rds) > 0L) {
  stop(
    "Missing RDS data: ",
    paste(missing_rds, collapse = ", ")
  )
}


# =============================================================================
# STEP 3 - Verify the reproducibility of gamma-based BMF estimates
# =============================================================================

# Load the full comparison of BMF scenarios.
bmf_results <- readRDS(RDS[1])

# Load the independent gamma reference check.
reference_check <- readRDS(RDS[3])

# Extract the maximum-likelihood censoring scenario.
bmf_mle <- dplyr::filter(
  bmf_results,
  scenario == "mle"
)

# Verify that:
#   - the independent check has the expected number of rows;
#   - every absolute difference is finite;
#   - differences remain below the specified numerical tolerance.
#
# The tolerance is relative to the BMF magnitude, with a minimum
# reference scale of 1:
#
#   tolerance = 1e-8 * max(1, reference BMF)
#
# This allows for negligible floating-point differences without
# accepting material discrepancies between implementations.

if (
  nrow(reference_check) != nrow(bmf_mle) ||
  any(!is.finite(reference_check$abs_diff)) ||
  any(
    reference_check$abs_diff >
    1e-8 * pmax(1, reference_check$bmf_main)
  )
) {
  stop(
    "Gamma BMFs do not reproduce script 03 results"
  )
}


# =============================================================================
# STEP 4 - Check that all expected analytical scenarios are present
# =============================================================================

# The BMF calculations must contain four treatments of censored values:
#
#   zero     : censored concentration replaced by zero
#   loq_half : censored concentration replaced by half the LOQ
#   loq      : censored concentration replaced by the LOQ
#   mle      : model-based gamma maximum-likelihood treatment

expected_scenarios <- c(
  "zero",
  "loq_half",
  "loq",
  "mle"
)

# Ensure that the set of scenarios matches the expected set.
if (
  !setequal(
    unique(bmf_results$scenario),
    expected_scenarios
  )
) {
  stop("Unexpected censoring scenarios")
}


# =============================================================================
# STEP 5 - Check the expected number of BMF estimates
# =============================================================================

# Count distinct contaminant compounds, including their family and unit.
n_compounds <- bmf_results |>
  dplyr::distinct(family, compound_name, unit) |>
  nrow()

# Each retained compound should have:
#
#   4 censoring treatments
# x 3 dietary weighting methods
# = 12 BMF estimates
#
# This check detects missing or additional combinations at the global
# dataset level.

expected_rows <- 4L * 3L * n_compounds

if (nrow(bmf_results) != expected_rows) {
  stop(
    "Unexpected number of BMF rows: ",
    nrow(bmf_results),
    " instead of ",
    expected_rows
  )
}

message(
  "Analytical checks passed for ",
  n_compounds,
  " compounds."
)


# =============================================================================
# STEP 6 - Define the snapshot destination
# =============================================================================

# Use date and time to identify the snapshot uniquely.
#
# Example:
#   20261126_103025
#
# This avoids collisions when several snapshots are created
# during the same day.

TAG <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

# Locate the parent directory of the project.
#
# For example, if the project is located at:
#
#   /home/user/projects/chopin.BMF/
#
# its parent directory is:
#
#   /home/user/projects/

PROJECT_PARENT <- dirname(
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
)

# Store snapshots outside the package directory.
#
# This prevents validated archives from being accidentally included
# in the package build or committed to the Git repository.

SNAPSHOT_ROOT <- file.path(
  PROJECT_PARENT,
  "CHOPIN_submission_snapshots"
)

# Create the snapshot root directory if it does not exist.
dir.create(
  SNAPSHOT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

# Construct the destination for the current snapshot.
OUT <- file.path(
  SNAPSHOT_ROOT,
  paste0("CHOPIN_validated_", TAG)
)

message("Snapshot destination: ", OUT)

# Never overwrite an existing validated snapshot.
#
# This safeguards previously frozen versions.
if (dir.exists(OUT)) {
  stop("Snapshot already exists: ", OUT)
}


# =============================================================================
# STEP 7 - Select the files to preserve
# =============================================================================

# First, list essential project files explicitly.
#
# These include:
#   - package metadata and documentation;
#   - manuscript and Supplementary sources;
#   - bibliography and Word templates;
#   - the original Gogina reference dataset.
#
# Paths are relative to the project root.

essential_paths <- c(
  "DESCRIPTION",
  "README.md",
  "manuscript/Article.Rmd",
  "manuscript/Article_SI.Rmd",
  "manuscript/references.bib",
  "manuscript/template.docx",
  "manuscript/template_SI.docx",
  "inst/extdata/Gogina_data-2021-0002-01.xlsx"
)

# Verify that every essential file is available.
#
# In the original version, missing files were silently excluded.
# Here, missing mandatory files cause an explicit error instead.

missing_essential <- essential_paths[
  !file.exists(essential_paths)
]

if (length(missing_essential) > 0L) {
  stop(
    "Missing essential snapshot files: ",
    paste(missing_essential, collapse = ", ")
  )
}

# Extend the inventory with generated files:
#
#   data-raw/           : R analysis and reproducibility scripts
#   results/            : CSV and RDS outputs
#   manuscript/tables/  : supplementary and manuscript tables
#   manuscript/figures/ : generated figures
#   data/               : package datasets
#
# list.files() selects individual files rather than copying directories.

paths <- unique(c(
  essential_paths,

  list.files(
    "data-raw",
    full.names = TRUE,
    pattern = "\\.R$"
  ),

  list.files(
    "results",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.(csv|rds)$"
  ),

  list.files(
    "manuscript/tables",
    full.names = TRUE
  ),

  list.files(
    "manuscript/figures",
    full.names = TRUE
  ),

  list.files(
    "data",
    full.names = TRUE,
    pattern = "\\.rda$"
  )
))

# Keep only paths pointing to existing files.
#
# In particular, folders are not copied recursively.
# Therefore, .git/ and .Rproj.user/ are not included.

paths <- paths[
  file.exists(paths) &
    !dir.exists(paths)
]

message(
  "Files selected for snapshot: ",
  length(paths)
)


# =============================================================================
# STEP 8 - Copy files while preserving the project directory structure
# =============================================================================

# Create the snapshot destination directory.
dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# Copy each file into the snapshot.
#
# Example:
#
#   results/bmf_compare.rds
#
# becomes:
#
#   CHOPIN_validated_YYYYMMDD_HHMMSS/results/bmf_compare.rds
#
# The directory structure is preserved so the archived project remains
# easy to inspect and compare with the original repository.

for (f in paths) {

  # Define the file destination.
  target <- file.path(OUT, f)

  # Create any missing parent directories.
  dir.create(
    dirname(target),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Copy without overwriting an existing file.
  success <- file.copy(
    from = f,
    to = target,
    overwrite = FALSE
  )

  # Stop if a copy operation fails.
  if (!success) {
    stop("Cannot freeze file: ", f)
  }
}


# =============================================================================
# STEP 9 - Record file checksums and data-package version
# =============================================================================

# Calculate an MD5 checksum for each copied file.
#
# A checksum is a fingerprint of file contents.
#
# If an archived file is altered after snapshot creation,
# its recalculated MD5 will normally differ from the stored value.
#
# MD5 is used here for file-integrity comparison, not for security
# or authentication purposes.

md5 <- tools::md5sum(
  file.path(OUT, paths)
)

# Create the manifest.
#
# Each line records:
#   - relative file path;
#   - MD5 fingerprint;
#   - installed chopin.data package version.
#
# Note: a package version is not the same as a Git commit SHA.
# The detailed provenance information is generated by script 11.

manifest <- tibble(
  file = paths,
  MD5 = unname(md5),
  source_chopin_data_version = as.character(
    utils::packageVersion("chopin.data")
  )
)

# Export the checksum manifest inside the snapshot.
readr::write_csv(
  manifest,
  file.path(OUT, "manifest_MD5.csv")
)


# =============================================================================
# STEP 10 - Record the R session information
# =============================================================================

# Save:
#   - R version;
#   - operating system;
#   - attached packages;
#   - package versions;
#   - session-specific configuration.
#
# This helps reviewers understand the computational environment.

writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "sessionInfo.txt")
)


# =============================================================================
# STEP 11 - Confirm successful snapshot creation
# =============================================================================

# The snapshot now contains:
#
#   - validated analytical outputs;
#   - selected R scripts;
#   - manuscript source files;
#   - supporting data;
#   - file checksums;
#   - R session information.
#
# It is a frozen copy of selected files, not a complete repository
# clone or a security-protected immutable storage system.
#
# A ZIP can be created separately if needed for reviewer distribution.

message(
  "Validated submission snapshot created: ",
  OUT
)

# =============================================================================
# STEP 12 - Create a ZIP archive of the validated snapshot
# =============================================================================

zip_file <- paste0(OUT, ".zip")

zip::zipr(
  zipfile = zip_file,
  files = list.files(
    OUT,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  ),
  root = OUT
)

message("Reviewer ZIP created: ", zip_file)
