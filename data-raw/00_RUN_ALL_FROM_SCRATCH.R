# CHOPIN | Regenerate all package data and manuscript results from a clean checkout.
# Run at the root of chopin.BMF. No derived data/ or results/ are prerequisites.
# Requires the source datasets in installed {chopin.data} and the Gogina XLSX.

if (!file.exists("DESCRIPTION") ||
    !any(readLines("DESCRIPTION") == "Package: chopin.BMF")) {
  stop("Run this script from the chopin.BMF package root")
}

steps <- c(
  "00_CHOPIN_preflight.R",
  "01_CHOPIN_contamination_gamma.R",
  "02_CHOPIN_diet.R",
  "03_CHOPIN_BMF_gamma.R",
  "03b_CHOPIN_zero_safe_sensitivity.R",
  "04_CHOPIN_gamma_vs_lognormal.R",
  "05_CHOPIN_sample_inventory_supplementary.R",
  "06_CHOPIN_BMF_mass_quantiles.R",
  "07_CHOPIN_leave_one_prey_out.R",
  "08_CHOPIN_prey_season_sensitivity.R",
  "09_CHOPIN_HBCDD_isomer_bounds.R"
)
run_all <- function(include_snapshot = FALSE) {
  for (step in steps) {
    path <- file.path("data-raw", step)
    if (!file.exists(path)) stop("Missing pipeline script: ", path)
    message("\n>>> Running ", step)
    # Independent environments prevent accidental reliance on stale R objects.
    source(path, local = new.env(parent = globalenv()), echo = FALSE)
    message("<<< Finished ", step)
  }
  message("\n>>> Post-run source/results consistency check")
  source("data-raw/10_CHOPIN_postrun_provenance.R",
         local = new.env(parent = globalenv()), echo = FALSE)
  message("<<< Consistency check passed")
  if (isTRUE(include_snapshot)) {
    message("\n>>> Creating submission snapshot")
    source("data-raw/11_CHOPIN_checked_submission_snapshot.R",
           local = new.env(parent = globalenv()), echo = FALSE)
  }
  message("\nCHOPIN pipeline completed. Inspect results/, data/, and manuscript/.")
  invisible(TRUE)
}

run_all(include_snapshot = FALSE)  # Change to TRUE only for a frozen submission run.
