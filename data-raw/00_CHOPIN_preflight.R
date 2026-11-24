# CHOPIN | 00 - PREFLIGHT: only SOURCE data are required.
# No results/*.rds, data/*.rda, figures or tables are needed at this stage.
# Run from the project root, before scripts 01-09.
source("data-raw/00_CHOPIN_config.R", local = TRUE)

if (!file.exists("DESCRIPTION") ||
    !any(readLines("DESCRIPTION") == "Package: chopin.BMF")) {
  stop("Open the chopin.BMF project root before running the scripts.")
}
required_files <- c(
  "data-raw/01_CHOPIN_contamination_gamma.R",
  "data-raw/02_CHOPIN_diet.R",
  "data-raw/03_CHOPIN_BMF_gamma.R",
  "data-raw/04_CHOPIN_gamma_vs_lognormal.R",
  "inst/extdata/Gogina_data-2021-0002-01.xlsx"
)
if (any(!file.exists(required_files))) {
  stop("Missing SOURCE files: ",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}

required_packages <- c(
  "chopin.data", "tidyverse", "Btrait", "readxl", "janitor",
  "survival", "usethis", "ggrepel", "scales", "tidytext", "lubridate"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop("Install these R packages first: ", paste(missing_packages, collapse = ", "),
       ". Source package chopin.data is installed separately; see README.")
}
if (getRversion() < "4.1.0" || utils::packageVersion("dplyr") < "1.1.0") {
  stop("This pipeline requires R >= 4.1 and dplyr >= 1.1.")
}
if (!SOLE_SEASON %in% c("spring", "autumn")) {
  stop("Set SOLE_SEASON in data-raw/00_CHOPIN_config.R.")
}

suppressPackageStartupMessages(library(dplyr))
for (nm in c("soles_contam", "benthos_contam", "soles_stomach",
             "soles_G0_pools_biometrics", "species_classif", "contam_info")) {
  e <- new.env(parent = emptyenv())
  utils::data(list = nm, package = "chopin.data", envir = e)
  if (!exists(nm, envir = e, inherits = FALSE)) {
    stop("Missing SOURCE dataset in chopin.data: ", nm)
  }
}
e <- new.env(parent = emptyenv())
utils::data("soles_contam", package = "chopin.data", envir = e)
sole_seasons <- e$soles_contam |>
  filter(class == "G0") |>
  group_by(season) |>
  summarise(n_pools = n_distinct(sample_TAG), .groups = "drop")
print(sole_seasons)
if (!any(sole_seasons$season == SOLE_SEASON & sole_seasons$n_pools > 0)) {
  stop("No age-0 sole chemistry in SOLE_SEASON=", SOLE_SEASON)
}
if (isTRUE(SAVE_TO_PACKAGE) && !requireNamespace("usethis", quietly = TRUE)) {
  stop("usethis required when SAVE_TO_PACKAGE = TRUE")
}
for (d in c("data", "results", "manuscript/figures", "manuscript/tables")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
write.csv(sole_seasons, "results/source_sole_seasons.csv", row.names = FALSE)
writeLines(c(
  paste0("SOLE_SEASON=", SOLE_SEASON),
  paste0("chopin.data_version=", as.character(utils::packageVersion("chopin.data"))),
  paste0("R_version=", as.character(getRversion())),
  paste0("SAVE_TO_PACKAGE=", SAVE_TO_PACKAGE)
), "results/source_configuration.txt")
message("00 preflight passed: raw data and tools available; no prior results required.")
