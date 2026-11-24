# ============================================================================
# CHOPIN | 05 - Sampling inventory for Supplementary Materials
# One record per sample_TAG x analysis type x contaminant family.
# Sources: {chopin.data} raw datasets; no inference from sample_TAG identifiers.
# Run AFTER scripts 01 and 02, from the root of the chopin.BMF project.
# ============================================================================

suppressPackageStartupMessages({
  library(chopin.data)
  library(tidyverse)
  library(lubridate)
})

OUT <- "manuscript/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

FOCUS_AGE <- "G0"              # NULL = include all age classes
CHEM_SCOPE <- "BMF"            # "BMF" = samples/families in scenarios.rds; "all" = raw assays
DIET_SCOPE <- "analysed"       # "analysed" = same prey filters as script 02; "all" = all G0 records
INCLUDE_BIOMETRICS <- FALSE    # TRUE = add separate biometric-analysis rows

for (nm in c("soles_contam", "benthos_contam", "soles_stomach",
             "soles_G0_pools_biometrics", "contam_info", "species_classif")) {
  data(list = nm, package = "chopin.data", envir = environment())
  if (!exists(nm, inherits = FALSE)) stop("Dataset missing: ", nm)
}

# ---- 1. Helpers: metadata are read, never guessed from sample IDs ----------

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NULL")] <- NA_character_
  x
}

get_field <- function(d, candidates) {
  hit <- intersect(candidates, names(d))
  if (length(hit)) clean_chr(d[[hit[1]]]) else rep(NA_character_, nrow(d))
}

parse_sample_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x <- clean_chr(x)
  out <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    out[bad] <- as.Date(lubridate::parse_date_time(
      x[bad], orders = c("ymd", "dmy", "mdy", "ymd HMS", "dmy HMS"),
      quiet = TRUE
    ))
  }
  out
}

collapse_values <- function(x) {
  x <- sort(unique(clean_chr(x)))
  x <- x[!is.na(x)]
  if (length(x)) paste(x, collapse = "; ") else NA_character_
}

filter_age <- function(d) {
  age_col <- intersect(c("class", "age_class", "age_classification"), names(d))
  if (is.null(FOCUS_AGE) || !length(age_col)) return(d)
  d[clean_chr(d[[age_col[1]]]) == FOCUS_AGE &
      !is.na(d[[age_col[1]]]), , drop = FALSE]
}

sample_metadata <- function(d, source) {
  date_col <- intersect(c("date", "sampling_date", "sample_date", "Date"), names(d))
  dates <- if (length(date_col)) parse_sample_date(d[[date_col[1]]]) else
    as.Date(rep(NA_character_, nrow(d)))
  raw_season <- get_field(d, c("season", "saison"))
  raw_campaign <- get_field(d, c("campaign", "campagne", "sampling_campaign",
                                  "campaign_id", "sampling_event"))
  # When no season is recorded, reproduce the operational seasonal split of
  # script 02 (August-December = autumn; otherwise spring).
  derived_season <- case_when(
    is.na(dates) ~ NA_character_,
    lubridate::month(dates) %in% 8:12 ~ "autumn",
    TRUE ~ "spring"
  )
  # YYYY-MM is a DATE-DERIVED grouping, NOT a recorded campaign identifier.
  derived_campaign <- if_else(is.na(dates), NA_character_, format(dates, "%Y-%m"))

  tibble(
    .row_id = seq_len(nrow(d)),
    sample_TAG = get_field(d, c("sample_TAG", "sample_tag", "sample_id")),
    date = as.character(dates),
    season = coalesce(raw_season, derived_season),
    season_source = case_when(
      !is.na(raw_season) ~ "Recorded in source",
      !is.na(derived_season) ~ "Derived using script 02 month rule",
      TRUE ~ "Missing"
    ),
    age_class = get_field(d, c("class", "age_class", "age_classification")),
    campaign = coalesce(raw_campaign, derived_campaign),
    campaign_source = case_when(
      !is.na(raw_campaign) ~ "Recorded in source",
      !is.na(derived_campaign) ~ "Derived from date (YYYY-MM)",
      TRUE ~ "Missing"
    ),
    source_dataset = source
  )
}

summarise_inventory <- function(x) {
  x |>
    filter(!is.na(sample_TAG)) |>
    group_by(sample_TAG, analysis_type, contaminant_family,
             source_dataset) |>
    summarise(
      date = collapse_values(date),
      season = collapse_values(season),
      season_source = collapse_values(season_source),
      age_class = collapse_values(age_class),
      campaign = collapse_values(campaign),
      campaign_source = collapse_values(campaign_source),
      n_compounds = if (all(is.na(n_compounds))) NA_integer_
                    else max(n_compounds, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- 2. Determine chemical families from actually populated measurements ---
# Both quantified concentrations and the recorded value 0 (censored marker)
# count as analytical results. Merely having a chemical column is insufficient.

chemical_dictionary <- contam_info |>
  filter(family %in% c("PCB", "HBCDD", "PFAS")) |>
  distinct(chemical, family)

chemical_inventory <- function(d, source, analysis_type) {
  d <- if (source == "soles_contam") filter_age(d) else d
  dictionary <- tibble(column = names(d)) |>
    filter(str_detect(column, "_ng_g(lw|dw|ww)(_censored)?$")) |>
    mutate(chemical = str_remove(column, "_ng_g(lw|dw|ww)(_censored)?$")) |>
    inner_join(chemical_dictionary, by = "chemical", relationship = "many-to-one")
  if (!nrow(dictionary)) stop("No known chemical measurement columns in ", source)

  chemical_presence <- d |>
    mutate(.row_id = row_number()) |>
    select(.row_id, all_of(dictionary$column)) |>
    pivot_longer(-.row_id, names_to = "column", values_to = "reading",
                 values_transform = list(reading = as.numeric)) |>
    inner_join(dictionary, by = "column", relationship = "many-to-one") |>
    filter(!is.na(reading)) |>
    group_by(.row_id, family) |>
    summarise(n_compounds = n_distinct(chemical), .groups = "drop")

  chemical_presence |>
    left_join(sample_metadata(d, source), by = ".row_id",
              relationship = "many-to-one") |>
    transmute(sample_TAG, date, season, season_source, age_class, campaign,
              campaign_source, source_dataset, analysis_type,
              contaminant_family = family, n_compounds) |>
    summarise_inventory()
}

sole_chemical <- chemical_inventory(
  soles_contam, "soles_contam", "Sole chemistry"
)
prey_chemical <- chemical_inventory(
  benthos_contam, "benthos_contam", "Benthic chemistry"
)

# ---- 3. Add stomach-contents and (optionally) biometric analyses -----------
# In soles_stomach, 'species' denotes the PREY taxon, not the sole's identity.

stomach_data <- filter_age(soles_stomach)
if (DIET_SCOPE == "analysed") {
  # Exact inclusion logic in script 02: G0, Aug-Dec, positive prey count,
  # non-Chordata taxa, observed in at least three fish.
  stomach_data <- stomach_data |>
    mutate(.month = lubridate::month(parse_sample_date(date)),
           .season_analysis = if_else(.month %in% 8:12, "autumn", "spring"),
           N_tract = suppressWarnings(as.numeric(N_tract))) |>
    filter(.season_analysis == "autumn", !is.na(N_tract), N_tract > 0)
  diet_taxa <- stomach_data |>
    distinct(species) |>
    left_join(species_classif, by = c("species" = "scientific_name")) |>
    filter(!is.na(Phylum), Phylum != "Chordata") |>
    pull(species)
  stomach_data <- stomach_data |>
    filter(species %in% diet_taxa) |>
    group_by(species) |>
    filter(n_distinct(sample_TAG) >= 3) |>
    ungroup()
} else if (DIET_SCOPE != "all") {
  stop("DIET_SCOPE must be 'analysed' or 'all'")
}

stomach_inventory <- sample_metadata(stomach_data, "soles_stomach") |>
  mutate(analysis_type = "Stomach contents",
         contaminant_family = NA_character_, n_compounds = NA_integer_) |>
  summarise_inventory()

# Biometrics are also used to recover sole sampling dates if the chemical
# table has no date column. Multiple dates in one pool stay visible.
biometrics_metadata <- sample_metadata(
  soles_G0_pools_biometrics, "soles_G0_pools_biometrics"
) |>
  filter(!is.na(sample_TAG)) |>
  group_by(sample_TAG) |>
  summarise(date_bio = collapse_values(date),
            season_bio = collapse_values(season), .groups = "drop")

# Explicitly flag discrepancies rather than silently resolving them.
qa_sole_pool_seasons <- sole_chemical |>
  left_join(biometrics_metadata, by = "sample_TAG", relationship = "many-to-one") |>
  filter(!is.na(season), !is.na(season_bio), season != season_bio) |>
  select(sample_TAG, season_chemical = season, season_biometrics = season_bio,
         date_chemical = date, date_biometrics = date_bio)
readr::write_csv(qa_sole_pool_seasons,
                 file.path(OUT, "QA_sole_season_disagreements.csv"))

sole_chemical <- sole_chemical |>
  left_join(biometrics_metadata, by = "sample_TAG", relationship = "many-to-one") |>
  mutate(
    date = coalesce(date, date_bio),
    season_source = if_else(is.na(season) & !is.na(season_bio),
                            "Recovered from sole pool biometrics", season_source),
    season = coalesce(season, season_bio),
    campaign_from_bio = if_else(!is.na(date) & !str_detect(date, fixed(";")),
                                str_sub(date, 1, 7), NA_character_),
    campaign_source = if_else(is.na(campaign) & !is.na(campaign_from_bio),
                              "Derived from biometric sampling date", campaign_source),
    campaign = coalesce(campaign, campaign_from_bio)
  ) |>
  select(-date_bio, -season_bio, -campaign_from_bio)

# If documenting the BMF calculation, keep only chemistry sample/family
# combinations actually used by script 01 (rather than all raw assays).
chemical_rows <- bind_rows(sole_chemical, prey_chemical)
if (CHEM_SCOPE == "BMF") {
  if (!file.exists("results/scenarios.rds")) {
    stop("Run script 01 first: results/scenarios.rds is required for CHEM_SCOPE='BMF'")
  }
  bmf_keys <- readRDS("results/scenarios.rds") |>
    transmute(
      sample_TAG = as.character(sample_id),
      analysis_type = if_else(compartment == "Sole", "Sole chemistry", "Benthic chemistry"),
      contaminant_family = family
    ) |>
    distinct()
  chemical_rows <- chemical_rows |>
    semi_join(bmf_keys, by = c("sample_TAG", "analysis_type", "contaminant_family"))
} else if (CHEM_SCOPE != "all") {
  stop("CHEM_SCOPE must be 'BMF' or 'all'")
}

if (!nrow(chemical_rows)) stop("No chemical samples matched the requested scope")
if (!nrow(stomach_inventory)) warning("No stomach samples matched DIET_SCOPE")

inventory <- bind_rows(chemical_rows, stomach_inventory)
if (INCLUDE_BIOMETRICS) {
  biometric_inventory <- sample_metadata(
    soles_G0_pools_biometrics, "soles_G0_pools_biometrics"
  ) |>
    mutate(analysis_type = "Sole biometrics",
           contaminant_family = NA_character_, n_compounds = NA_integer_) |>
    summarise_inventory()
  inventory <- bind_rows(inventory, biometric_inventory)
}

# ---- 4. Export detailed and compact Supplementary tables ------------------

inventory <- inventory |>
  arrange(analysis_type, campaign, date, sample_TAG, contaminant_family)

inventory_compact <- inventory |>
  group_by(sample_TAG, analysis_type, source_dataset) |>
  summarise(
    date = collapse_values(date),
    season = collapse_values(season),
    age_class = collapse_values(age_class),
    campaign = collapse_values(campaign),
    contaminant_family = collapse_values(contaminant_family),
    n_families = n_distinct(contaminant_family, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(analysis_type, campaign, date, sample_TAG)

# Requested seven columns: one row per sample ID x type of analysis.
table_supp <- inventory_compact |>
  transmute(
    `Sample ID` = sample_TAG,
    `Sampling date(s)` = date,
    `Season` = season,
    `Age class` = age_class,
    `Campaign` = campaign,
    `Analysis type` = analysis_type,
    `Contaminant families` = replace_na(contaminant_family, "Not applicable")
  )

readr::write_csv(table_supp, file.path(OUT, "Table_Sx_sampling_inventory.csv"), na = "")
readr::write_csv(inventory, file.path(OUT, "Dataset_Sx_sampling_inventory_detailed.csv"), na = "")
saveRDS(inventory, "results/sample_inventory.rds")

# QA: do not conceal missing or conflicting source metadata.
qa_metadata <- inventory |>
  filter(is.na(date) | is.na(campaign) | str_detect(date, fixed(";")) |
           str_detect(season, fixed(";")) | str_detect(campaign, fixed(";")))
readr::write_csv(qa_metadata, file.path(OUT, "QA_sampling_metadata_to_check.csv"))

qa_counts <- inventory |>
  distinct(sample_TAG, source_dataset, analysis_type, contaminant_family) |>
  count(source_dataset, analysis_type, contaminant_family, name = "n_sample_TAG")
readr::write_csv(qa_counts, file.path(OUT, "QA_sampling_counts.csv"))

# Optional Word table. Exact source detail remains in the CSV files.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- flextable::flextable(table_supp) |>
    flextable::theme_booktabs() |>
    flextable::fontsize(size = 8, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
  flextable::save_as_docx(
    "Sampling inventory" = ft,
    path = file.path(OUT, "Table_Sx_sampling_inventory.docx")
  )
}

print(qa_counts, n = Inf)
message("Rows in compact table: ", nrow(table_supp))
message("Rows needing metadata verification: ", nrow(qa_metadata))
message("Sole chemistry/biometric season disagreements: ", nrow(qa_sole_pool_seasons))
message("Outputs written to: ", OUT)
