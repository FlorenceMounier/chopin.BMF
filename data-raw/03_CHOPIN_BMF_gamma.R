# ================================================================
# CHOPIN | 03 - BMFs, censoring, dietary methods and sensitivities
# Run AFTER scripts 01_CHOPIN_contamination_gamma.R and 02_CHOPIN_diet.R
# ================================================================
# Requires data object 'scenarios' prepared in contamination script 01.
# Does not recalculate censoring imputation or alter analytical data.
# Core outputs: BMF by prey, phylum, dietary method; 4 censoring
# scenarios; 3 diet weightings; pairwise comparisons; sensitivity
# hierarchy; bivalve siphon sensitivity (2-5%); manuscript figures.
# ================================================================

suppressPackageStartupMessages({
  library(chopin.data)
  library(tidyverse)
  library(ggrepel)
  library(scales)
})

OUTPUT_DIR <- "results"
FIGURE_DIR <- "manuscript/figures"
REFERENCE_SIPHON <- 0.0375
SIPHON_FACTORS <- sort(unique(c(seq(0.02, 0.05, by = 0.005), REFERENCE_SIPHON)))
# Main reference for Fig.3 and siphon sensitivities is now GAMMA-MLE.
# Change to "loq_half" only to reproduce the older reference figures.
FOCUS_SCENARIO <- "mle"
source("data-raw/00_CHOPIN_config.R", local = TRUE)
for (d in c(OUTPUT_DIR, FIGURE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Load the output of script 02 -- independent of what is in workspace.
prop_diet <- readRDS(file.path(OUTPUT_DIR, "prop_diet.rds"))
diet_sensitivity <- readRDS(file.path(OUTPUT_DIR, "diet_sensitivity.rds"))

# Read the EXACT gamma-MLE scenarios created by script 01; do not use
# data("scenarios", package=...) because an installed package can be stale.
scenario_file <- file.path(OUTPUT_DIR, "scenarios.rds")
if (!file.exists(scenario_file)) stop("Run script 01 first: ", scenario_file)
scenarios <- readRDS(scenario_file)
if (!nrow(scenarios)) stop("No analytical scenarios in ", scenario_file)

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}
safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}
safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else unname(quantile(x, p))
}
safe_extreme <- function(x, fun = min) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else fun(x)
}
safe_ratio <- function(num, den) {
  if_else(is.finite(num) & is.finite(den) & den > 0,
          num / den, NA_real_)
}
check_unique <- function(data, by) {
  duplicates <- data |> count(across(all_of(by))) |> filter(n > 1)
  if (nrow(duplicates)) stop("Duplicate keys: ", paste(by, collapse = ", "))
  invisible(TRUE)
}
write_result <- function(x, name) {
  write_csv(x, file.path(OUTPUT_DIR, paste0(name, ".csv")))
  saveRDS(x, file.path(OUTPUT_DIR, paste0(name, ".rds")))
}

# ---- 1. Chemistry data: four censoring scenarios -----------------
value_cols <- c("value_zero", "value_loq_half", "value_loq", "value_mle")
required_cols <- c("family", "compound_name", "unit", "species", "sample_id", value_cols)
missing_cols <- setdiff(required_cols, names(scenarios))
if (length(missing_cols)) stop("Missing columns in scenarios: ",
                              paste(missing_cols, collapse = ", "))

stopifnot(FOCUS_SCENARIO %in% c("zero", "loq_half", "loq", "mle"))
cat("03 analytical input:", nrow(scenarios), "sample-compound records;",
    "MLE = gamma from script 01\n")

chem_data <- scenarios |>
  pivot_longer(all_of(value_cols), names_to = "scenario", values_to = "value") |>
  mutate(scenario = str_remove(scenario, "^value_"),
         value = if_else(is.finite(value), value, NA_real_))

taxa_mapping <- tribble(
  ~taxa,          ~Phylum,
  "Polychaetes",  "Annelida",
  "Crustaceans",   "Arthropoda",
  "Bivalves",     "Mollusca"
)
sole_data <- chem_data |> filter(species == "Solea_solea")
prey_data <- chem_data |>
  filter(species != "Solea_solea") |>
  left_join(taxa_mapping, by = "taxa", relationship = "many-to-one")

diet_data <- prop_diet |>
  transmute(Phylum, p_numeric = p_numeric_phylum,
            p_biomass = p_biomass_phylum)

# Derive censoring audit from the same analytical data; no additional
# package object required. Collapse to one record per biological pool.
censor_labels <- scenarios |>
  mutate(group = if_else(species == "Solea_solea", "sole", "prey")) |>
  distinct(group, family, compound_name, species, sample_id, censored) |>
  group_by(group, family, compound_name) |>
  summarise(pct_censored = 100 * mean(censored %in% TRUE),
            .groups = "drop") |>
  pivot_wider(names_from = group, values_from = pct_censored) |>
  mutate(
    censor_max = pmax(sole, prey, na.rm = TRUE),
    censor_max = if_else(is.infinite(censor_max), NA_real_, censor_max)
  )

# ---- 2. Generic BMF function ------------------------------------
# Equal: phyla equally weighted (NOT prey species equally weighted).
# Within each phylum: each chemically measured species weighted equally.
# For any compound: renormalize on phyla with finite prey concentration.
calc_bmfs <- function(sole_data, prey_data, diet_data) {
  keys <- c("scenario", "family", "compound_name", "unit")
  check_unique(diet_data, "Phylum")

  sole_mean <- sole_data |>
    group_by(across(all_of(keys))) |>
    summarise(c_sole = safe_mean(value),
              n_sole = sum(is.finite(value)), .groups = "drop")

  prey_mean <- prey_data |>
    group_by(across(all_of(keys)), species, Phylum) |>
    summarise(c_prey = safe_mean(value),
              n_prey_samples = sum(is.finite(value)), .groups = "drop")

  by_prey <- prey_mean |>
    inner_join(sole_mean, by = keys, relationship = "many-to-one") |>
    mutate(bmf_prey = safe_ratio(c_sole, c_prey))

  by_phylum <- prey_mean |>
    filter(is.finite(c_prey)) |>
    group_by(across(all_of(c(keys, "Phylum")))) |>
    summarise(c_phylum = mean(c_prey),
              n_species = n_distinct(species), .groups = "drop") |>
    inner_join(diet_data, by = "Phylum", relationship = "many-to-one")

  by_compound <- by_phylum |>
    group_by(across(all_of(keys))) |>
    summarise(
      n_phyla = n_distinct(Phylum),
      n_species = sum(n_species),
      coverage_numeric = sum(p_numeric),
      coverage_biomass = sum(p_biomass),
      c_equal = mean(c_phylum),
      c_numeric = if (sum(p_numeric) > 0) {
        sum(p_numeric * c_phylum) / sum(p_numeric)
      } else NA_real_,
      c_biomass = if (sum(p_biomass) > 0) {
        sum(p_biomass * c_phylum) / sum(p_biomass)
      } else NA_real_,
      .groups = "drop"
    ) |>
    inner_join(sole_mean, by = keys, relationship = "one-to-one")

  comparison <- by_compound |>
    pivot_longer(c(c_equal, c_numeric, c_biomass),
                 names_to = "diet_method", values_to = "c_diet") |>
    mutate(diet_method = str_remove(diet_method, "^c_"),
           bmf = safe_ratio(c_sole, c_diet),
           threshold_above_one = if_else(is.finite(bmf), bmf > 1, NA))
  list(prey = by_prey, phylum = by_phylum, comparison = comparison)
}


# ---- 3. Main 4 x 3 BMF analyses ---------------------------------
results <- calc_bmfs(sole_data, prey_data, diet_data)
bmf_prey <- results$prey
bmf_phylum <- results$phylum
bmf_compare <- results$comparison

# Stability of BMF classification across 12 scenarios

bmf_threshold <- bmf_compare |>
  group_by(family, compound_name, unit) |>
  summarise(
    n_scenarios = n(),
    n_valid = sum(is.finite(bmf)),

    n_above_1 = sum(bmf > 1, na.rm = TRUE),
    n_below_1 = sum(bmf < 1, na.rm = TRUE),
    n_equal_1 = sum(bmf == 1, na.rm = TRUE),

    bmf_min = if (all(!is.finite(bmf))) NA_real_
    else min(bmf[is.finite(bmf)]),

    bmf_max = if (all(!is.finite(bmf))) NA_real_
    else max(bmf[is.finite(bmf)]),

    .groups = "drop"
  ) |>
  mutate(
    # TRUE if at least one scenario is below 1
    # and another is above 1
    flip = n_above_1 > 0 & n_below_1 > 0,

    classification = case_when(
      n_valid != 12 ~ "Incomplete",
      flip ~ "Method-dependent",
      n_above_1 == 12 ~ "Always > 1",
      n_below_1 == 12 ~ "Always < 1",
      n_equal_1 == 12 ~ "Always = 1",
      TRUE ~ "Includes BMF = 1"
    )
  )

bmf_threshold <- bmf_threshold |> select(-unit, -n_scenarios, -n_valid)


bmf_threshold |>
  filter(flip) |>
  select(
    family,
    compound_name,
    n_above_1,
    n_below_1,
    bmf_min,
    bmf_max
  )

keys <- c("scenario", "family", "compound_name", "unit")
check_unique(bmf_compare, c(keys, "diet_method"))

bmf_effects <- bmf_compare |>
  select(all_of(keys), diet_method, bmf) |>
  pivot_wider(names_from = diet_method, values_from = bmf,
              names_prefix = "bmf_") |>
  mutate(
    delta_numeric_equal = if_else(bmf_numeric > 0 & bmf_equal > 0,
                                  log2(bmf_numeric / bmf_equal), NA_real_),
    delta_biomass_equal = if_else(bmf_biomass > 0 & bmf_equal > 0,
                                  log2(bmf_biomass / bmf_equal), NA_real_),
    delta_biomass_numeric = if_else(bmf_biomass > 0 & bmf_numeric > 0,
                                    log2(bmf_biomass / bmf_numeric), NA_real_),
    delta_biomass_numeric_pct = if_else(bmf_numeric > 0 & is.finite(bmf_biomass),
                                        100 * (bmf_biomass / bmf_numeric - 1), NA_real_),
    flip_biomass_numeric = is.finite(bmf_numeric) & is.finite(bmf_biomass) &
      ((bmf_numeric > 1 & bmf_biomass < 1) |
       (bmf_numeric < 1 & bmf_biomass > 1))
  ) |>
  left_join(
    bmf_compare |>
      distinct(across(all_of(keys)), n_phyla, n_species,
               coverage_numeric, coverage_biomass),
    by = keys, relationship = "one-to-one"
  )

family_summary <- bmf_effects |>
  group_by(scenario, family) |>
  summarise(n_compounds = sum(is.finite(delta_biomass_numeric)),
            median_delta = safe_median(delta_biomass_numeric),
            n_flip = sum(flip_biomass_numeric, na.rm = TRUE),
            .groups = "drop")

# ---- 4. Predator-prey concentration ratio distribution -----------
bmf_pairs <- bmf_prey |>
  filter(!is.na(Phylum), Phylum %in% diet_data$Phylum) |>
  mutate(bmf_pair = safe_ratio(c_sole, c_prey),
         pair_valid = is.finite(bmf_pair))

bmf_pair_summary <- bmf_pairs |>
  group_by(across(all_of(keys))) |>
  summarise(
    n_species_measured = n_distinct(species[is.finite(c_prey)]),
    n_species_valid = n_distinct(species[pair_valid]),
    n_phyla = n_distinct(Phylum[pair_valid]),
    bmf_pair_median = safe_median(bmf_pair),
    bmf_pair_q25 = safe_quantile(bmf_pair, .25),
    bmf_pair_q75 = safe_quantile(bmf_pair, .75),
    bmf_pair_min = safe_extreme(bmf_pair),
    bmf_pair_max = safe_extreme(bmf_pair, max), .groups = "drop")

bmf_pair_diet_compare <- bmf_compare |>
  left_join(bmf_pair_summary, by = keys, relationship = "many-to-one") |>
  mutate(
    delta_pair_diet_pct = if_else(bmf_pair_median > 0 & is.finite(bmf),
                                   100 * (bmf / bmf_pair_median - 1), NA_real_),
    delta_pair_diet_log2 = if_else(bmf > 0 & bmf_pair_median > 0,
                                    log2(bmf / bmf_pair_median), NA_real_),
    diet_position = case_when(
      !is.finite(bmf) | !is.finite(bmf_pair_q25) |
        !is.finite(bmf_pair_q75) ~ NA_character_,
      bmf < bmf_pair_q25 ~ "Below Q25",
      bmf > bmf_pair_q75 ~ "Above Q75",
      TRUE ~ "Within IQR"))

pair_rank <- bmf_pairs |>
  filter(is.finite(bmf_pair)) |>
  inner_join(bmf_compare |> select(all_of(keys), diet_method, bmf),
             by = keys, relationship = "many-to-many") |>
  group_by(across(all_of(keys)), diet_method) |>
  summarise(pct_pairs_below_diet = if (all(!is.finite(bmf))) NA_real_
            else 100 * mean(bmf_pair < bmf), .groups = "drop")
bmf_pair_diet_compare <- bmf_pair_diet_compare |>
  left_join(pair_rank, by = c(keys, "diet_method"), relationship = "one-to-one")

# ---- 5. Sensitivity hierarchy -----------------------------------
diet_methods <- c("equal", "numeric", "biomass")
censor_methods <- c("zero", "loq_half", "loq", "mle")
compound_keys <- c("family", "compound_name", "unit")

sensitivity_data <- bmf_effects |>
  left_join(censor_labels, by = c("family", "compound_name"),
            relationship = "many-to-one") |>
  mutate(scenario = factor(scenario, levels = censor_methods))

bmf_long <- sensitivity_data |>
  select(all_of(compound_keys), scenario, bmf_equal, bmf_numeric,
         bmf_biomass, censor_max) |>
  pivot_longer(c(bmf_equal, bmf_numeric, bmf_biomass),
               names_to = "diet_method", names_prefix = "bmf_",
               values_to = "bmf") |>
  mutate(scenario = as.character(scenario))
check_unique(bmf_long, c(compound_keys, "scenario", "diet_method"))


# Completeness of the 12 BMF combinations
completeness <- bmf_long |>
  group_by(across(all_of(compound_keys))) |>
  summarise(
    n_combinations = n(),
    n_nonnegative = sum(is.finite(bmf) & bmf >= 0),
    n_positive = sum(is.finite(bmf) & bmf > 0),
    complete_nonnegative =
      n_combinations == 12 & n_nonnegative == 12,
    complete_positive =
      n_combinations == 12 & n_positive == 12,
    contains_zero = any(is.finite(bmf) & bmf == 0),
    .groups = "drop"
  ) |>
  mutate(
    sensitivity_status = case_when(
      complete_positive ~ "Complete: log-ratio eligible",
      complete_nonnegative & contains_zero ~
        "Complete: zero BMF prevents finite log-ratios",
      TRUE ~ "Incomplete or invalid BMF values"
    )
  )

write_result(completeness, "sensitivity_completeness")

# All valid BMFs, including zero values
bmf_all_valid <- bmf_long |>
  inner_join(
    completeness |>
      filter(complete_nonnegative) |>
      select(all_of(compound_keys)),
    by = compound_keys
  )

# Strictly positive BMFs for log2 sensitivity calculations
bmf_complete <- bmf_long |>
  inner_join(
    completeness |>
      filter(complete_positive) |>
      select(all_of(compound_keys)),
    by = compound_keys
  )


bmf_family <- bmf_all_valid |>
  filter(scenario == FOCUS_SCENARIO & diet_method == "biomass") |>
  group_by(family) |>
  summarise(
    min_BMF = min(bmf),
    median_BMF = median(bmf),
    max_BMF = max(bmf),
    .groups = "drop"
    )

sensitivity_diet <- bmf_complete |>
  group_by(across(all_of(compound_keys)), scenario) |>
  summarise(S_diet = log2(max(bmf) / min(bmf)), .groups = "drop")
sensitivity_censor <- bmf_complete |>
  group_by(across(all_of(compound_keys)), diet_method) |>
  summarise(S_censor = log2(max(bmf) / min(bmf)), .groups = "drop")

diet_summary <- sensitivity_diet |>
  group_by(across(all_of(compound_keys))) |>
  summarise(S_diet = median(S_diet), S_diet_max = max(S_diet), .groups = "drop")
censor_summary <- sensitivity_censor |>
  group_by(across(all_of(compound_keys))) |>
  summarise(S_censor = median(S_censor), S_censor_max = max(S_censor),
            .groups = "drop")

sensitivity_hierarchy <- diet_summary |>
  inner_join(censor_summary, by = compound_keys, relationship = "one-to-one") |>
  left_join(sensitivity_data |>
              distinct(across(all_of(compound_keys)), censor_max),
            by = compound_keys, relationship = "one-to-one") |>
  mutate(dominant_source = case_when(
    S_diet > S_censor ~ "Diet weighting",
    S_diet < S_censor ~ "Censoring treatment",
    TRUE ~ "Equal sensitivity"),
    fold_diet = 2^S_diet, fold_censor = 2^S_censor,
    delta_S = S_diet - S_censor) |>
  arrange(desc(abs(delta_S)))

# Targeted sensitivity: numeric vs biomass (not confounded by equal weighting)
sensitivity_numeric_biomass <- bmf_complete |>
  filter(diet_method %in% c("numeric", "biomass")) |>
  group_by(across(all_of(compound_keys)), scenario) |>
  summarise(S_numeric_biomass = log2(max(bmf) / min(bmf)), .groups = "drop") |>
  group_by(across(all_of(compound_keys))) |>
  summarise(S_numeric_biomass = median(S_numeric_biomass), .groups = "drop")
sensitivity_hierarchy <- sensitivity_hierarchy |>
  left_join(sensitivity_numeric_biomass, by = compound_keys,
            relationship = "one-to-one")

family_sensitivity <- sensitivity_hierarchy |>
  group_by(family) |>
  summarise(n_compounds = n(),
            median_S_diet = median(S_diet),
            median_S_censor = median(S_censor),
            median_S_numeric_biomass = median(S_numeric_biomass),
            n_diet_dominant = sum(dominant_source == "Diet weighting"),
            n_censor_dominant = sum(dominant_source == "Censoring treatment"),
            .groups = "drop")

# ---- 6. Bivalve siphon correction sensitivity --------------------
# IMPORTANT: use diet_sensitivity from script 02 (NOT sensitivity_diet).
# Function argument deliberately named 'factor_value' to avoid collision
# with the column 'factor_siphon'.
check_unique(diet_sensitivity, c("factor_siphon", "Phylum"))

calc_bmf_siphon <- function(factor_value) {
  diet_test <- diet_sensitivity |>
    filter(dplyr::near(.data$factor_siphon, .env$factor_value)) |>
    transmute(Phylum, p_numeric, p_biomass)
  check_unique(diet_test, "Phylum")
  if (!nrow(diet_test)) stop("No diet at siphon factor ", factor_value)

  calc_bmfs(sole_data, prey_data, diet_test)$comparison |>
    filter(scenario == FOCUS_SCENARIO, diet_method == "biomass") |>
    mutate(factor_siphon = factor_value)
}

bmf_siphon <- purrr::map_dfr(SIPHON_FACTORS, calc_bmf_siphon)
check_unique(bmf_siphon, c("family", "compound_name", "unit", "factor_siphon"))

bmf_reference <- bmf_siphon |>
  filter(near(factor_siphon, REFERENCE_SIPHON)) |>
  select(family, compound_name, unit, bmf_reference = bmf)
check_unique(bmf_reference, c("family", "compound_name", "unit"))

bmf_siphon <- bmf_siphon |>
  left_join(bmf_reference, by = c("family", "compound_name", "unit"),
            relationship = "many-to-one") |>
  mutate(delta_pct = if_else(bmf_reference > 0 & is.finite(bmf),
                              100 * (bmf / bmf_reference - 1), NA_real_))

# Exact methodological consistency: 3.75% reproduces the main BMF.
siphon_check <- bmf_reference |>
  left_join(bmf_compare |>
              filter(scenario == FOCUS_SCENARIO, diet_method == "biomass") |>
              select(family, compound_name, unit, bmf_original = bmf),
            by = c("family", "compound_name", "unit"),
            relationship = "one-to-one") |>
  mutate(abs_diff = abs(bmf_reference - bmf_original))
if (any(!is.finite(siphon_check$abs_diff) |
        siphon_check$abs_diff > 1e-8 * pmax(1, siphon_check$bmf_original))) {
  stop("Reference siphon factor does not reproduce original biomass BMFs")
}

siphon_summary <- bmf_siphon |>
  group_by(family, compound_name, unit) |>
  summarise(bmf_min = safe_extreme(bmf),
            bmf_max = safe_extreme(bmf, max),
            max_abs_delta_pct = safe_extreme(abs(delta_pct), max),
            .groups = "drop") |>
  mutate(fold_range = safe_ratio(bmf_max, bmf_min),
         S_siphon = if_else(fold_range > 0, log2(fold_range), NA_real_))

siphon_family <- siphon_summary |>
  group_by(family) |>
  summarise(n_compounds = n(),
            median_fold = safe_median(fold_range),
            median_S = safe_median(S_siphon),
            median_max_delta_pct = safe_median(max_abs_delta_pct),
            .groups = "drop")

sensitivity_extended <- sensitivity_hierarchy |>
  left_join(siphon_summary |>
              select(family, compound_name, unit, S_siphon),
            by = compound_keys, relationship = "one-to-one")

# ---- 7. Export all analysis tables -------------------------------
export_tables <- list(
  bmf_prey = bmf_prey,
  bmf_phylum = bmf_phylum,
  bmf_compare = bmf_compare,
  bmf_effects = bmf_effects,
  family_summary = family_summary,
  bmf_pairs = bmf_pairs,
  bmf_pair_summary = bmf_pair_summary,
  bmf_pair_diet_compare = bmf_pair_diet_compare,
  bmf_threshold = bmf_threshold,
  bmf_family = bmf_family,
  sensitivity_completeness = completeness,
  sensitivity_numeric_biomass = sensitivity_numeric_biomass,
  sensitivity_data = sensitivity_data,
  sensitivity_diet = sensitivity_diet,
  sensitivity_censor = sensitivity_censor,
  sensitivity_hierarchy = sensitivity_hierarchy,
  family_sensitivity = family_sensitivity,
  bmf_siphon = bmf_siphon,
  siphon_summary = siphon_summary,
  siphon_family = siphon_family,
  sensitivity_extended = sensitivity_extended,
  siphon_check = siphon_check,
  censor_labels = censor_labels,
  bmf_pair_rank = pair_rank
)
purrr::iwalk(export_tables, ~ write_result(.x, .y))

if (SAVE_TO_PACKAGE) {
  usethis::use_data(bmf_prey, bmf_phylum, bmf_compare, bmf_effects,
                    family_summary, bmf_pairs, bmf_pair_summary,
                    bmf_pair_diet_compare, sensitivity_data,
                    sensitivity_hierarchy, family_sensitivity,
                    bmf_siphon, siphon_summary, siphon_family,
                    bmf_threshold, bmf_family,
                    overwrite = TRUE)
}



# ---- 8. Manuscript figures --------------------------------------


# ---- 8A. Effect of numerical vs biomass weighting (reference: gamma-MLE) ----
p_delta <- bmf_effects |>
  filter(scenario == FOCUS_SCENARIO, is.finite(delta_biomass_numeric_pct)) |>
  ggplot(aes(x = reorder(compound_name, delta_biomass_numeric_pct),
             y = delta_biomass_numeric_pct)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2.3) + coord_flip() +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Compound", y = "Change in BMF: biomass vs numerical (%)") + theme_bw()
ggsave(file.path(FIGURE_DIR, "Fig.S5_diet_weighting_delta_pct.png"),
       p_delta, width = 11, height = 8, dpi = 300)


# ---- 8B. Censoring sensitivity: biomass-weighted BMF only ----
plot_sens <- sensitivity_data |>
  filter(family != "PCB", is.finite(bmf_biomass), bmf_biomass > 0)
label_sens <- plot_sens |>
  arrange(scenario) |>
  group_by(family, compound_name) |> slice_tail(n = 1) |> ungroup()
p_sensitivity <- ggplot(plot_sens, aes(x = scenario, y = bmf_biomass,
                                       group = compound_name, color = censor_max)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.8) +
  geom_line(linewidth = .8, alpha = .8) + geom_point(size = 2.2) +
  geom_text_repel(data = label_sens, aes(label = compound_name),
                  direction = "y", hjust = 0, nudge_x = .25,
                  size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_y_log10() +
  scale_color_gradient(low = "#2A9D8F", high = "#D95F02",
                       limits = c(0, 100), oob = scales::squish,
                       na.value = "grey70", name = "Maximum censoring (%)") +
  facet_wrap(~ family, scales = "free_y") +
  scale_x_discrete(expand = expansion(add = c(.5, 1.8))) +
  coord_cartesian(clip = "off") +
  labs(x = "Censoring scenario", y = "Biomass-weighted BMF",
       title = "Sensitivity of biomass-weighted BMFs to censoring treatment") +
  theme_bw() + theme(plot.margin = margin(5.5, 90, 5.5, 5.5))
ggsave(file.path(FIGURE_DIR, "Fig.4_censor_sensitivity_biomass.png"),
       p_sensitivity, width = 12, height = 6, dpi = 300)


# ---- 8C. Predator-prey distributions and the three diet-BMF methods ----
pair_plot <- bmf_pairs |>
  filter(scenario == FOCUS_SCENARIO, is.finite(bmf_pair), bmf_pair > 0)
summary_plot <- bmf_pair_summary |>
  filter(scenario == FOCUS_SCENARIO, is.finite(bmf_pair_median),
         bmf_pair_median > 0)
diet_plot <- bmf_compare |>
  filter(scenario == FOCUS_SCENARIO, is.finite(bmf), bmf > 0) |>
  mutate(diet_method = factor(diet_method,
                              levels = c("equal", "numeric", "biomass")))
p_pair_diet <- ggplot() +
  geom_point(data = pair_plot, aes(x = compound_name, y = bmf_pair),
             color = "grey65", alpha = .45,
             position = position_jitter(width = .15, height = 0), size = 1.4) +
  geom_linerange(data = summary_plot,
                 aes(x = compound_name, ymin = bmf_pair_q25, ymax = bmf_pair_q75),
                 color = "black", linewidth = .9) +
  geom_point(data = summary_plot,
             aes(x = compound_name, y = bmf_pair_median),
             color = "black", shape = 18, size = 2.8) +
  geom_point(data = diet_plot, aes(x = compound_name, y = bmf,
                                  color = diet_method),
             position = position_dodge(width = .6), size = 2.3) +
  geom_hline(yintercept = 1, linetype = 2) + scale_y_log10() + coord_flip() +
  facet_wrap(~ family, scales = "free") +
  labs(x = "Compound", y = "Biomagnification factor (BMF)",
       color = "Diet weighting",
       title = paste("Predator-prey versus diet-based BMFs:",
                     if (FOCUS_SCENARIO == "mle") "Gamma-MLE" else FOCUS_SCENARIO)) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIGURE_DIR, "Fig.3_predator_prey_vs_diet.png"),
       p_pair_diet, width = 12, height = 9, dpi = 300)


# ---- 8D. Compound-level methodological hierarchy ----
hierarchy_plot <- sensitivity_hierarchy |>
  filter(family %in% c("HBCDD", "PFAS")) |>
  mutate(
    family = factor(family, levels = c("HBCDD", "PFAS"))
  )
p_hierarchy <- ggplot(
  hierarchy_plot,
  aes(
    x = S_diet,
    y = S_censor,
    color = censor_max
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_text_repel(
    aes(label = compound_name),
    size = 3.2,
    max.overlaps = Inf,
    show.legend = FALSE,
    box.padding = 0.4,
    point.padding = 0.2
  ) +
  scale_color_gradient(
    low = "#2A9D8F",
    high = "#D95F02",
    limits = c(0, 100),
    oob = scales::squish,
    name = "Maximum censoring (%)"
  ) +
  facet_wrap(
    ~ family,
    ncol = 1
  ) +
  coord_equal() +
  labs(
    x = "Sensitivity to diet weighting (log2 fold range)",
    y = "Sensitivity to censoring (log2 fold range)",
    title = "Relative importance of methodological choices",
    subtitle = paste(
      "Above diagonal: censoring dominates;",
      "below: dietary weighting dominates"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_hierarchy

ggsave(file.path(FIGURE_DIR, "Fig.5_methodological_sensitivity_hierarchy.png"),
       p_hierarchy, width = 11, height = 5.5, dpi = 300)


# ---- 8E. Family-level methodological hierarchy ----
p_family <- family_sensitivity |>
  select(family, median_S_diet, median_S_censor) |>
  pivot_longer(-family, names_to = "Source", values_to = "Sensitivity") |>
  mutate(Source = recode(Source,
                         median_S_diet = "Diet weighting",
                         median_S_censor = "Censoring treatment")) |>
  ggplot(aes(x = family, y = Sensitivity, fill = Source)) +
  geom_col(position = position_dodge(width = .7), width = .6) +
  scale_fill_manual(values = c("Diet weighting" = "#275662",
                               "Censoring treatment" = "#D95F02")) +
  labs(x = "Contaminant family", y = "Median sensitivity (log2 fold range)",
       fill = "Methodological factor") +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(FIGURE_DIR, "Fig.S4_methodological_sensitivity_hierarchy_family.png"),
       p_family, width = 8, height = 5, dpi = 300)


# ---- 8F. Bivalve siphon factor: BMF (%) relative to 3.75% ----
p_siphon_delta <- bmf_siphon |>
  filter(is.finite(delta_pct)) |>
  ggplot(aes(x = 100 * factor_siphon, y = delta_pct,
             group = compound_name)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_vline(xintercept = 100 * REFERENCE_SIPHON,
             linetype = 3, color = "grey60") +
  geom_line(alpha = .7) + geom_point(size = 1.5) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Assumed bivalve siphon mass (%)",
       y = "Change in biomass-weighted BMF relative to 3.75% (%)") +
  theme_bw()
ggsave(file.path(FIGURE_DIR, "Fig.S2_siphon_sensitivity.png"),
       p_siphon_delta, width = 11, height = 6, dpi = 300)


# ---- 8G. Bivalve siphon factor: absolute BMFs ----
# Labels for the three HBCDD isomers
labels_hbcdd <- bmf_siphon |>
  filter(
    family == "HBCDD",
    is.finite(bmf),
    bmf > 0
  ) |>
  group_by(compound_name) |>
  slice_max(factor_siphon, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    label = case_when(
      str_detect(compound_name, "^a\\.") ~ "α-HBCDD",
      str_detect(compound_name, "^b\\.") ~ "β-HBCDD",
      str_detect(compound_name, "^g\\.") ~ "γ-HBCDD",
      TRUE ~ compound_name
    )
  )

# Main figure
p_siphon_bmf <- bmf_siphon |>
  filter(is.finite(bmf), bmf > 0) |>
  ggplot(aes(
    x = 100 * factor_siphon,
    y = bmf,
    group = compound_name
  )) +

  # Reference thresholds
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 100 * REFERENCE_SIPHON,
    linetype = "dotted",
    color = "grey40",
    linewidth = 0.6
  ) +

  # Compound trajectories
  geom_line(
    color = "grey30",
    alpha = 0.75,
    linewidth = 0.55
  ) +
  geom_point(
    color = "black",
    size = 1.6
  ) +

  # HBCDD labels at the right end of each trajectory
  ggrepel::geom_text_repel(
    data = labels_hbcdd,
    aes(
      x = 100 * factor_siphon,
      y = bmf,
      label = label
    ),
    inherit.aes = FALSE,
    direction = "y",
    hjust = 0,
    nudge_x = 0.15,
    size = 3.2,
    color = "black",
    segment.color = "grey60",
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +

  scale_y_log10() +

  scale_x_continuous(
    breaks = seq(2, 5, 0.5),
    limits = c(2, 5.9)
  ) +

  facet_wrap(
    ~ family,
    scales = "free_y",
    nrow = 1
  ) +

  labs(
    x = "Assumed bivalve siphon mass (%)",
    y = "Biomass-weighted BMF"
  ) +

  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "grey90",
      color = "grey40"
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

p_siphon_bmf

ggsave(file.path(FIGURE_DIR, "Fig.S3_siphon_sensitivity_bmf.png"),
       p_siphon_bmf, width = 11, height = 6, dpi = 300)


# ---- 8H. Sensitivity to numeric vs biomass diet weighting ----
p_numeric_biomass <- sensitivity_numeric_biomass |>
  mutate(
    # Convert log2 sensitivity into multiplicative difference
    fold_change = 2^S_numeric_biomass,

    # Order compounds within each family
    compound_plot = tidytext::reorder_within(
      compound_name,
      S_numeric_biomass,
      family
    )
  ) |>
  ggplot(aes(
    x = compound_plot,
    y = S_numeric_biomass,
    fill = family
  )) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(
    ~ family,
    scales = "free_y",
    ncol = 1
  ) +
  tidytext::scale_x_reordered() +
  geom_text(
    aes(label = paste0("x", round(fold_change, 2))),
    hjust = -0.15,
    size = 3
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  scale_fill_manual(
    values = c(
      PCB = "#275662",
      HBCDD = "#D95F02",
      PFAS = "#2A9D8F"
    )
  ) +
  labs(
    x = NULL,
    y = "Sensitivity to dietary weighting (log2 fold range)",
    title = "Sensitivity of BMFs to numerical versus biomass weighting",
    subtitle = "Median absolute difference across four censoring scenarios (gamma-MLE)"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

p_numeric_biomass

ggsave(
  "manuscript/figures/Fig.S6_numeric_biomass_sensitivity.png",
  p_numeric_biomass,
  width = 10,
  height = 11,
  dpi = 300
)

cat("\n03 completed |", n_distinct(bmf_compare$compound_name),
    "compounds;", sum(completeness$complete),
    "have complete 12-scenario grids.\n")
cat("Reference factor check: maximum abs difference = ",
    max(siphon_check$abs_diff), "\n", sep = "")
cat("Outputs: ", OUTPUT_DIR, " and ", FIGURE_DIR, "\n", sep = "")

