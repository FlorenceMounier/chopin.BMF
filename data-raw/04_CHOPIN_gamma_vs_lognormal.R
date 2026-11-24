# ============================================================================
# CHOPIN | 04 - Sensitivity to MLE distribution: gamma vs lognormal
# Run AFTER 01, 02 and 03.  No refitting here: the 2 fits are computed ONCE
# with identical settings in 01 and passed to 04 via results/imputed_data.rds.
#
# Analyses from the former script 04 are all preserved:
#   * Fit-status counts and paired fit outcomes (494 groups in prior run)
#   * Individual censored concentrations, relative differences and summaries
#   * Predator / measured-prey and diet-based BMFs (equal/numeric/biomass)
#   * BMF paired differences, maxima and crossings of the unity threshold
#   * Two supplementary scatterplots and complete CSV/RDS exports
# NEW: a strict numeric check that gamma BMFs equal the main script 03 BMFs.
# ============================================================================

suppressPackageStartupMessages(library(tidyverse))
OUTPUT_DIR <- "results"
OUT <- file.path(OUTPUT_DIR, "gamma_vs_lognormal")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIT_KEYS <- c("family", "compound_name", "compartment", "group")
BMF_KEYS <- c("scenario", "family", "compound_name", "unit")

needed <- file.path(OUTPUT_DIR, c(
  "imputed_data.rds", "fit_comparison.rds", "prop_diet.rds", "bmf_compare.rds"
))
if (any(!file.exists(needed))) {
  stop("Run scripts 01-03 first; missing: ", paste(needed[!file.exists(needed)], collapse = ", "))
}
imp_data <- readRDS(needed[1])
fit_comparison <- readRDS(needed[2])
prop_diet <- readRDS(needed[3])
bmf_reference <- readRDS(needed[4])
stopifnot(setequal(unique(imp_data$distribution), c("gamma", "lognormal")))

# ---- 1. Compare model feasibility and concentrations ----------------------
fit_summary <- fit_comparison |>
  count(distribution, fit_status, name = "n_groups")
fit_pairs <- fit_comparison |>
  select(all_of(FIT_KEYS), distribution, fit_status) |>
  pivot_wider(names_from = distribution, values_from = fit_status,
              names_prefix = "status_") |>
  count(status_lognormal, status_gamma, name = "n_groups")

# One pair per censored analytical record. Fallback records are explicitly
# identified, so they can be distinguished from successful MLE imputations.
status_pairs <- fit_comparison |>
  select(all_of(FIT_KEYS), distribution, fit_status) |>
  pivot_wider(names_from = distribution, values_from = fit_status,
              names_prefix = "status_")
imputation_compare <- imp_data |>
  filter(censored) |>
  select(all_of(FIT_KEYS), sample_id, loq, distribution, value) |>
  pivot_wider(names_from = distribution, values_from = value) |>
  left_join(status_pairs, by = FIT_KEYS, relationship = "many-to-one") |>
  mutate(both_fitted = status_gamma == "fitted" &
           status_lognormal == "fitted",
         delta_pct = if_else(lognormal > 0,
                             100 * (gamma / lognormal - 1), NA_real_),
         abs_delta_pct = abs(delta_pct),
         log2_ratio = if_else(gamma > 0 & lognormal > 0,
                              log2(gamma / lognormal), NA_real_))

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

imputation_summary <- imputation_compare |>
  group_by(family, compound_name) |>
  summarise(n_censored = n(),
            n_fitted_observations = sum(both_fitted),
            median_delta_pct = safe_median(delta_pct),
            max_abs_delta_pct = safe_max(abs_delta_pct),
            .groups = "drop")
# Among genuinely fitted pairs only (not diluted by identical LOQ/2 fallbacks):
imputation_summary_fitted <- imputation_compare |>
  filter(both_fitted) |>
  group_by(family, compound_name) |>
  summarise(n_censored = n(),
            median_abs_delta_pct = safe_median(abs_delta_pct),
            max_abs_delta_pct = safe_max(abs_delta_pct),
            .groups = "drop")

# ---- 2. Recalculate BMFs identically to script 03 -------------------------
# Species means -> phylum means (equal weight for measured species) ->
# across-phylum arithmetic mean, numeric weights or biomass weights.
# The dietary denominator is a weighted mean of NORMALISED concentrations,
# not a physical concentration in a reconstructed whole-prey tissue mixture.

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}
safe_ratio <- function(num, den) {
  if_else(is.finite(num) & is.finite(den) & den > 0,
          num / den, NA_real_)
}
check_unique <- function(dat, by) {
  dup <- dat |> count(across(all_of(by))) |> filter(n > 1)
  if (nrow(dup)) stop("Duplicate key: ", paste(by, collapse = ", "))
  invisible(TRUE)
}

diet_data <- prop_diet |>
  transmute(Phylum, p_numeric = p_numeric_phylum,
            p_biomass = p_biomass_phylum)
taxa_mapping <- tribble(
  ~taxa,         ~Phylum,
  "Polychaetes", "Annelida",
  "Crustaceans",  "Arthropoda",
  "Bivalves",    "Mollusca"
)
chem_data <- imp_data |>
  mutate(scenario = distribution,
         value = if_else(is.finite(value), value, NA_real_))
sole_data <- chem_data |> filter(compartment == "Sole")
prey_data <- chem_data |>
  filter(compartment == "Prey") |>
  left_join(taxa_mapping, by = "taxa", relationship = "many-to-one")
if (anyNA(prey_data$Phylum)) stop("Unmapped prey taxa in 04")

# Same function body and handling of missing concentrations as in 03.
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

results <- calc_bmfs(sole_data, prey_data, diet_data)
prey_bmf <- results$prey
phylum_bmf <- results$phylum
diet_bmf <- results$comparison |>
  rename(distribution = scenario)

# Gamma in script 04 MUST match the published/main pipeline in script 03.
join_keys <- c("family", "compound_name", "unit", "diet_method")
gamma_keys <- diet_bmf |> filter(distribution == "gamma") |> select(all_of(join_keys))
main_keys <- bmf_reference |> filter(scenario == "mle") |> select(all_of(join_keys))
if (nrow(anti_join(gamma_keys, main_keys, by = join_keys)) ||
    nrow(anti_join(main_keys, gamma_keys, by = join_keys))) {
  stop("Gamma/main BMF compound-method keys differ (scripts 03 and 04)")
}
reference_check <- diet_bmf |>
  filter(distribution == "gamma") |>
  select(family, compound_name, unit, diet_method, bmf_gamma = bmf) |>
  full_join(bmf_reference |>
              filter(scenario == "mle") |>
              select(family, compound_name, unit, diet_method,
                     bmf_main = bmf),
            by = c("family", "compound_name", "unit", "diet_method"),
            relationship = "one-to-one") |>
  mutate(abs_diff = abs(bmf_gamma - bmf_main))
if (any(is.na(reference_check$bmf_gamma) != is.na(reference_check$bmf_main)) ||
    any(reference_check$abs_diff >
          1e-8 * pmax(1, abs(reference_check$bmf_main)), na.rm = TRUE)) {
  write_csv(reference_check, file.path(OUT, "gamma_main_mismatch.csv"))
  stop("Gamma BMFs from 04 do not reproduce script 03; see gamma_main_mismatch.csv")
}

# ---- 3. Paired BMF differences and threshold sensitivity ------------------
bmf_comparison <- diet_bmf |>
  select(family, compound_name, unit, diet_method,
         distribution, bmf, c_sole, c_diet, n_phyla) |>
  pivot_wider(names_from = distribution,
              values_from = c(bmf, c_sole, c_diet, n_phyla)) |>
  mutate(delta_pct = if_else(bmf_lognormal > 0,
                             100 * (bmf_gamma / bmf_lognormal - 1), NA_real_),
         abs_delta_pct = abs(delta_pct),
         log2_ratio = if_else(bmf_lognormal > 0 & bmf_gamma > 0,
                              log2(bmf_gamma / bmf_lognormal), NA_real_),
         threshold_change = if_else(is.finite(bmf_lognormal) & is.finite(bmf_gamma),
                                    (bmf_lognormal > 1) != (bmf_gamma > 1), NA))

bmf_summary <- bmf_comparison |>
  group_by(family, diet_method) |>
  summarise(n_bmf = sum(is.finite(delta_pct)),
            median_abs_delta_pct = safe_median(abs_delta_pct),
            max_abs_delta_pct = safe_max(abs_delta_pct),
            n_threshold_changes = sum(threshold_change, na.rm = TRUE),
            .groups = "drop")

# ---- 4. Save ALL comparison data and diagnostic figures -------------------
export_tables <- list(
  fit_summary = fit_summary, fit_pairs = fit_pairs,
  fit_comparison = fit_comparison,
  imputation_compare = imputation_compare,
  imputation_summary = imputation_summary,
  imputation_summary_fitted = imputation_summary_fitted,
  prey_bmf = prey_bmf, phylum_bmf = phylum_bmf,
  diet_bmf = diet_bmf, bmf_comparison = bmf_comparison,
  bmf_summary = bmf_summary, reference_check = reference_check
)
purrr::iwalk(export_tables, function(tbl, name) {
  write_csv(tbl, file.path(OUT, paste0(name, ".csv")), na = "")
  saveRDS(tbl, file.path(OUT, paste0(name, ".rds")))
})
saveRDS(list(fit_audit = fit_comparison, imputed = imp_data,
             diet_bmf = diet_bmf, bmf_comparison = bmf_comparison),
        file.path(OUT, "gamma_lognormal_full_results.rds"))

p_imputation <- imputation_compare |>
  ggplot(aes(lognormal, gamma, color = family)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~family, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(x = "Lognormal MLE (or LOQ/2 fallback)",
       y = "Gamma MLE (or LOQ/2 fallback)",
       title = "Individual left-censored concentrations")

p_bmf <- bmf_comparison |>
  ggplot(aes(bmf_lognormal, bmf_gamma, color = family)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(alpha = 0.8) +
  facet_wrap(~diet_method, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(x = "BMF: lognormal MLE (with fallback)",
       y = "BMF: gamma MLE (with fallback)",
       title = "Effect of MLE distribution on diet-based BMFs")

ggsave(file.path(OUT, "imputation_gamma_vs_lognormal.png"),
       p_imputation, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "bmf_gamma_vs_lognormal.png"),
       p_bmf, width = 11, height = 5.5, dpi = 300)

print(fit_summary)
print(fit_pairs)
print(bmf_summary)
cat("04 completed | gamma main consistency check passed; exports:", OUT, "\n")
