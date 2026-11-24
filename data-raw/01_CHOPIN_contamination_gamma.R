# ============================================================================
# CHOPIN | 01 - Chemical preparation, gamma-MLE scenarios, censoring audit
# Run first, from the CHOPIN project root (R >= 4.1; dplyr >= 1.1).
# Outputs: results/*.rds; manuscript/tables/*; manuscript/figures/*.
#
# SCIENTIFIC CONVENTIONS
# * Age-0 sole chemistry: explicitly configured in 00_CHOPIN_config.R.
#   Use the verified autumn G0 pools (set SOLE_SEASON="autumn" in config).
# * PCB/HBCDD: ng/g lipid weight; PFAS: ng/g dry weight.
# * PCB/PFAS: exclude when >=80% censored in either compartment.
# * HBCDD: retain all 3 isomers (including beta/gamma), flag heavy censoring.
#   HBCDD ratios dominated by substitutions must be interpreted exploratorily.
# * Four main scenarios: zero, LOQ/2, LOQ and GAMMA-MLE.
# * Gamma and lognormal are fitted ONCE here, using identical criteria,
#   and saved for 04 (so comparisons cannot drift from the main analysis).
# * For MLE, fit compound x biological group; >=3 samples and >=2 detected.
#   Fallback LOQ/2 only for groups needing imputation and unable to fit.
# ============================================================================

suppressPackageStartupMessages({
  library(chopin.data)
  library(tidyverse)
})

FIG_DIR <- "manuscript/figures"
OUT_DIR <- "manuscript/tables"
RESULT_DIR <- "results"
# Single source of truth for provenance choices and package-data export.
source("data-raw/00_CHOPIN_config.R", local = TRUE)

CENSOR_LIMIT <- 80
MIN_N <- 3L
MIN_UNCENSORED <- 2L
FORBIDDEN_PFAS <- c("MeFOSA", "EtFOSA") # Outside the biotic BMF analytical scope
for (d in c(FIG_DIR, OUT_DIR, RESULT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

for (nm in c("soles_contam", "benthos_contam",
             "soles_G0_pools_biometrics", "contam_info")) {
  data(list = nm, package = "chopin.data", envir = environment())
}

# Fail rather than silently mixing the current code with outdated outputs.
if (length(SOLE_SEASON) != 1L || is.na(SOLE_SEASON) ||
    !SOLE_SEASON %in% c("spring", "autumn")) {
  message("G0 source sole pools by recorded season:")
  print(soles_contam |> filter(class == "G0") |> count(season))
  stop("Invalid SOLE_SEASON in central configuration; choose spring or autumn.")
}

# 1. Chemical metadata and long-format concentrations ----------------------

chem <- contam_info |>
  filter(family %in% c("PCB", "HBCDD", "PFAS"),
         !(family == "PFAS" & chemical %in% FORBIDDEN_PFAS)) |>
  distinct(chemical, family, .keep_all = TRUE) |>
  transmute(chemical, family, loq_reference = if_else(family == "PCB", 0, LOQ_biote))

chem_cols <- chem |>
  filter(family != "PFAS") |>
  mutate(unit = "ng_glw") |>
  transmute(column = paste0(chemical, "_", unit), family, chemical) |>
  bind_rows(chem |>
    filter(family != "PCB") |>
    mutate(unit = if_else(family == "PFAS", "ng_gdw", "ng_glw")) |>
    transmute(column = paste0(chemical, "_", unit, "_censored"), family, chemical))

prepare_chem <- function(dat, compartment) {
  is_sole <- compartment == "Sole"
  dat <- if (is_sole) {
    dat |> filter(class == "G0" & season == SOLE_SEASON) |>
      rename(dry_fraction = dry_content_percent, lipid_fraction = lip_dw_percent) |>
      mutate(taxa = "Sole", feeding_mode = NA_character_)
  } else {
    dat |> rename(dry_fraction = mass_dry_percent, lipid_fraction = lip_dw_percent)
  }

  dat |>
    select(sample_TAG, species, taxa, feeding_mode, dry_fraction, lipid_fraction,
           all_of(chem_cols$column)) |>
    pivot_longer(all_of(chem_cols$column), names_to = "column", values_to = "reading") |>
    left_join(chem_cols, by = c("column" = "column"), relationship = "many-to-one") |>
    mutate(
      sample_id = sample_TAG,
      compartment = compartment,
      group = if_else(compartment == "Sole", "Solea solea", as.character(species)),
      unit = if_else(family == "PFAS", "ng_gdw", "ng_glw"),
      version = if_else(str_ends(column, "_censored"), "flag", "original")
    ) |>
    select(sample_id, compartment, group, species, taxa, feeding_mode,
           dry_fraction, lipid_fraction, family, compound_name = chemical,
           unit, version, reading) |>
    pivot_wider(names_from = version, values_from = reading) |>
    # HBCDD: 'flag' == 0 is censored, original value is the LOQ.
    # PFAS: 'flag' == 0 is censored, use the reference LOQ.
    # PCB: all observations uncensored.
    left_join(chem, by = c("compound_name" = "chemical", "family" = "family"),
              relationship = "many-to-one") |>
    mutate(
      censored = case_when(
        family == "PCB" ~ FALSE,
        is.na(flag) ~ NA,
        family == "HBCDD" ~ flag == 0,
        family == "PFAS" ~ flag == 0
      ),
      concentration = case_when(
        family == "PCB" ~ original,
        family == "HBCDD" ~ original,
        family == "PFAS" & censored %in% TRUE ~ loq_reference,
        family == "PFAS" ~ flag
      ),
      loq = case_when(
        family == "HBCDD" & censored %in% TRUE ~ original,
        family == "HBCDD" & censored %in% FALSE ~ 0,
        family == "PFAS" ~ loq_reference,
        family == "PCB" ~ 0
      ),
      compound = paste0(compound_name, "_", unit)
    ) |>
    # As in the previous pipeline, missing analytical concentrations are not used.
    filter(!is.na(concentration)) |>
    select(sample_id, compartment, group, species, taxa, feeding_mode,
           dry_fraction, lipid_fraction, family, compound, compound_name,
           unit, concentration, censored, loq) |>
    arrange(sample_id, family, compound_name)
}

# Auditable pool-level provenance for the selected chemical subset.
sole_pool_provenance <- soles_contam |>
  filter(class == "G0", season == SOLE_SEASON) |>
  distinct(sample_TAG, date, season, class) |>
  arrange(sample_TAG)
readr::write_csv(sole_pool_provenance,
                 file.path(OUT_DIR, "sole_pool_provenance_current_run.csv"))
sole <- prepare_chem(soles_contam, "Sole")
prey <- prepare_chem(benthos_contam, "Prey")
if (!nrow(sole) || !nrow(prey) ||
    !all(c("PCB", "HBCDD", "PFAS") %in% unique(sole$family))) {
  stop("No valid source chemical data in one or more compartments/families for SOLE_SEASON=", SOLE_SEASON)
}

audit_data <- bind_rows(sole, prey)

soles_G0_pools_stats <- soles_G0_pools_biometrics |>
  group_by(sample_TAG, season) |>
  summarise(n = n(), mean_length_TL_cm = mean(length_TL_cm),
            sd_length_TL_cm = sd(length_TL_cm),
            mean_mass_tot_gww = mean(mass_tot_gww),
            sd_mass_tot_gww = sd(mass_tot_gww), .groups = "drop")


# 2. Quality checks and family-specific compound exclusion ----------------------------

censor_audit <- audit_data |>
  group_by(compartment, family, compound_name) |>
  summarise(n_total = n(), n_censored = sum(censored),
            pct_censored = 100 * n_censored / n_total, .groups = "drop")

# The 80% exclusion applies to PCB and PFAS only. The three HBCDD isomers
# are retained for the biologically motivated, explicitly exploratory analysis.
# Compounds absent from either compartment remain excluded for all families.
compound_decisions <- censor_audit |>
  group_by(family, compound_name) |>
  summarise(
    n_compartments = n_distinct(compartment),
    high_censoring = any(pct_censored >= CENSOR_LIMIT),
    all_censored_in_compartment = any(n_censored == n_total),
    censored_in_sole_pct = pct_censored[match("Sole", compartment)],
    censored_in_prey_pct = pct_censored[match("Prey", compartment)],
    .groups = "drop"
  ) |>
  mutate(
    excluded = n_compartments != 2L |
      (family != "HBCDD" & high_censoring),
    exclusion_reason = case_when(
      n_compartments != 2L ~ "Missing compartment",
      family == "HBCDD" & high_censoring ~ "Retained: HBCDD exception to 80% rule",
      high_censoring & censored_in_sole_pct >= CENSOR_LIMIT &
        censored_in_prey_pct >= CENSOR_LIMIT ~ "Both compartments",
      high_censoring & censored_in_sole_pct >= CENSOR_LIMIT ~ "Sole",
      high_censoring & censored_in_prey_pct >= CENSOR_LIMIT ~ "Prey",
      TRUE ~ "Retained"
    ),
    evidence_tier = case_when(
      excluded ~ "Excluded",
      family == "HBCDD" & all_censored_in_compartment ~
        "Exploratory: fully censored compartment",
      family == "HBCDD" & high_censoring ~
        "Exploratory: highly censored compartment",
      TRUE ~ "Standard"
    )
  )
selected_compound_keys <- compound_decisions |>
  filter(!excluded) |>
  select(family, compound_name)
selected_compounds <- selected_compound_keys |> pull(compound_name)

# Stop early if any HBCDD isomer has disappeared before fitting and plotting.
expected_hbcdd <- c("a.HBCDD", "b.HBCDD", "g.HBCDD")
missing_hbcdd <- setdiff(expected_hbcdd,
                         selected_compound_keys$compound_name[
                           selected_compound_keys$family == "HBCDD"])
if (length(missing_hbcdd)) {
  stop("Missing HBCDD isomer(s) after selection: ",
       paste(missing_hbcdd, collapse = ", "),
       ". Inspect audit_data and source chemistry.")
}


# 3. Gamma-MLE and lognormal sensitivity fits -------------------------------
# Each censored observation contributes log[P(X <= its own LOQ)] to the
# likelihood; each uncensored observation contributes log[density(X)].
# Parameters are fitted on a rescaled concentration axis for stability.
# Gamma: shape/scale; lognormal: meanlog/sdlog.  Bounds and multi-start
# optimisation are identical to the previously supplied comparison script 04.
#
# The fitted value is E[X | X <= LOQ], NOT the unconditional model mean.
# Source observations and the three substitution scenarios are unchanged.

FIT_KEYS <- c("family", "compound_name", "compartment", "group")

fit_censored_group <- function(df, distribution) {
  stopifnot(distribution %in% c("gamma", "lognormal"), nrow(df) > 0)
  cc <- which(df$censored %in% TRUE)
  uu <- which(df$censored %in% FALSE)
  n <- nrow(df)
  status <- "not_needed"
  fit_error <- NA_character_
  p1 <- p2 <- NA_real_
  value <- df$concentration
  value[cc] <- df$loq[cc] / 2

  if (length(cc)) {
    if (any(!is.finite(df$concentration[uu]) | df$concentration[uu] <= 0)) {
      status <- "failed"
      fit_error <- "Nonpositive/invalid uncensored concentrations"
    } else if (n < MIN_N || length(uu) < MIN_UNCENSORED) {
      status <- "insufficient_data"
    } else {
      unit_scale <- median(df$concentration[uu])
      y <- df$concentration[uu] / unit_scale
      L <- df$loq[cc] / unit_scale

      starts <- if (distribution == "lognormal") {
        mu0 <- mean(log(y)); s0 <- max(sd(log(y)), 0.3)
        list(c(mu0, log(s0)), c(mu0 - 1, log(s0 + 0.5)),
             c(mu0 + 1, log(s0 + 0.5)))
      } else {
        m <- mean(y); v <- max(var(y), m^2 / 4)
        k <- max(0.1, m^2 / v)
        list(log(c(k, m / k)), log(c(1, m)), log(c(2, m / 2)))
      }

      neg_loglik <- function(par) {
        if (distribution == "lognormal") {
          a <- dlnorm(y, meanlog = par[1], sdlog = exp(par[2]), log = TRUE)
          b <- plnorm(L, meanlog = par[1], sdlog = exp(par[2]), log.p = TRUE)
        } else {
          a <- dgamma(y, shape = exp(par[1]), scale = exp(par[2]), log = TRUE)
          b <- pgamma(L, shape = exp(par[1]), scale = exp(par[2]), log.p = TRUE)
        }
        if (any(!is.finite(c(a, b)))) return(1e100)
        -sum(a) - sum(b)
      }
      trials <- lapply(starts, function(start) {
        tryCatch(suppressWarnings(optim(
          start, neg_loglik, method = "L-BFGS-B",
          lower = c(-20, -15), upper = c(20, 15),
          control = list(maxit = 2000)
        )), error = function(e) NULL)
      })
      trials <- Filter(function(x) !is.null(x) && x$convergence == 0 &&
                         is.finite(x$value) && x$value < 1e99, trials)

      if (!length(trials)) {
        status <- "failed"
        fit_error <- "Likelihood optimisation failed"
      } else {
        best <- trials[[which.min(vapply(trials, `[[`, numeric(1), "value"))]]
        par <- best$par
        if (distribution == "lognormal") {
          mu <- par[1]; sig <- exp(par[2]); z <- (log(L) - mu) / sig
          estimates <- exp(mu + sig^2 / 2 +
                             pnorm(z - sig, log.p = TRUE) -
                             pnorm(z, log.p = TRUE)) * unit_scale
          p1 <- mu + log(unit_scale)  # meanlog on original concentration axis
          p2 <- sig                    # sdlog
        } else {
          shape <- exp(par[1]); scale <- exp(par[2]) * unit_scale
          # E[X|X <= LOQ] = shape*scale*F_gamma(shape+1)/F_gamma(shape)
          estimates <- shape * scale * exp(
            pgamma(df$loq[cc], shape = shape + 1, scale = scale, log.p = TRUE) -
              pgamma(df$loq[cc], shape = shape, scale = scale, log.p = TRUE))
          p1 <- shape; p2 <- scale
        }
        if (all(is.finite(estimates) & estimates >= 0 &
                estimates <= df$loq[cc] * (1 + 1e-8))) {
          value[cc] <- estimates
          status <- "fitted"
        } else {
          status <- "failed"
          fit_error <- "Invalid conditional expectation"
        }
      }
    }
  }

  df$distribution <- distribution
  df$value <- value
  df$imputation_method <- if_else(
    df$censored,
    if (status == "fitted") "mle" else "loq_half_fallback",
    "observed"
  )
  audit <- df |> slice(1) |>
    transmute(across(all_of(FIT_KEYS)), distribution, n = n,
              n_censored = length(cc), n_uncensored = length(uu),
              pct_censored = 100 * length(cc) / n,
              fit_status = status,
              fallback_method = if (status %in% c("failed", "insufficient_data"))
                "loq_half" else NA_character_,
              fit_error = fit_error,
              parameter_1 = p1, parameter_2 = p2,
              n_imputed = length(cc),
              n_imputed_mle = if (status == "fitted") length(cc) else 0L,
              n_imputed_fallback = if (status %in% c("failed", "insufficient_data"))
                length(cc) else 0L,
              n_unresolved = sum(!is.finite(df$value)),
              imputed_min = if (length(cc)) min(value[cc]) else NA_real_,
              imputed_median = if (length(cc)) median(value[cc]) else NA_real_,
              imputed_max = if (length(cc)) max(value[cc]) else NA_real_)
  list(data = df, audit = audit)
}

retained <- audit_data |>
  semi_join(selected_compound_keys, by = c("family", "compound_name"))

# Exactly the same analytical records and eligibility criteria for both fits.
groups <- retained |> group_by(across(all_of(FIT_KEYS))) |> group_split()
model_fits <- purrr::map(c("gamma", "lognormal"), function(dist) {
  purrr::map(groups, fit_censored_group, distribution = dist)
}) |> purrr::flatten()

imputed_data <- purrr::map_dfr(model_fits, "data")
fit_comparison <- purrr::map_dfr(model_fits, "audit")
fit_audit <- fit_comparison |> filter(distribution == "gamma")

gamma_data <- imputed_data |> filter(distribution == "gamma")
scenarios <- gamma_data |>
  mutate(value_zero = if_else(censored, 0, concentration),
         value_loq_half = if_else(censored, loq / 2, concentration),
         value_loq = if_else(censored, loq, concentration),
         value_mle = value) |>
  select(-value, -distribution)
sole_scenarios <- scenarios |> filter(compartment == "Sole")
prey_scenarios <- scenarios |> filter(compartment == "Prey")
data_imputed <- scenarios

stopifnot(nrow(scenarios) == nrow(retained),
          !anyNA(scenarios$value_mle),
          all(scenarios$value_mle[!scenarios$censored] ==
                scenarios$concentration[!scenarios$censored]))

# 4. Full chemical audit, compact supplementary audit and summary -----------

censoring_audit <- audit_data |>
  group_by(family, compound_name, compartment, group) |>
  summarise(n = n(), n_censored = sum(censored),
            n_uncensored = sum(!censored), pct_censored = 100 * n_censored / n,
            loq_min = if (any(censored)) min(loq[censored]) else NA_real_,
            loq_max = if (any(censored)) max(loq[censored]) else NA_real_,
            .groups = "drop") |>
  left_join(compound_decisions, by = c("family", "compound_name"))

audit_final <- censoring_audit |>
  left_join(fit_audit |>
              select(-n, -n_censored, -n_uncensored, -pct_censored),
            by = FIT_KEYS, relationship = "one-to-one") |>
  mutate(mle_distribution = if_else(excluded, NA_character_, "gamma"),
         fit_status = if_else(excluded, "excluded", fit_status),
         gamma_shape = if_else(fit_status == "fitted", parameter_1, NA_real_),
         gamma_scale = if_else(fit_status == "fitted", parameter_2, NA_real_),
         unit = if_else(family == "PFAS", "ng/g dry weight", "ng/g lipid weight"),
         fallback_reason = case_when(
           coalesce(n_imputed_fallback, 0L) == 0 ~ NA_character_,
           n < MIN_N ~ "Fewer than 3 observations",
           n_uncensored < MIN_UNCENSORED ~ "Fewer than 2 uncensored observations",
           fit_status == "failed" ~ "MLE fitting failed",
           TRUE ~ "Other"
         )) |>
  select(-parameter_1, -parameter_2) |>
  arrange(family, compound_name, compartment, group)

fallback_summary <- fit_audit |>
  summarise(n_total = n(),
            n_requiring_imputation = sum(n_censored > 0),
            n_fitted = sum(fit_status == "fitted"),
            n_fallback = sum(n_imputed_fallback > 0)) |>
  mutate(fallback_pct_total = 100 * n_fallback / n_total,
         fallback_pct_censored = 100 * n_fallback / n_requiring_imputation)
fallback_reasons <- audit_final |>
  filter(!is.na(fallback_reason)) |>
  count(fallback_reason, name = "n")
audit_summary <- tibble(
  name = c(names(fallback_summary), fallback_reasons$fallback_reason),
  value = c(as.numeric(unlist(fallback_summary[1, ], use.names = FALSE)),
            fallback_reasons$n)
)

# Exactly one row per compound, including HBCDD exception decisions.
audit_supp <- censor_audit |>
  select(family, compound_name, compartment, n_total, pct_censored) |>
  pivot_wider(names_from = compartment, values_from = c(n_total, pct_censored)) |>
  left_join(compound_decisions, by = c("family", "compound_name"),
            relationship = "one-to-one") |>
  transmute(Family = family, Compound = compound_name,
            `Sole n` = n_total_Sole,
            `Sole censored (%)` = round(pct_censored_Sole, 1),
            `Prey n` = n_total_Prey,
            `Prey censored (%)` = round(pct_censored_Prey, 1),
            Excluded = if_else(excluded, "Yes", "No")) |>
  arrange(Family, Compound)

write_csv(audit_final, file.path(OUT_DIR, "Table_Sx_Censoring_MLE_audit.csv"), na = "")
write_csv(audit_supp, file.path(OUT_DIR, "Table_Sx_Censoring_summary.csv"), na = "")
write_csv(censor_audit, file.path(OUT_DIR, "censoring_by_compartment.csv"), na = "")
write_csv(compound_decisions, file.path(OUT_DIR, "compound_exclusion_decisions.csv"), na = "")
write_csv(audit_summary, file.path(OUT_DIR, "censoring_audit_summary.csv"), na = "")
write_csv(fit_comparison, file.path(OUT_DIR, "gamma_lognormal_fit_audit.csv"), na = "")

# Complete audit in Excel; the CSV remains available without writexl.
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(`Printable Table Sx` = audit_supp,
         `Full gamma MLE audit` = audit_final,
         `Gamma and lognormal fits` = fit_comparison,
         `Exclusion decisions` = compound_decisions,
         `Audit summary` = audit_summary),
    file.path(OUT_DIR, "Supplementary_Censoring_MLE_audit.xlsx"))
} else {
  warning("Install 'writexl' to generate the Excel audit; CSV exports are saved.")
}

# This is the authoritative handoff to scripts 02, 03 and 04.
for (nm in c("sole", "prey", "audit_data", "retained", "scenarios",
             "imputed_data", "fit_comparison",
             "fit_audit", "audit_final", "audit_supp", "audit_summary",
             "censor_audit", "censoring_audit", "compound_decisions",
             "selected_compound_keys", "soles_G0_pools_stats")) {
  saveRDS(get(nm), file.path(RESULT_DIR, paste0(nm, ".rds")))
}
if (SAVE_TO_PACKAGE) {
  usethis::use_data(sole, prey, scenarios, fit_audit, audit_final,
                    audit_supp, audit_summary, compound_decisions,
                    overwrite = TRUE)
}

# 5. Supplementary contamination figure (LOQ/2) ---------------------------

# Use the analysed compounds; all values are within-family comparable units.
contam_all <- retained |>
  mutate(value = if_else(censored, loq / 2, concentration),
         species_label = str_replace_all(as.character(species), "_", " "),
         plot_group = case_when(
           compartment == "Sole" ~ "Sole",
           taxa == "Bivalves" ~ "Bivalves",
           taxa == "Polychaetes" ~ "Polychaetes",
           taxa == "Crustaceans" ~ "Crustaceans",
           TRUE ~ as.character(taxa)))

n_expected <- selected_compound_keys |> count(family, name = "n_expected")
contam_sample <- contam_all |>
  group_by(sample_id, species_label, plot_group, family) |>
  summarise(n_valid = sum(is.finite(value)), total = sum(value, na.rm = TRUE),
            .groups = "drop") |>
  left_join(n_expected, by = "family") |>
  filter(n_valid == n_expected, is.finite(total), total > 0)

contam_summary <- contam_sample |>
  group_by(family, plot_group, species_label) |>
  summarise(n = n(), mean = mean(total), q25 = quantile(total, .25),
            q75 = quantile(total, .75), .groups = "drop")
sp_order <- c(sort(setdiff(unique(contam_sample$species_label), "Solea solea")),
              "Solea solea")
contam_sample <- contam_sample |>
  mutate(species_label = factor(species_label, levels = sp_order))
contam_summary <- contam_summary |>
  mutate(species_label = factor(species_label, levels = sp_order))

p_contamination <- ggplot(contam_sample,
                          aes(total, species_label, color = plot_group)) +
  geom_point(position = position_jitter(width = 0, height = .12, seed = 123),
             alpha = .4, size = 1.7) +
  geom_segment(data = contam_summary,
               aes(x = q25, xend = q75, y = species_label, yend = species_label),
               inherit.aes = FALSE, linewidth = .9) +
  geom_point(data = contam_summary, aes(x = mean, y = species_label),
             inherit.aes = FALSE, shape = 18, size = 3) +
  facet_wrap(~family, scales = "free_x", nrow = 1,
             labeller = labeller(family = c(
               PCB = "PCBs (ng/g lipid)", HBCDD = "HBCDDs (ng/g lipid)",
               PFAS = "PFASs (ng/g dry)"))) +
  scale_x_log10() +
  scale_color_manual(values = c(Sole = "#222222", Bivalves = "#D95F02",
                                Polychaetes = "#275662", Crustaceans = "#2A9D8F")) +
  labs(x = "Sum of compound concentrations (log scale)", y = NULL,
       color = "Organism group") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 9),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "S_contamination_levels.png"), p_contamination,
       width = 13, height = 8, dpi = 300)


# 6. Supplementary HBCDD isomer profile (LOQ/2) ----------------------------

# Work in long format to avoid missing columns after pivot_wider().
# A complete profile requires exactly one finite value per isomer and pool.
hbcdd_data <- contam_all |>
  filter(family == "HBCDD", compound_name %in% expected_hbcdd) |>
  mutate(isomer = recode(compound_name,
                         "a.HBCDD" = "alpha", "b.HBCDD" = "beta",
                         "g.HBCDD" = "gamma")) |>
  select(sample_id, species, plot_group, isomer, value)

hbcdd_complete <- hbcdd_data |>
  group_by(sample_id, species, plot_group) |>
  filter(n() == 3L, n_distinct(isomer) == 3L,
         all(is.finite(value)), sum(value) > 0) |>
  ungroup()

if (!nrow(hbcdd_complete)) {
  stop("No complete HBCDD profiles for the supplementary figure. ",
       "Inspect hbcdd_data and censored concentrations.")
}

hbcdd_profiles <- hbcdd_complete |>
  group_by(plot_group, species, isomer) |>
  summarise(mean_conc = mean(value), n_samples = n_distinct(sample_id),
            .groups = "drop") |>
  group_by(plot_group, species) |>
  mutate(percentage = 100 * mean_conc / sum(mean_conc)) |>
  ungroup() |>
  mutate(species_label = paste0(str_replace_all(as.character(species), "_", " "),
                                " (n=", n_samples, ")"),
         plot_group = factor(plot_group,
                             levels = c("Sole", "Polychaetes", "Crustaceans", "Bivalves")),
         isomer = factor(isomer, levels = c("alpha", "beta", "gamma")))
sp_order_hbcdd <- hbcdd_profiles |>
  distinct(plot_group, species, species_label) |>
  arrange(plot_group, species) |> pull(species_label)
hbcdd_profiles <- hbcdd_profiles |>
  mutate(species_label = factor(species_label, levels = rev(sp_order_hbcdd)))

p_hbcdd_profiles <- ggplot(hbcdd_profiles,
                           aes(x = percentage, y = species_label, fill = isomer)) +
  geom_col(orientation = "y", position = position_stack(reverse = TRUE),
           width = .75, color = "white", linewidth = .2) +
  facet_grid(plot_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_fill_manual(values = c(alpha = "#275662", beta = "#E9C46A", gamma = "#D95F02"),
                    labels = c(alpha = "α-HBCDD", beta = "β-HBCDD", gamma = "γ-HBCDD")) +
  labs(x = "Relative isomer contribution (%)", y = NULL, fill = "HBCDD isomer") +
  theme_bw(base_size = 11) +
  theme(strip.placement = "outside", strip.background = element_rect(fill = "grey95"),
        strip.text.y.left = element_text(angle = 90, face = "bold"),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8.5, face = "italic"),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "S_HBCDD_isomer_profiles.png"), p_hbcdd_profiles,
       width = 11, height = 9, dpi = 300)


# 7. Final checks ----------------------------------------------------------

print(compound_decisions |> count(family, excluded, exclusion_reason))
print(compound_decisions |>
        filter(family == "HBCDD") |>
        select(compound_name, high_censoring, all_censored_in_compartment,
               excluded, evidence_tier))
print(fit_audit |> count(family, fit_status, fallback_method))
print(audit_final |> summarise(rows = n(), excluded = sum(excluded),
                                unresolved = sum(n_unresolved, na.rm = TRUE)))
print(audit_summary)
cat("01 completed | main MLE distribution: gamma; retained compounds:",
    nrow(selected_compound_keys), "\n")
# Objects available: sole, prey, censor_audit, compound_decisions,
# selected_compounds, selected_compound_keys, sole_scenarios, prey_scenarios, scenarios,
# fit_audit, audit_final, contam_all, contam_sample, contam_summary,
# hbcdd_profiles, p_contamination, p_hbcdd_profiles.
