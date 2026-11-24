# CHOPIN | 09 - HBCDD isomer-specific exploratory analysis
# Run AFTER scripts 01 (chemistry) and 02 (diet); does not use precomputed BMFs.
# All three HBCDD are included in this targeted analysis independently of 80% censoring.
# The main script 01 should use the documented HBCDD exception to permit QA matching.
# Beta/gamma are evaluated separately, including analytically defensible bounds.
# Units: ng/g total lipid for all HBCDDs.
# No values are inferred below LOQ beyond the explicit sensitivity assumptions.

suppressPackageStartupMessages(library(tidyverse))
DIR <- "results/hbcdd_isomers"
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists("results/audit_data.rds"),
          file.exists("results/prop_diet.rds"))
MIN_N <- 3L
MIN_UNCENSORED <- 2L
FIT_KEYS <- c("family", "compound_name", "compartment", "group")
ISOMERS <- c("a.HBCDD", "b.HBCDD", "g.HBCDD")

chem <- readRDS("results/audit_data.rds") |>
  filter(family == "HBCDD", compound_name %in% ISOMERS) |>
  mutate(censored = as.logical(censored),
         concentration = as.numeric(concentration), loq = as.numeric(loq))
stopifnot(all(c("sample_id", "species", "taxa", "unit", FIT_KEYS,
                "concentration", "censored", "loq") %in% names(chem)),
          !anyDuplicated(chem[c("compartment", "sample_id", "compound_name")]),
          !anyNA(chem$censored),
          all(is.finite(chem$concentration)),
          all(!chem$censored | (is.finite(chem$loq) & chem$loq > 0)),
          all(chem$unit == "ng_glw"),
          setequal(unique(chem$compound_name), ISOMERS),
          setequal(unique(chem$compartment), c("Sole", "Prey")))

weights <- readRDS("results/prop_diet.rds") |>
  transmute(Phylum, p_numeric = p_numeric_phylum,
            p_biomass = p_biomass_phylum)
stopifnot(!anyDuplicated(weights["Phylum"]))
taxa_mapping <- tibble(taxa = c("Polychaetes", "Crustaceans", "Bivalves"),
                      Phylum = c("Annelida", "Arthropoda", "Mollusca"))

# Reuse the SAME gamma model specification as script 01 (duplicate here so
# this analysis also works from a clean R session with only saved source data).
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


# Fit each isomer x species separately, without excluding high-censoring groups.
fits <- chem |>
  group_by(across(all_of(FIT_KEYS))) |>
  group_split() |>
  purrr::map(fit_censored_group, distribution = "gamma")
fit_detail <- purrr::map_dfr(fits, "audit")
chem_fit <- purrr::map_dfr(fits, "data")

# Four nominal scenarios plus true lower and upper analytical bounds.
# Zero and LOQ are the boundary values for each LEFT-CENSORED observation;
# they are NOT alternative estimates or confidence intervals.
chem_scenarios <- chem_fit |>
  mutate(zero = if_else(censored, 0, concentration),
         loq_half = if_else(censored, loq / 2, concentration),
         loq = if_else(censored, loq, concentration),
         mle_gamma = value,
         bound_lower = zero,
         bound_upper = loq) |>
  select(-value, -distribution) |>
  pivot_longer(c(zero, loq_half, loq, mle_gamma,
                 bound_lower, bound_upper),
               names_to = "scenario", values_to = "value")

# Same aggregation as calc_bmfs() in script 03:
# sample -> species means -> equal mean across chemical species within phylum
# -> renormalized equal/numerical/biomass weighting across represented phyla.
keys <- c("compound_name", "scenario")
sole_means <- chem_scenarios |>
  filter(compartment == "Sole") |>
  group_by(across(all_of(keys))) |>
  summarise(c_sole = mean(value), n_sole = n(), .groups = "drop")

prey_species <- chem_scenarios |>
  filter(compartment == "Prey") |>
  left_join(taxa_mapping, by = "taxa", relationship = "many-to-one")
if (anyNA(prey_species$Phylum)) stop("Unmapped chemical prey taxa")
prey_species <- prey_species |>
  group_by(across(all_of(keys)), species, Phylum) |>
  summarise(c_prey = mean(value), n_samples = n(), .groups = "drop")

prey_phyla <- prey_species |>
  group_by(across(all_of(c(keys, "Phylum")))) |>
  summarise(c_phylum = mean(c_prey), n_species = n(), .groups = "drop") |>
  inner_join(weights, by = "Phylum", relationship = "many-to-one")

# Preserve represented phyla; do not replace missing data by a zero.
diet_means <- prey_phyla |>
  group_by(across(all_of(keys))) |>
  summarise(
    n_phyla = n(), n_species = sum(n_species),
    c_equal = mean(c_phylum),
    c_numeric = if (sum(p_numeric) > 0)
      weighted.mean(c_phylum, p_numeric) else NA_real_,
    c_biomass = if (sum(p_biomass) > 0)
      weighted.mean(c_phylum, p_biomass) else NA_real_,
    .groups = "drop"
  ) |>
  pivot_longer(starts_with("c_"), names_to = "diet_method",
               values_to = "c_diet") |>
  mutate(diet_method = str_remove(diet_method, "^c_"))

hbcdd_ratios <- sole_means |>
  inner_join(diet_means, by = keys, relationship = "many-to-many") |>
  mutate(bmf = if_else(is.finite(c_diet) & c_diet > 0,
                       c_sole / c_diet, NA_real_))

# Ratio bounds: min numerator / max denominator and vice versa.
# If dietary lower bound is zero, the upper BMF bound is unbounded;
# do not report a finite value obtained by a replacement assumption.
hbcdd_bounds <- hbcdd_ratios |>
  filter(scenario %in% c("bound_lower", "bound_upper")) |>
  select(compound_name, diet_method, scenario, c_sole, c_diet) |>
  pivot_wider(names_from = scenario,
              values_from = c(c_sole, c_diet)) |>
  mutate(
    bmf_lower = if_else(c_diet_bound_upper > 0,
                        c_sole_bound_lower / c_diet_bound_upper, NA_real_),
    bmf_upper = case_when(
      is.finite(c_diet_bound_lower) & c_diet_bound_lower > 0 ~
        c_sole_bound_upper / c_diet_bound_lower,
      c_diet_bound_lower == 0 & c_sole_bound_upper > 0 ~ Inf,
      TRUE ~ NA_real_
    ),
    interpretation = case_when(
      is.na(bmf_lower) | is.na(bmf_upper) ~ "undefined",
      bmf_upper < 1 ~ "definitely_below_1_given_LOQs",
      bmf_lower > 1 ~ "definitely_above_1_given_LOQs",
      TRUE ~ "indeterminate_given_LOQs"
    ))

# For a readable per-isomer point-estimate table.
hbcdd_table <- hbcdd_ratios |>
  filter(!str_starts(scenario, "bound_")) |>
  select(compound_name, diet_method, scenario, c_sole, c_diet, bmf) |>
  pivot_wider(names_from = scenario,
              values_from = c(c_sole, c_diet, bmf)) |>
  left_join(hbcdd_bounds, by = c("compound_name", "diet_method"),
            relationship = "one-to-one") |>
  arrange(compound_name, diet_method)

# Cross-check the gamma point estimate against the standard script 03,
# when the main analysis has already been executed.
if (file.exists("results/bmf_compare.rds")) {
  reference <- readRDS("results/bmf_compare.rds") |>
    filter(family == "HBCDD", scenario == "mle") |>
    select(compound_name, diet_method, bmf_reference = bmf)
  qa <- hbcdd_ratios |>
    filter(scenario == "mle_gamma") |>
    select(compound_name, diet_method, bmf) |>
    inner_join(reference, by = c("compound_name", "diet_method"),
               relationship = "one-to-one") |>
    mutate(absolute_diff = abs(bmf - bmf_reference))
  write_csv(qa, file.path(DIR, "HBCDD_reference_crosscheck.csv"))
  if (nrow(qa) != 9L || any(qa$absolute_diff > 1e-8, na.rm = TRUE)) {
    stop("HBCDD gamma BMFs differ from main script 03: inspect crosscheck")
  }
}

# Quantification rates, with sample counts (NOT % of concentration).
hbcdd_censor_summary <- chem |>
  group_by(compound_name, compartment) |>
  summarise(n = n(), n_quantified = sum(!censored),
            n_censored = sum(censored),
            pct_censored = 100 * mean(censored),
            loq_min = if (any(censored)) min(loq[censored]) else NA_real_,
            loq_max = if (any(censored)) max(loq[censored]) else NA_real_,
            .groups = "drop")

# The isomeric proportions calculated using substituted values are
# scenario-dependent: they cannot demonstrate chemical conversion.
# Use biomass-weighted reconstructed prey concentrations for comparison.
hbcdd_profile <- hbcdd_ratios |>
  filter(diet_method == "biomass",
         scenario %in% c("zero", "loq_half", "loq", "mle_gamma")) |>
  select(compound_name, scenario, c_sole, c_diet) |>
  pivot_longer(c(c_sole, c_diet), names_to = "compartment",
               values_to = "concentration") |>
  group_by(scenario, compartment) |>
  mutate(total = sum(concentration),
         proportion = if_else(total > 0, concentration / total, NA_real_)) |>
  ungroup()

# Diagnostic: sample_ID, season/date provenance if present in source records.
# No unverified causal inference from profile dominance.
write_csv(hbcdd_censor_summary, file.path(DIR, "HBCDD_detection_summary.csv"))
write_csv(fit_detail, file.path(DIR, "HBCDD_gamma_fit_audit.csv"))
write_csv(hbcdd_table, file.path(DIR, "HBCDD_bmf_scenarios_and_bounds.csv"))
write_csv(hbcdd_profile, file.path(DIR, "HBCDD_isomer_proportions.csv"))
write_csv(chem |> select(any_of(c("sample_id", "compartment", "group",
                                  "compound_name", "censored", "loq",
                                  "concentration", "species", "taxa"))),
          file.path(DIR, "HBCDD_analytical_provenance.csv"))
saveRDS(list(censoring = hbcdd_censor_summary, fits = fit_detail,
             ratios = hbcdd_ratios, bounds = hbcdd_bounds,
             table = hbcdd_table, profiles = hbcdd_profile),
        file.path(DIR, "HBCDD_full_analysis.rds"))


# 4. Publication-ready summary (one row per isomer and weighting method) -----
# Detection fractions use *analytical observations*, not substituted values.
detection_wide <- hbcdd_censor_summary |>
  transmute(compound_name, compartment, n,
            n_quantified, n_censored, pct_censored) |>
  pivot_wider(names_from = compartment,
              values_from = c(n, n_quantified, n_censored, pct_censored),
              names_glue = "{.value}_{compartment}")

hbcdd_table_supp <- hbcdd_table |>
  left_join(detection_wide, by = "compound_name", relationship = "many-to-one") |>
  transmute(
    Isomer = factor(compound_name, levels = ISOMERS,
                    labels = c("alpha-HBCDD", "beta-HBCDD", "gamma-HBCDD")),
    `Diet weighting` = diet_method,
    `Soles quantified (n/N)` = paste0(n_quantified_Sole, "/", n_Sole),
    `Prey quantified (n/N)` = paste0(n_quantified_Prey, "/", n_Prey),
    `BMF zero` = bmf_zero,
    `BMF LOQ/2` = bmf_loq_half,
    `BMF LOQ` = bmf_loq,
    `BMF gamma-MLE` = bmf_mle_gamma,
    `Analytical lower bound` = bmf_lower,
    `Analytical upper bound` = bmf_upper,
    `Threshold interpretation` = interpretation
  ) |>
  arrange(Isomer, factor(`Diet weighting`, levels = c("equal", "numeric", "biomass")))

# Restrict the publication table to biomass weighting, retain all methods
# in the separate full table to avoid silently omitting results.
hbcdd_table_supp_biomass <- hbcdd_table_supp |>
  filter(`Diet weighting` == "biomass") |>
  select(-`Diet weighting`)

stopifnot(nrow(hbcdd_table_supp) == 9L,
          nrow(hbcdd_table_supp_biomass) == 3L,
          all(!is.na(hbcdd_table_supp_biomass$`BMF gamma-MLE`)))

write_csv(hbcdd_table_supp, file.path(DIR, "Table_Sx_HBCDD_all_weightings.csv"))
write_csv(hbcdd_table_supp_biomass,
          file.path("manuscript/tables", "Table_Sx_HBCDD_isomers_biomass.csv"))
saveRDS(hbcdd_table_supp, file.path(DIR, "HBCDD_supplementary_table.rds"))

# 5. Supplementary figures -------------------------------------------------
FIG_DIR <- "manuscript/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)

# 5a. Detection rate plot (directly measured evidence, no substitution).
detection_plot <- hbcdd_censor_summary |>
  mutate(isomer = factor(compound_name, ISOMERS,
                         c("alpha-HBCDD", "beta-HBCDD", "gamma-HBCDD")),
         pct_quantified = 100 - pct_censored,
         compartment = factor(compartment, levels = c("Sole", "Prey")))
p_hbcdd_detection <- ggplot(detection_plot,
                             aes(x = isomer, y = pct_quantified, fill = compartment)) +
  geom_col(position = position_dodge(width = .75), width = .68) +
  geom_text(aes(label = paste0(n_quantified, "/", n)),
            position = position_dodge(width = .75), vjust = -.4, size = 3.3) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "Quantified observations (%)", fill = "Compartment") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "Fig.S_HBCDD_detection_rates.png"),
       p_hbcdd_detection, width = 7.2, height = 4.8, dpi = 300)

# 5b. Gamma-MLE point estimates against censoring-derived analytical bounds.
# Use a linear scale because beta-HBCDD can have a zero LOWER bound.
plot_bounds <- hbcdd_table |>
  mutate(isomer = factor(compound_name, ISOMERS,
                         c("alpha-HBCDD", "beta-HBCDD", "gamma-HBCDD")),
         diet_method = factor(diet_method, levels = c("equal", "numeric", "biomass")),
         upper_finite = is.finite(bmf_upper),
         upper_plot = if_else(upper_finite, bmf_upper, NA_real_))
if (any(!plot_bounds$upper_finite)) {
  warning("Some upper bounds are unbounded and are not drawn; see output table")
}
p_hbcdd_bounds <- ggplot(plot_bounds, aes(y = diet_method)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = bmf_lower, xend = upper_plot, yend = diet_method),
               linewidth = 1.25, colour = "#275662", na.rm = TRUE) +
  geom_point(aes(x = bmf_mle_gamma), size = 2.8, colour = "#D95F02") +
  facet_wrap(~ isomer, scales = "free_x") +
  labs(x = "Field-based BMF: analytical bounds and gamma-MLE point estimate",
       y = "Diet weighting") +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "Fig.S_HBCDD_BMF_analytical_bounds.png"),
       p_hbcdd_bounds, width = 10.5, height = 4.5, dpi = 300)

# 5c. Zoomed visualisation of the SAME analytical BMF bounds. ---------
# The overview plot above shows the position relative to BMF = 1, but all
# isomers lie far below one and the small analytical ranges are nearly invisible.
# This zoom intentionally does not include the BMF = 1 line so facets can
# use compound-specific x-scales. Classification vs unity remains in the table.
plot_bounds_zoom <- plot_bounds |>
  mutate(
    diet_method = factor(as.character(diet_method),
                         levels = c("equal", "numeric", "biomass")),
    label_bounds = if_else(
      is.finite(bmf_lower) & is.finite(bmf_upper),
      sprintf("%.4f–%.4f", bmf_lower, bmf_upper),
      "unbounded or unavailable"
    )
  )

p_hbcdd_bounds_zoom <- ggplot(
  plot_bounds_zoom,
  aes(y = diet_method)
) +
  geom_segment(aes(x = bmf_lower, xend = upper_plot,
                   yend = diet_method),
               linewidth = 2.3, colour = "#275662", na.rm = TRUE) +
  geom_point(aes(x = bmf_lower), shape = 124, size = 4,
             colour = "#275662", na.rm = TRUE) +
  geom_point(aes(x = upper_plot), shape = 124, size = 4,
             colour = "#275662", na.rm = TRUE) +
  geom_point(aes(x = bmf_mle_gamma), shape = 21, size = 3.3,
             fill = "#D95F02", colour = "white", stroke = 0.6) +
  facet_wrap(~ isomer, scales = "free_x", nrow = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.13)),
                     labels = scales::label_number(accuracy = 0.0001)) +
  labs(
    x = "Field BMF (zoomed within each isomer; different x-axis scales)",
    y = "Diet weighting",
    caption = paste(
      "Horizontal segments: min Csole / max Cdiet to max Csole / min Cdiet.",
      "Orange points: gamma-MLE (LOQ/2 fallback where needed).",
      "All intervals remain below BMF = 1; see overview figure and supplementary table."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8.5))

ggsave(file.path(FIG_DIR, "Fig.S_HBCDD_BMF_analytical_bounds_ZOOM.png"),
       p_hbcdd_bounds_zoom, width = 12, height = 5.1, dpi = 300)

# 5c. Print-ready Word table is optional: full precision preserved in CSV.
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  display_table <- hbcdd_table_supp_biomass |>
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  tab <- flextable::flextable(display_table) |>
    flextable::autofit()
  doc <- officer::read_docx() |>
    officer::body_add_par("Table Sx. HBCDD isomer-specific BMFs and analytical bounds",
                          style = "heading 1") |>
    flextable::body_add_flextable(tab)
  print(doc, target = file.path("manuscript/tables",
                                "Table_Sx_HBCDD_isomers_biomass.docx"))
}

print(hbcdd_censor_summary)
print(fit_detail |> select(compound_name, compartment, group,
                           fit_status, n_censored, n_uncensored))
print(hbcdd_table |> select(compound_name, diet_method,
                            bmf_loq_half, bmf_mle_gamma,
                            bmf_lower, bmf_upper, interpretation))
cat("HBCDD outputs saved to", DIR, "\n")
