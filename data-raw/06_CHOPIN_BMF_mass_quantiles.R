# CHOPIN | 06 - Propagate prey mass Q25 / median / Q75 to diet-based BMFs
# Run from package root AFTER 01 (chemistry), 02 (diet) and 03 (BMFs).
# Reference: gamma-MLE censoring scenario and biomass diet weighting.
# Q25 and Q75 are simultaneously applied to all prey taxa.
# These are DETERMINISTIC SENSITIVITY SCENARIOS, NOT 95% confidence intervals.
# They do not propagate analytical variance, sample uncertainty, or covariance.

suppressPackageStartupMessages(library(tidyverse))
OUT <- "results/mass_quantiles"
FIG <- "manuscript/figures"
TAB <- "manuscript/tables"
for (d in c(OUT, FIG, TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

required <- c("results/bmf_phylum.rds", "results/bmf_compare.rds",
              "results/prop_diet.rds")
missing <- required[!file.exists(required)]
if (length(missing)) stop("Run scripts 01-03 first. Missing: ", paste(missing, collapse=", "))

keys <- c("family", "compound_name", "unit")
phylum <- readRDS(required[1]) |>
  filter(scenario == "mle") |>
  select(all_of(keys), Phylum, c_phylum)
reference <- readRDS(required[2]) |>
  filter(scenario == "mle", diet_method == "biomass") |>
  select(all_of(keys), c_sole, bmf_reference = bmf)
weights <- readRDS(required[3]) |>
  transmute(Phylum,
            Q25 = p_biomass_q25_phylum,
            Median = p_biomass_phylum,
            Q75 = p_biomass_q75_phylum) |>
  pivot_longer(c(Q25, Median, Q75),
               names_to = "mass_assumption", values_to = "weight")

# Validate phylum weights independently for each mass hypothesis.
stopifnot(!anyDuplicated(phylum[c(keys, "Phylum")]),
          !anyDuplicated(reference[keys]),
          !anyDuplicated(weights[c("Phylum", "mass_assumption")]),
          all(is.finite(weights$weight)), all(weights$weight >= 0))
weight_check <- weights |>
  group_by(mass_assumption) |>
  summarise(weight_sum = sum(weight), .groups = "drop")
if (any(abs(weight_check$weight_sum - 1) > 1e-6)) {
  warning("Some mass-quantile phylum weights do not sum to one; check diet script 02")
}
if (!nrow(reference)) stop("No gamma-MLE + biomass reference BMFs found")

# Same denominator as script 03: equal mean over measured chemical prey
# species within phylum, then weighted average over represented phyla.
# If some phyla lack chemical data for a compound, renormalize the weights.
scenario <- phylum |>
  inner_join(weights, by = "Phylum", relationship = "many-to-many") |>
  group_by(across(all_of(keys)), mass_assumption) |>
  summarise(
    n_phyla = n_distinct(Phylum),
    weight_coverage = sum(weight),
    c_diet = if (sum(weight) > 0) sum(weight * c_phylum) / sum(weight)
             else NA_real_,
    .groups = "drop"
  ) |>
  left_join(reference, by = keys, relationship = "many-to-one") |>
  mutate(
    bmf = if_else(is.finite(c_diet) & c_diet > 0, c_sole / c_diet, NA_real_),
    delta_from_reference_pct = if_else(is.finite(bmf_reference) & bmf_reference > 0,
                                       100 * (bmf / bmf_reference - 1), NA_real_),
    above_unity = if_else(is.finite(bmf), bmf > 1, NA)
  )

# Integrity check: quantile = Median reproduces the original reference.
med <- scenario |> filter(mass_assumption == "Median")
if (nrow(med) != nrow(reference) || anyNA(med$bmf) ||
    any(abs(med$bmf - med$bmf_reference) >
        1e-8 * pmax(1, abs(med$bmf_reference)))) {
  stop("Mass median scenario differs from script 03: check aggregation / keys")
}
if (anyNA(scenario$bmf)) warning("Some BMF mass-quantile scenarios are undefined")

# Q25 and Q75 do NOT guarantee an ordered BMF interval: increasing all
# individual masses changes phylum weights non-monotonically.

summary <- scenario |>
  group_by(across(all_of(keys))) |>
  summarise(
    bmf_q25 = bmf[mass_assumption == "Q25"],
    bmf_median = bmf[mass_assumption == "Median"],
    bmf_q75 = bmf[mass_assumption == "Q75"],

    n_valid = sum(is.finite(bmf)),
    complete = n() == 3L &&
      n_distinct(mass_assumption) == 3L &&
      n_valid == 3L,

    bmf_sensitivity_min = if (complete) min(bmf)
    else NA_real_,

    bmf_sensitivity_max = if (complete) max(bmf)
    else NA_real_,

    .groups = "drop"
  ) |>
  mutate(
    fold_range = if_else(
      complete & bmf_sensitivity_min > 0,
      bmf_sensitivity_max / bmf_sensitivity_min,
      NA_real_
    ),
    crosses_unity = if_else(
      complete,
      bmf_sensitivity_min < 1 &
        bmf_sensitivity_max > 1,
      NA
    ),
    touches_unity = if_else(
      complete,
      bmf_sensitivity_min <= 1 &
        bmf_sensitivity_max >= 1,
      NA
    ),
    reference_class = case_when(
      bmf_median > 1 ~ "Above 1",
      bmf_median < 1 ~ "Below 1",
      bmf_median == 1 ~ "Equal to 1",
      TRUE ~ "Undefined"
    )
  )


family_summary <- summary |>
  group_by(family) |>
  summarise(n_compounds = n(),
            median_fold_range = median(fold_range, na.rm = TRUE),
            maximum_fold_range = max(fold_range, na.rm = TRUE),
            n_crossing_unity = sum(crosses_unity, na.rm = TRUE),
            n_touching_unity = sum(touches_unity, na.rm = TRUE),
            .groups = "drop")

# Concise supplementary table: all compounds, including stable ones.
table_supp <- summary |>
  transmute(Family = family, Compound = compound_name,
            `BMF Q25` = bmf_q25, `BMF median mass` = bmf_median,
            `BMF Q75` = bmf_q75,
            `Sensitivity minimum` = bmf_sensitivity_min,
            `Sensitivity maximum` = bmf_sensitivity_max,
            `Fold range` = fold_range,
            `Crosses BMF = 1` = crosses_unity) |>
  arrange(Family, Compound)

write_csv(scenario, file.path(OUT, "BMF_mass_quantile_scenarios.csv"))
write_csv(summary, file.path(OUT, "BMF_mass_quantile_per_compound.csv"))
write_csv(family_summary, file.path(OUT, "BMF_mass_quantile_family.csv"))
write_csv(table_supp, file.path(TAB, "Table_Sx_BMF_mass_quantiles.csv"))
saveRDS(list(scenarios=scenario, per_compound=summary,
             family=family_summary, supplementary=table_supp),
        file.path(OUT, "BMF_mass_quantiles.rds"))

# Fig Sx-1: sensitivity ranges (three deterministic scenarios) by compound.
# Segment joins min/max among Q25, median and Q75; circle = median mass.
plot_data <- summary |>
  mutate(compound_name = factor(compound_name,
                                levels = rev(sort(unique(compound_name)))))

if (any(summary$bmf_sensitivity_min <= 0, na.rm = TRUE)) {
  warning(
    "Some sensitivity ranges include zero and cannot be ",
    "fully displayed on a log10 axis."
  )
}

p_mass_range <- ggplot(plot_data, aes(y = compound_name)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = bmf_sensitivity_min, xend = bmf_sensitivity_max,
                   yend = compound_name), linewidth = .95, colour = "#275662") +
  geom_point(aes(x = bmf_median, shape = crosses_unity), size = 2.2,
             colour = "#D95F02") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  labs(x = "Biomass-weighted gamma-MLE BMF (log scale)",
       y = NULL, shape = "Crosses BMF = 1") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "Fig.S_mass_quantile_BMF_ranges.png"),
       p_mass_range, width = 11, height = 10.5, dpi = 300)

# Fig Sx-2: effect of the three assumptions for compounds changing
# position relative to unity. Do not invent values if there are no switches.
focus <- summary |> filter(touches_unity) |> select(all_of(keys))
if (nrow(focus)) {
  focus_data <- scenario |> semi_join(focus, by = keys) |>
    mutate(mass_assumption = factor(mass_assumption,
                                    levels = c("Q25", "Median", "Q75")))
  p_mass_focus <- ggplot(focus_data,
                         aes(x = mass_assumption, y = bmf, group = compound_name)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
    geom_line(alpha = .6) + geom_point(size = 2) +
    facet_wrap(~ compound_name, scales = "free_y") +
    labs(x = "Individual prey mass assumption",
         y = "Biomass-weighted gamma-MLE BMF") +
    theme_bw(base_size = 11)
  ggsave(file.path(FIG, "Fig.S_mass_quantile_threshold_sensitive.png"),
         p_mass_focus, width = 10, height = 7, dpi = 300)
} else {
  # Prevent a figure from a previous run from being mistaken for current results.
  old_figure <- file.path(FIG, "Fig.S_mass_quantile_threshold_sensitive.png")
  if (file.exists(old_figure)) file.remove(old_figure)
}

# Optional publication-ready Word table (CSV retains full precision).
if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {
  ft <- table_supp |>
    mutate(across(where(is.numeric), ~ round(.x, 3))) |>
    flextable::flextable() |>
    flextable::autofit()
  doc <- officer::read_docx() |>
    officer::body_add_par("Table Sx. Prey-mass quantile sensitivity", style = "heading 1") |>
    flextable::body_add_flextable(ft)
  print(doc, target = file.path(TAB, "Table_Sx_BMF_mass_quantiles.docx"))
}

print(family_summary)
print(summary |> filter(touches_unity))
cat("Mass-quantile results saved in", OUT, "\n")
