# CHOPIN | 03b - Eligibility and zero-compatible methodological sensitivity
# Run after 03_CHOPIN_BMF_gamma.R. No pre-existing results beyond script 03.
# A zero BMF is valid descriptively, but cannot enter log2(max/min).
# Bounded contrast D=(max-min)/(max+min): 0=no difference, 1=one extreme is zero.
# If both extremes are zero, D=0 by explicit convention (no variation).
# This is a deterministic sensitivity score, not a confidence interval.

suppressPackageStartupMessages(library(tidyverse))
input <- "results/bmf_compare.rds"
if (!file.exists(input)) stop("Run script 03 first: missing ", input)
B <- readRDS(input)
keys <- c("family", "compound_name", "unit")
expected_scenarios <- c("zero", "loq_half", "loq", "mle")
expected_diets <- c("equal", "numeric", "biomass")
stopifnot(all(c(keys, "scenario", "diet_method", "bmf") %in% names(B)))
if (anyDuplicated(B[c(keys, "scenario", "diet_method")])) {
  stop("Duplicate BMF combinations; resolve before sensitivity calculations")
}
if (any(!B$scenario %in% expected_scenarios) ||
    any(!B$diet_method %in% expected_diets)) {
  stop("Unexpected censoring or diet method")
}

eligibility <- B |>
  group_by(across(all_of(keys))) |>
  summarise(
    n_combinations = n(),
    n_nonnegative = sum(is.finite(bmf) & bmf >= 0),
    n_positive = sum(is.finite(bmf) & bmf > 0),
    n_zero = sum(is.finite(bmf) & bmf == 0),
    n_invalid = sum(!is.finite(bmf) | bmf < 0),
    full_grid = n() == 12L &&
      n_distinct(paste(scenario, diet_method, sep = "::")) == 12L,
    .groups = "drop"
  ) |>
  mutate(
    metric_estimable = full_grid & n_positive == 12L,
    bounded_metric_estimable = full_grid & n_nonnegative == 12L,
    non_calculation_reason = case_when(
      metric_estimable ~ "Estimable",
      !full_grid ~ "Incomplete 4 x 3 combinations",
      n_invalid > 0L ~ "Non-finite, missing or negative BMF",
      n_zero > 0L ~ "Zero BMF: logarithmic ratio undefined",
      TRUE ~ "Other"
    )
  )

# The denominator max + min is nonzero unless both BMFs are zero.
# In that exceptional case, declare zero contrast because neither BMF varies.
bounded_contrast <- function(x) {
  if (!length(x) || any(!is.finite(x) | x < 0)) return(NA_real_)
  lo <- min(x); hi <- max(x)
  if (hi == 0) return(0)
  (hi - lo) / (hi + lo)
}

# Comparable complete grids only. A zero is allowed, but NA/Inf/negative isn't.
B_valid <- B |>
  inner_join(eligibility |> filter(bounded_metric_estimable) |>
               select(all_of(keys)), by = keys)

by_censor_scenario <- B_valid |>
  group_by(across(all_of(keys)), scenario) |>
  summarise(D_diet = bounded_contrast(bmf), .groups = "drop")
by_diet_weight <- B_valid |>
  group_by(across(all_of(keys)), diet_method) |>
  summarise(D_censor = bounded_contrast(bmf), .groups = "drop")

zero_safe_sensitivity <- by_censor_scenario |>
  group_by(across(all_of(keys))) |>
  summarise(D_diet = median(D_diet), D_diet_max = max(D_diet),
            .groups = "drop") |>
  inner_join(by_diet_weight |>
    group_by(across(all_of(keys))) |>
    summarise(D_censor = median(D_censor),
              D_censor_max = max(D_censor), .groups = "drop"),
    by = keys, relationship = "one-to-one") |>
  left_join(eligibility, by = keys, relationship = "one-to-one") |>
  mutate(dominant_source_bounded = case_when(
    D_diet > D_censor ~ "Diet weighting",
    D_censor > D_diet ~ "Censoring treatment",
    TRUE ~ "Equal sensitivity"))

# Explicit status row for every retained chemical (including non-estimable ones)
all_compounds <- eligibility |>
  left_join(zero_safe_sensitivity |>
    select(all_of(keys), D_diet, D_diet_max, D_censor, D_censor_max,
           dominant_source_bounded),
    by = keys, relationship = "one-to-one")

family_eligibility <- all_compounds |>
  group_by(family) |>
  summarise(
    n_compounds = n(),
    n_log2_eligible = sum(metric_estimable),
    n_log2_not_estimable = sum(!metric_estimable),
    n_bounded_eligible = sum(bounded_metric_estimable),
    n_bounded_not_estimable = sum(!bounded_metric_estimable),
    n_zero_related_exclusions = sum(!metric_estimable & n_zero > 0),
    median_D_diet = if (all(is.na(D_diet))) NA_real_ else median(D_diet, na.rm = TRUE),
    median_D_censor = if (all(is.na(D_censor))) NA_real_ else median(D_censor, na.rm = TRUE),
    .groups = "drop"
  )

for (pair in list(
  list(all_compounds, "sensitivity_metric_eligibility"),
  list(zero_safe_sensitivity, "sensitivity_zero_safe"),
  list(family_eligibility, "sensitivity_family_eligibility"),
  list(by_censor_scenario, "sensitivity_zero_safe_by_scenario"),
  list(by_diet_weight, "sensitivity_zero_safe_by_diet_method")
)) {
  name <- pair[[2]]; dat <- pair[[1]]
  readr::write_csv(dat, file.path("results", paste0(name, ".csv")))
  saveRDS(dat, file.path("results", paste0(name, ".rds")))
}

readr::write_csv(
  sensitivity_family_eligibility,
  "manuscript/tables/Table_Sx_Sensitivity_eligibility.csv"
)

# Optional supplementary visual, no forced figure in the main manuscript.
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)
p_zero_safe <- ggplot(zero_safe_sensitivity,
                      aes(x = D_diet, y = D_censor)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(aes(shape = !metric_estimable), size = 2.4) +
  facet_wrap(~family) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Log2 eligible",
                                `TRUE` = "Zero prevents log2"),
                     name = "Log2 sensitivity") +
  labs(x = "Diet weighting sensitivity (bounded contrast)",
       y = "Censoring sensitivity (bounded contrast)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("manuscript/figures/Fig.S_zero_safe_sensitivity.png",
       p_zero_safe, width = 10, height = 5, dpi = 300)

print(family_eligibility)
print(all_compounds |>
        filter(family == "HBCDD") |>
        select(compound_name, metric_estimable, non_calculation_reason,
               bounded_metric_estimable, D_diet, D_censor))
cat("03b completed. Outputs: results/sensitivity_* and Fig.S_zero_safe_sensitivity.png\n")
