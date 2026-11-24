# CHOPIN | 07 - Leave-one-chemically-analysed-prey-species-out BMF sensitivity
# Gamma-MLE, biomass weighting; means: pooled samples -> species -> phylum.
# Removing the only analysed species in a phylum drops that phylum, and
# renormalises remaining dietary weights. Such removals are explicitly flagged.
suppressPackageStartupMessages(library(tidyverse))
OUT <- 'results/reviewer_checks'; dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
prey <- readRDS('results/bmf_prey.rds') |>
  filter(scenario=='mle',is.finite(c_prey))
base <- readRDS('results/bmf_compare.rds') |>
  filter(scenario=='mle',diet_method=='biomass') |>
  select(family,compound_name,unit,c_sole,bmf_reference=bmf)
w <- readRDS('results/prop_diet.rds') |>
  transmute(Phylum,weight=p_biomass_phylum)
keys <- c('family','compound_name','unit')
if (anyDuplicated(base[keys]) || anyDuplicated(w['Phylum'])) stop('Non-unique keys')

# Use exactly the same reference denominator as script 03.
recalc_bmf <- function(species_table,sole) {
  by_phylum <- species_table |>
    group_by(Phylum) |>
    summarise(c_phylum=mean(c_prey),.groups='drop') |>
    left_join(w,by='Phylum',relationship='many-to-one')
  if (!nrow(by_phylum) || any(!is.finite(by_phylum$weight))) return(NA_real_)
  c_diet <- sum(by_phylum$c_phylum*by_phylum$weight)/sum(by_phylum$weight)
  if (!is.finite(c_diet) || c_diet<=0) NA_real_ else sole/c_diet
}
parts <- prey |> group_by(across(all_of(keys))) |> group_split()
res <- purrr::map_dfr(parts,function(d) {
  b <- base |> semi_join(d |> distinct(across(all_of(keys))),by=keys)
  if (nrow(b)!=1) stop('Cannot uniquely identify reference BMF')
  ref <- recalc_bmf(d,b$c_sole)
  if (!isTRUE(all.equal(ref,b$bmf_reference,tolerance=1e-8)))
    stop('LOO baseline differs from script 03 BMF: ',b$compound_name)
  purrr::map_dfr(seq_len(nrow(d)),function(i) {
    n_in_phylum <- sum(d$Phylum==d$Phylum[i])
    left <- d[-i,,drop=FALSE]
    trial <- recalc_bmf(left,b$c_sole)
    tibble(family=b$family,compound_name=b$compound_name,unit=b$unit,
           removed_species=as.character(d$species[i]),
           removed_phylum=as.character(d$Phylum[i]),
           phylum_removed=n_in_phylum==1L,
           n_species_remaining=nrow(left),
           bmf_reference=b$bmf_reference,bmf_without_species=trial,
           delta_pct=100*(trial/b$bmf_reference-1),
           threshold_change=is.finite(trial) & ((trial>1)!=(b$bmf_reference>1)))
  })
})
by_compound <- res |>
  group_by(across(all_of(keys))) |>
  summarise(n_removals=n(),max_abs_delta_pct=max(abs(delta_pct),na.rm=TRUE),
            n_threshold_changes=sum(threshold_change),
            n_phylum_drops=sum(phylum_removed),.groups='drop')
by_family <- by_compound |>
  group_by(family) |>
  summarise(
    n_compounds = n(),
    n_removal_tests = sum(n_removals),

    # Must be computed before aggregating n_threshold_changes
    n_compounds_with_flip = sum(n_threshold_changes > 0),

    n_threshold_changes = sum(n_threshold_changes),
    n_phylum_drops = sum(n_phylum_drops),
    max_abs_delta_pct = max(max_abs_delta_pct),
    .groups = "drop"
  )
readr::write_csv(res,file.path(OUT,'LOO_species_detailed.csv'))
readr::write_csv(by_compound,file.path(OUT,'LOO_species_by_compound.csv'))
readr::write_csv(by_family,file.path(OUT,'LOO_species_by_family.csv'))
print(by_family)
print(res |> filter(threshold_change) |> arrange(family,compound_name))


# --------------------------------------------------------------------------
# 07b. SI-ready quantitative summaries and figures (reference = gamma MLE,
# biomass weighting). This is leave-one-CHEMICALLY-ANALYSED-species-out,
# not a resampling analysis of identified stomach-content prey.
# --------------------------------------------------------------------------
TAB_DIR <- 'manuscript/tables'
FIG_DIR <- 'manuscript/figures'
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Quality checks and study-wide counts
stopifnot(nrow(res) == sum(by_family$n_removal_tests))
stopifnot(sum(res$threshold_change) == sum(by_family$n_threshold_changes))
stopifnot(res |> filter(threshold_change) |> distinct(family, compound_name) |> nrow() ==
          sum(by_family$n_compounds_with_flip))

loo_overall <- tibble(
  n_compounds = nrow(by_compound),
  n_removal_tests = nrow(res),
  n_threshold_changes = sum(res$threshold_change),
  n_compounds_with_flip = res |>
    filter(threshold_change) |>
    distinct(family, compound_name) |>
    nrow(),
  n_phylum_drops = sum(res$phylum_removed)
)

# One row per chemically analysed species: number of tests, class flips,
# and effect sizes. A test here is one compound-by-species deletion.
loo_species_summary <- res |>
  mutate(abs_delta_pct = abs(delta_pct)) |>
  group_by(removed_species, removed_phylum) |>
  summarise(
    n_tests = n(),
    n_threshold_changes = sum(threshold_change),
    n_compounds_with_flip = n_distinct(compound_name[threshold_change]),
    n_phylum_drops = sum(phylum_removed),
    median_abs_delta_pct = median(abs_delta_pct, na.rm = TRUE),
    max_abs_delta_pct = max(abs_delta_pct, na.rm = TRUE),
    .groups = 'drop'
  ) |>
  arrange(desc(n_threshold_changes), desc(median_abs_delta_pct))

# An explicit table of cases where the LAST chemically analysed prey
# species of a phylum is removed. The phylum disappears from the chemical
# denominator; weights for remaining phyla are renormalised.
loo_phylum_drops <- res |>
  filter(phylum_removed) |>
  select(family, compound_name, removed_species, removed_phylum,
         bmf_reference, bmf_without_species, delta_pct, threshold_change)

loo_compound_summary <- by_compound |>
  left_join(base, by = c('family', 'compound_name', 'unit'),
            relationship = 'one-to-one') |>
  mutate(threshold_sensitive = n_threshold_changes > 0) |>
  arrange(family, desc(max_abs_delta_pct))

readr::write_csv(loo_overall, file.path(OUT, 'LOO_overall.csv'))
readr::write_csv(loo_species_summary, file.path(OUT, 'LOO_species_influence.csv'))
readr::write_csv(loo_phylum_drops, file.path(OUT, 'LOO_phylum_drops.csv'))
readr::write_csv(loo_compound_summary, file.path(OUT, 'LOO_compound_SI.csv'))

# Stable SI paths under manuscript/tables: can be read from Article_SI.Rmd
# using tables/..., as in all other flextable chunks.
readr::write_csv(by_family, file.path(TAB_DIR, 'Table_Sx_LOO_by_family.csv'))
readr::write_csv(loo_species_summary, file.path(TAB_DIR, 'Table_Sx_LOO_species_influence.csv'))
readr::write_csv(loo_phylum_drops, file.path(TAB_DIR, 'Table_Sx_LOO_phylum_drops.csv'))
readr::write_csv(loo_compound_summary, file.path(TAB_DIR, 'Table_Sx_LOO_by_compound.csv'))

# Figure 1: which analytically characterised species drive threshold flips?
p_loo_species <- loo_species_summary |>
  mutate(species_label = stringr::str_replace_all(removed_species, '_', ' '),
         species_label = forcats::fct_reorder(species_label, n_threshold_changes)) |>
  ggplot(aes(x = n_threshold_changes, y = species_label, fill = removed_phylum)) +
  geom_col(width = .72) +
  geom_text(aes(label = n_threshold_changes), hjust = -.15, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = 'Number of BMF = 1 classification changes',
       y = 'Chemically analysed prey species', fill = 'Phylum') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = 'bottom')
ggsave(file.path(FIG_DIR, 'Fig.S_LOO_species_threshold_changes.png'),
       p_loo_species, width = 9, height = 6.5, dpi = 300)

# Figure 2: distribution of change magnitude by compound; use percentage
# change, which is descriptive and not a confidence interval.
p_loo_compound <- loo_compound_summary |>
  mutate(compound_label = forcats::fct_reorder(compound_name, max_abs_delta_pct)) |>
  ggplot(aes(x = max_abs_delta_pct, y = compound_label,
             shape = threshold_sensitive)) +
  geom_point(size = 2.5) +
  facet_grid(family ~ ., scales = 'free_y', space = 'free_y') +
  labs(x = 'Maximum absolute BMF change after one species removal (%)',
       y = 'Compound', shape = 'Crosses BMF = 1') +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
ggsave(file.path(FIG_DIR, 'Fig.S_LOO_compound_max_change.png'),
       p_loo_compound, width = 10, height = 10, dpi = 300)

print(loo_overall)
print(by_family)
print(head(loo_species_summary, 12))
cat('07 LOO: SI tables and two figures exported.\n')
