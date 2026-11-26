# CHOPIN | 08 - Exploratory comparison of prey seasons for the gamma-MLE scenario
# IMPORTANT: uses the full-dataset imputations from script 01 WITHOUT refitting
# the censored distribution in each season; interpretation is descriptive.
# This script does not repair non-contemporaneous sole/diet sampling.
suppressPackageStartupMessages({library(chopin.data);library(tidyverse)})
OUT <- 'results/reviewer_checks';dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
data('benthos_contam',package='chopin.data')
chem <- readRDS('results/scenarios.rds')
base <- readRDS('results/bmf_compare.rds') |>
  filter(scenario=='mle',diet_method=='biomass') |>
  select(family,compound_name,unit,c_sole,bmf_pooled=bmf)
w <- readRDS('results/prop_diet.rds') |>
  transmute(Phylum,weight=p_biomass_phylum)
map <- tribble(~taxa,~Phylum,'Polychaetes','Annelida','Crustaceans','Arthropoda','Bivalves','Mollusca')
meta <- benthos_contam |>
  transmute(sample_id=as.character(sample_TAG),season=as.character(season)) |>
  distinct()
if (anyDuplicated(meta['sample_id'])) stop('One prey sample has multiple season labels; resolve metadata')
prey <- chem |>
  filter(compartment=='Prey') |>
  left_join(meta,by='sample_id',relationship='many-to-one') |>
  left_join(map,by='taxa',relationship='many-to-one') |>
  filter(season %in% c('spring','autumn'))
if (anyNA(prey$Phylum)) stop('Missing phylum mapping')
keys <- c('family','compound_name','unit')
seasonal <- prey |>
  group_by(across(all_of(keys)),season,Phylum,species) |>
  summarise(c_species=mean(value_mle),n_pools=n(),.groups='drop') |>
  group_by(across(all_of(keys)),season,Phylum) |>
  summarise(c_phylum=mean(c_species),n_species=n_distinct(species),
            n_pools=sum(n_pools),.groups='drop') |>
  left_join(w,by='Phylum',relationship='many-to-one') |>
  group_by(across(all_of(keys)),season) |>
  summarise(c_diet=sum(weight*c_phylum)/sum(weight),
            n_phyla=n(),n_species=sum(n_species),n_pools=sum(n_pools),
            weight_coverage=sum(weight),.groups='drop') |>
  left_join(base,by=keys,relationship='many-to-one') |>
  mutate(bmf_season=c_sole/c_diet,
         pct_diff_vs_pooled=100*(bmf_season/bmf_pooled-1),
         threshold_change=(bmf_season>1)!=(bmf_pooled>1))
readr::write_csv(seasonal,file.path(OUT,'BMF_seasonal_prey_exploratory.csv'))

# Supplementary summary: season-restricted prey chemistry

dir.create(
  "manuscript/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

season_summary <- seasonal |>
  group_by(family, season) |>
  summarise(
    n_compounds = n_distinct(compound_name),
    n_class_changes = sum(threshold_change, na.rm = TRUE),
    n_species_min = min(n_species),
    n_species_max = max(n_species),
    n_pools_min = min(n_pools),
    n_pools_max = max(n_pools),
    n_phyla_min = min(n_phyla),
    n_phyla_max = max(n_phyla),
    weight_coverage_min = min(weight_coverage),
    weight_coverage_max = max(weight_coverage),
    .groups = "drop"
  ) |>
  arrange(family, season)

readr::write_csv(
  season_summary,
  "manuscript/tables/Table_Sx_prey_season_sensitivity.csv"
)

print(seasonal |> group_by(family,season) |>
        summarise(n_compounds=n(),n_flips=sum(threshold_change),
                  max_abs_pct=max(abs(pct_diff_vs_pooled)),.groups='drop'))
# Do not infer season from sampling dates if the supplied season field is missing.
