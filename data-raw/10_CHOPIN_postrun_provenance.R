# CHOPIN | 10 - POST-RUN provenance: verify generated pools vs raw data.
# This script runs AFTER 01-08; it is not a bootstrap prerequisite.
# Run BEFORE trusting manuscript estimates. The source package chopin.data is required.
# This script does not modify or overwrite any analytical result.
suppressPackageStartupMessages({library(chopin.data); library(tidyverse)})
OUT <- 'results/reviewer_checks'; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
required <- c('results/scenarios.rds', 'results/bmf_compare.rds',
              'data-raw/01_CHOPIN_contamination_gamma.R')
if (any(!file.exists(required))) stop('Missing: ', paste(required[!file.exists(required)],collapse=', '))
data('soles_contam', package='chopin.data')
data('benthos_contam', package='chopin.data')
dat <- readRDS('results/scenarios.rds')
keys <- c('compartment','sample_id','family','compound_name')
stopifnot(all(c(keys,'censored','value_mle') %in% names(dat)))
if (anyDuplicated(dat[keys])) stop('Duplicated analytical records in saved scenarios')

# Read the same configuration used in script 01 (not a stale script string).
source('data-raw/00_CHOPIN_config.R', local = TRUE)
SCRIPT_SOLE_SEASON <- SOLE_SEASON
if (!SCRIPT_SOLE_SEASON %in% c('spring','autumn'))
  stop('Invalid SOLE_SEASON in 00_CHOPIN_config.R')

source_meta <- soles_contam |>
  transmute(sample_id = as.character(sample_TAG),
            source_season = as.character(season),
            source_date = as.character(date), source_age = as.character(class)) |>
  distinct() |>
  group_by(sample_id) |>
  summarise(source_season=paste(sort(unique(source_season)),collapse='; '),
            source_date=paste(sort(unique(source_date)),collapse='; '),
            source_age=paste(sort(unique(source_age)),collapse='; '), .groups='drop')
sole_inventory <- dat |>
  filter(compartment=='Sole') |>
  distinct(sample_id) |>
  left_join(source_meta, by='sample_id', relationship='many-to-one') |>
  mutate(matches_current_script = source_season == SCRIPT_SOLE_SEASON,
         missing_source_metadata = is.na(source_season))
readr::write_csv(sole_inventory, file.path(OUT,'01_sole_pool_provenance.csv'))

sample_seasons <- benthos_contam |>
  transmute(sample_id=as.character(sample_TAG), source_season=as.character(season)) |>
  distinct() |>
  count(source_season, name='n_prey_pools')
readr::write_csv(sample_seasons, file.path(OUT,'02_prey_pool_seasons.csv'))

source_script_time <- file.info('data-raw/01_CHOPIN_contamination_gamma.R')$mtime
scenarios_time <- file.info('results/scenarios.rds')$mtime
manifest <- tibble(item=c('source_script_01','stored_scenarios','source_chopin.data_version','R_version'),
                   value=c(as.character(source_script_time),as.character(scenarios_time),
                           as.character(utils::packageVersion('chopin.data')),
                           as.character(getRversion())))
readr::write_csv(manifest, file.path(OUT,'03_provenance_manifest.csv'))
writeLines(capture.output(sessionInfo()), file.path(OUT,'04_sessionInfo.txt'))
print(sole_inventory)
print(sample_seasons)
if (source_script_time > scenarios_time) {
  warning('Script 01 is newer than scenarios.rds: saved results may be stale.')
}
if (!nrow(sole_inventory) || any(sole_inventory$missing_source_metadata) ||
    any(!sole_inventory$matches_current_script, na.rm=TRUE)) {
  stop('PROVENANCE FAILURE: stored sole samples do not match the season in the current script. ',
       'Inspect results/reviewer_checks/01_sole_pool_provenance.csv. ',
       'Resolve this BEFORE rerunning the analysis and rewriting the manuscript.')
}
message('Provenance gate passed for saved sole pools. This does not verify full analytical reproducibility.')
