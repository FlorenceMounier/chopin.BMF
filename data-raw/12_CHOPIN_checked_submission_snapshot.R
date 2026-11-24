# CHOPIN | 12 - New submission snapshot only after passing provenance checks
# Run only AFTER rerunning the original production pipeline + reviewer checks.
suppressPackageStartupMessages(library(tidyverse))
source('data-raw/10_CHOPIN_postrun_provenance.R', local=TRUE) # Check regenerated vs source pools
RDS <- c('results/bmf_compare.rds','results/scenarios.rds',
         'results/gamma_vs_lognormal/reference_check.rds')
if (any(!file.exists(RDS))) stop('Missing RDS data: ',paste(RDS[!file.exists(RDS)],collapse=', '))
c <- readRDS(RDS[1]);rc <- readRDS(RDS[3])
if (nrow(rc)!=nrow(dplyr::filter(c, scenario=='mle')) || any(!is.finite(rc$abs_diff)) ||
    any(rc$abs_diff > 1e-8*pmax(1,rc$bmf_main)))
  stop('Gamma BMFs do not reproduce script 03 results')
expected <- c('zero','loq_half','loq','mle')
if (!setequal(unique(c$scenario),expected)) stop('Unexpected censoring scenarios')
n_compounds <- nrow(dplyr::distinct(c, family, compound_name, unit))
if (nrow(c) != 12L * n_compounds)
  stop('Expected 4 censoring treatments x 3 diet weights per retained compound')

## SNAPSHOT PATH

TAG <- format(Sys.Date(),'%Y%m%d')
# Parent directory of the chopin.BMF project
PROJECT_PARENT <- dirname(
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
)
# Dedicated external directory for frozen snapshots
SNAPSHOT_ROOT <- file.path(
  PROJECT_PARENT,
  "CHOPIN_submission_snapshots"
)
# Create archive directory if necessary
dir.create(
  SNAPSHOT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)
# Snapshot destination
OUT <- file.path(
  SNAPSHOT_ROOT,
  paste0("CHOPIN_validated_", TAG)
)
# Do not overwrite an existing snapshot
if (dir.exists(OUT)) {
  stop(
    "Snapshot already exists: ", OUT,
    "\nChoose another date or archive name."
  )
}
message("Snapshot destination: ", OUT)
if (dir.exists(OUT)) stop('Snapshot already exists: ',OUT)


paths <- c('DESCRIPTION','README.md','README_BUILD_FROM_ZERO.md','manuscript/Article.Rmd','manuscript/Article_SI.Rmd',
           'manuscript/references.bib','manuscript/template.docx',
           'manuscript/template_SI.docx','inst/extdata/Gogina_data-2021-0002-01.xlsx')
paths <- unique(c(paths,list.files('data-raw',full.names=TRUE,pattern='\\.R$'),
                  list.files('results',recursive=TRUE,full.names=TRUE,pattern='\\.(csv|rds)$'),
                  list.files('manuscript/tables',full.names=TRUE),
                  list.files('manuscript/figures',full.names=TRUE),
                  list.files('data', full.names=TRUE, pattern='\\.rda$')))
paths <- paths[file.exists(paths) & !dir.exists(paths)]
dir.create(OUT,recursive=TRUE)
for (f in paths) {
 target <- file.path(OUT,f)
 dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
 if (!file.copy(f,target,overwrite=FALSE)) stop('Cannot freeze ',f)
}
md5 <- tools::md5sum(file.path(OUT,paths))
manifest <- tibble(file=paths,MD5=unname(md5),
                   source_chopin_data_version=as.character(utils::packageVersion('chopin.data')))
write_csv(manifest,file.path(OUT,'manifest_MD5.csv'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
message('Immutable new snapshot created: ',OUT)
