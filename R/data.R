# Dataset documentation for the 31 packaged analysis objects.
# Automatically prepared from CHOPIN scripts, actual packaged RData
# column names and supplied minimal Rd documentation.
# Keep this source in R/data.R and run devtools::document() after editing.

#' Detailed compound-level censoring and gamma-fit audit
#'
#' Audit by compound, biological compartment and analytical group, combining the censoring inventory, selection decisions, fitted gamma parameters and fallback diagnostics.
#'
#' @format A data frame with 746 rows and 34 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{compartment}}{Biological compartment, Sole or Prey.}
#'   \item{\code{group}}{Chemical fitting stratum: sole group or chemically analysed prey species.}
#'   \item{\code{n}}{Number of analytical sample-compound records for a fitting stratum.}
#'   \item{\code{n_censored}}{Number of censored analytical records.}
#'   \item{\code{n_uncensored}}{Number of uncensored analytical records.}
#'   \item{\code{pct_censored}}{Percentage of censored analytical records, 100 times n_censored divided by n.}
#'   \item{\code{loq_min}}{Minimum LOQ among censored records in the fitting group (same compound-specific concentration units).}
#'   \item{\code{loq_max}}{Maximum LOQ among censored records in the fitting group (same compound-specific concentration units).}
#'   \item{\code{n_compartments}}{Number of represented biological compartments (expected maximum: two).}
#'   \item{\code{high_censoring}}{TRUE when at least one compartment has a censoring percentage at or above the configured 80-percent threshold.}
#'   \item{\code{all_censored_in_compartment}}{TRUE when at least one represented compartment has no quantified observations for the compound.}
#'   \item{\code{censored_in_sole_pct}}{Percentage censored among age-zero sole chemical pools.}
#'   \item{\code{censored_in_prey_pct}}{Percentage censored among benthic prey chemical samples.}
#'   \item{\code{excluded}}{TRUE when the compound is excluded from current BMF estimates under the selection rules.}
#'   \item{\code{exclusion_reason}}{Human-readable reason for exclusion or retention, including the HBCDD exception.}
#'   \item{\code{evidence_tier}}{Evidence designation: excluded, standard, or exploratory HBCDD retention due to censoring.}
#'   \item{\code{distribution}}{Fitted concentration distribution (gamma in the packaged fit_audit object).}
#'   \item{\code{fit_status}}{Status of likelihood optimisation: fitted, insufficient_data, failed or excluded, as appropriate.}
#'   \item{\code{fallback_method}}{Alternative imputation used when maximum-likelihood fitting is not available (LOQ/2).}
#'   \item{\code{fit_error}}{Diagnostic message explaining a failed/infeasible fit when applicable.}
#'   \item{\code{n_imputed}}{Number of censored records requiring imputation in the fitting group.}
#'   \item{\code{n_imputed_mle}}{Number of censored records imputed using successful gamma maximum-likelihood fitting.}
#'   \item{\code{n_imputed_fallback}}{Number of censored records assigned a fallback value (LOQ/2).}
#'   \item{\code{n_unresolved}}{Number of records with non-finite imputed concentrations.}
#'   \item{\code{imputed_min}}{Minimum imputed censored concentration (same concentration units as the compound).}
#'   \item{\code{imputed_median}}{Median imputed censored concentration (same concentration units as the compound).}
#'   \item{\code{imputed_max}}{Maximum imputed censored concentration (same concentration units as the compound).}
#'   \item{\code{mle_distribution}}{Name of imputation distribution (gamma for included compounds).}
#'   \item{\code{gamma_shape}}{Estimated gamma shape parameter (dimensionless), when the gamma fit succeeds.}
#'   \item{\code{gamma_scale}}{Estimated gamma scale parameter (same units as compound concentration), when fitting succeeds.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{fallback_reason}}{Reason for an LOQ/2 fallback, such as insufficient data or fit failure.}
#' }
#'
#' @details Each row corresponds to a chemical compound within a compartment and a biological group (the sole group or an analysed prey species). The audit starts from all prepared chemical records, including compounds later excluded from BMF computation. A compound is excluded if it is missing from a biological compartment or, for PCB and PFAS, at least 80 percent of observations in either compartment are censored. HBCDD isomers are an explicit exploratory exception to the 80 percent censoring rule. Gamma parameters are reported only when fitting succeeds; censored values otherwise receive the predefined LOQ/2 fallback. The count and percentage variables are based on analytical records, not independent fish.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"audit_final"

#' Summary of gamma-imputation audit outcomes
#'
#' Two-column key/value summary of the number and percentage of compound-by-group gamma fits requiring substitution, succeeding or falling back to LOQ/2.
#'
#' @format A data frame with 8 rows and 2 columns.
#' \describe{
#'   \item{\code{name}}{Name of a summary metric or imputation fallback reason.}
#'   \item{\code{value}}{Value of the named summary metric: a count or a percentage according to name.}
#' }
#'
#' @details The name column includes n_total, n_requiring_imputation, n_fitted, n_fallback, fallback_pct_total, fallback_pct_censored, and named fallback reasons when present. Values are counts except for variables prefixed or suffixed with pct, which are percentages. The counting unit is the compound-by-compartment-by-biological-group fit, not an individual analytical observation.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"audit_summary"

#' Supplementary table of censoring and compound retention
#'
#' Compact one-row-per-compound summary of analysed records and their censored proportions in sole and prey compartments.
#'
#' @format A data frame with 51 rows and 7 columns.
#' \describe{
#'   \item{\code{Family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{Compound}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{Sole n}}{Number of sole sample-compound observations considered for this compound.}
#'   \item{\code{Sole censored (\%)}}{Percentage of sole observations censored.}
#'   \item{\code{Prey n}}{Number of benthic prey sample-compound observations considered for this compound.}
#'   \item{\code{Prey censored (\%)}}{Percentage of prey observations censored.}
#'   \item{\code{Excluded}}{Selection decision rendered Yes or No according to excluded.}
#' }
#'
#' @details Counts are analytical sample-compound records. The Excluded field reflects the eligibility decision from compound_decisions, including the exploratory retention of HBCDD isomers.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"audit_supp"

#' Predator-to-diet BMF estimates across censoring and diet scenarios
#'
#' Core dataset containing compound-specific field-based normalised predator-to-diet concentration ratios for four censoring treatments and three diet-weighting methods.
#'
#' @format A data frame with 384 rows and 14 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{n_phyla}}{Number of prey phyla contributing chemically determined diet concentrations for the compound (count).}
#'   \item{\code{n_species}}{Total number of chemically characterised prey species in represented phyla for the compound.}
#'   \item{\code{coverage_numeric}}{Sum of the original numerical diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{coverage_biomass}}{Sum of the original biomass diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{c_sole}}{Mean sole-pool concentration in the compound-specific unit field.}
#'   \item{\code{n_sole}}{Number of sole chemical records with finite concentrations used for the mean.}
#'   \item{\code{diet_method}}{Phylum weight model: equal, numeric (prey counts) or biomass (prey count times estimated individual mass).}
#'   \item{\code{c_diet}}{Chemically reconstructed phylum-weighted diet concentration in the compound-specific normalisation unit (not a measured meal concentration).}
#'   \item{\code{bmf}}{Dimensionless operational predator-to-diet normalised concentration ratio, c_sole / c_diet for a positive denominator.}
#'   \item{\code{threshold_above_one}}{TRUE when the computed operational BMF ratio exceeds one.}
#' }
#'
#' @details For each compound and censoring scenario, c_sole is the mean normalised concentration across analysed sole pools. The diet denominator is a mean of chemically analysed prey-species means within each phylum, followed by an average across phyla. Equal weighting averages represented phyla equally; numeric and biomass weighting use reconstructed stomach-content proportions at phylum level, renormalised over phyla with valid chemical data. bmf = c_sole / c_diet when c_diet is positive. Four scenarios are zero, loq_half, loq and mle. Concentrations are ng/g lipid weight for PCB and HBCDD and ng/g dry weight for PFAS. These are operational normalised concentration ratios, not physically reconstructed contaminant concentrations of a consumed meal, transfer efficiencies or probabilities.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_compare"

#' Effects of dietary weighting on compound-specific BMF
#'
#' Wide-format diet-weighting comparisons under each chemical censoring scenario.
#'
#' @format A data frame with 128 rows and 16 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{bmf_equal}}{BMF with equal weights assigned to chemically represented phyla (dimensionless).}
#'   \item{\code{bmf_numeric}}{BMF using relative stomach-prey numerical abundances at phylum level (dimensionless).}
#'   \item{\code{bmf_biomass}}{BMF using reconstructed stomach-prey biomass proportions at phylum level (dimensionless).}
#'   \item{\code{delta_numeric_equal}}{log2(bmf_numeric / bmf_equal), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_equal}}{log2(bmf_biomass / bmf_equal), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_numeric}}{log2(bmf_biomass / bmf_numeric), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_numeric_pct}}{100 * (bmf_biomass / bmf_numeric - 1), percent relative to numeric weighting.}
#'   \item{\code{flip_biomass_numeric}}{TRUE when numeric and biomass BMFs lie strictly on opposite sides of one.}
#'   \item{\code{n_phyla}}{Number of prey phyla contributing chemically determined diet concentrations for the compound (count).}
#'   \item{\code{n_species}}{Total number of chemically characterised prey species in represented phyla for the compound.}
#'   \item{\code{coverage_numeric}}{Sum of the original numerical diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{coverage_biomass}}{Sum of the original biomass diet-weight proportions of chemically represented phyla (0-1).}
#' }
#'
#' @details The three BMF columns share the same compound and censoring treatment. Each delta_* without a pct suffix is the base-2 logarithm of the ratio between two strictly positive BMF values; a percentage column gives the relative change from numeric to biomass weighting. The flip flag is TRUE only when numeric and biomass estimates fall strictly on opposite sides of one; exact equality with one does not count as a flip.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_effects"

#' Family ranges for fully eligible reference BMFs
#'
#' Minimum, median and maximum reference BMF estimates by contaminant family among compounds eligible for all 12 log-based methodological combinations.
#'
#' @format A data frame with 3 rows and 4 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{min_BMF}}{Minimum gamma-MLE biomass-weighted reference BMF across eligible compounds in the family (dimensionless).}
#'   \item{\code{median_BMF}}{Median gamma-MLE biomass-weighted reference BMF across eligible compounds in the family (dimensionless).}
#'   \item{\code{max_BMF}}{Maximum gamma-MLE biomass-weighted reference BMF across eligible compounds in the family (dimensionless).}
#' }
#'
#' @details Eligibility requires 12 available and strictly positive BMF values across four censoring by three diet-weighting scenarios. Summary statistics are then calculated for the gamma-MLE and biomass-weighted reference configuration. Compounds not eligible for the full log-based sensitivity comparison, including compounds with a zero BMF under some scenario, do not contribute to this object. Do not treat it as an unconditional summary of every retained compound.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_family"

#' Species-level sole-to-prey concentration ratio observations
#'
#' One row per compound, scenario and chemically analysed prey species, with the corresponding sole-to-prey concentration ratio and its validity indicator.
#'
#' @format A data frame with 1792 rows and 13 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{species}}{Scientific name/identifier of the chemically analysed organism (not necessarily an observed stomach-prey identification).}
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{c_prey}}{Mean chemically analysed prey-species concentration (compound-specific units).}
#'   \item{\code{n_prey_samples}}{Number of prey chemical records contributing finite concentrations to the species mean.}
#'   \item{\code{c_sole}}{Mean sole-pool concentration in the compound-specific unit field.}
#'   \item{\code{n_sole}}{Number of sole chemical records with finite concentrations used for the mean.}
#'   \item{\code{bmf_prey}}{Operational sole-to-species concentration ratio, c_sole / c_prey for a positive denominator.}
#'   \item{\code{bmf_pair}}{Operational sole-to-species concentration ratio (same value as bmf_prey when valid).}
#'   \item{\code{pair_valid}}{TRUE when the species-level ratio is finite and its denominator is positive.}
#' }
#'
#' @details Each prey species is represented by its mean analytical concentration across samples; the sole numerator is the mean across analysed sole pools. pair_valid identifies ratios with a positive prey denominator and a finite result. A species here is a chemically characterised prey taxon; these rows are not independently observed predator-prey feeding events.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_pairs"

#' Comparison of individual prey and reconstructed diet BMFs
#'
#' Joins diet-based BMF estimates to the distribution of species-level predator-to-prey ratios for each compound and censoring scenario.
#'
#' @format A data frame with 384 rows and 26 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{n_phyla.x}}{Number of chemically represented phyla in the diet-based BMF calculation.}
#'   \item{\code{n_species}}{Total number of chemically characterised prey species represented in the diet denominator.}
#'   \item{\code{coverage_numeric}}{Sum of the original numerical diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{coverage_biomass}}{Sum of the original biomass diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{c_sole}}{Mean sole-pool concentration in the compound-specific unit field.}
#'   \item{\code{n_sole}}{Number of sole chemical records with finite concentrations used for the mean.}
#'   \item{\code{diet_method}}{Phylum weight model: equal, numeric (prey counts) or biomass (prey count times estimated individual mass).}
#'   \item{\code{c_diet}}{Chemically reconstructed phylum-weighted diet concentration in the compound-specific normalisation unit (not a measured meal concentration).}
#'   \item{\code{bmf}}{Dimensionless operational predator-to-diet normalised concentration ratio, c_sole / c_diet for a positive denominator.}
#'   \item{\code{threshold_above_one}}{TRUE when the computed operational BMF ratio exceeds one.}
#'   \item{\code{n_species_measured}}{Number of distinct analysed prey species with finite mean concentrations.}
#'   \item{\code{n_species_valid}}{Number of distinct analysed prey species with valid finite predator-to-prey concentration ratios.}
#'   \item{\code{n_phyla.y}}{Number of phyla with valid positive-denominator prey-species ratios.}
#'   \item{\code{bmf_pair_median}}{Median of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_q25}}{First quartile of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_q75}}{Third quartile of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_min}}{Smallest finite species-level BMF ratio (dimensionless).}
#'   \item{\code{bmf_pair_max}}{Largest finite species-level BMF ratio (dimensionless).}
#'   \item{\code{delta_pair_diet_pct}}{100 * (diet-based bmf / bmf_pair_median - 1), for positive species-ratio median (percent).}
#'   \item{\code{delta_pair_diet_log2}}{log2(diet-based bmf / bmf_pair_median), for positive values (dimensionless).}
#'   \item{\code{diet_position}}{Diet BMF position relative to prey-species Q25/Q75: Below Q25, Within IQR or Above Q75.}
#'   \item{\code{pct_pairs_below_diet}}{Percentage of finite chemically measured prey-species ratios below the reconstructed diet BMF.}
#' }
#'
#' @details The diet BMF is compared with the median and interquartile range of measured-prey species ratios, not with stomach observations paired to individual prey specimens. Suffixes .x and .y preserve the original columns from the diet BMF table and the number of phyla contributing valid species ratios, respectively.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_pair_diet_compare"

#' Distribution summaries of measured-prey species BMFs
#'
#' Per-compound summaries of predator-to-individual-analysed-prey concentration ratios, including median, quartiles, minima and maxima.
#'
#' @format A data frame with 128 rows and 12 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{n_species_measured}}{Number of distinct analysed prey species with finite mean concentrations.}
#'   \item{\code{n_species_valid}}{Number of distinct analysed prey species with valid finite predator-to-prey concentration ratios.}
#'   \item{\code{n_phyla}}{Number of prey phyla contributing chemically determined diet concentrations for the compound (count).}
#'   \item{\code{bmf_pair_median}}{Median of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_q25}}{First quartile of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_q75}}{Third quartile of finite species-level BMF ratios (dimensionless).}
#'   \item{\code{bmf_pair_min}}{Smallest finite species-level BMF ratio (dimensionless).}
#'   \item{\code{bmf_pair_max}}{Largest finite species-level BMF ratio (dimensionless).}
#' }
#'
#' @details Quantiles are computed over valid prey-species ratios with no further weighting for the number of biological samples per prey species; the counts describe chemical coverage and valid ratios.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_pair_summary"

#' Chemically reconstructed prey concentrations by phylum
#'
#' Phylum-level prey-chemistry means and stomach-derived weights used to compute diet-based BMF estimates.
#'
#' @format A data frame with 384 rows and 9 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{c_phylum}}{Arithmetic mean of chemically measured species-level mean prey concentrations within the phylum, in compound-specific units.}
#'   \item{\code{n_species}}{Number of chemically characterised species contributing to the unweighted within-phylum mean.}
#'   \item{\code{p_numeric}}{Proportion of prey counts belonging to the taxon/phylum (0-1), based on stomach contents.}
#'   \item{\code{p_biomass}}{Proportion of reconstructed prey biomass belonging to the taxon/phylum (0-1), among taxa with mass estimates.}
#' }
#'
#' @details c_phylum is the simple arithmetic mean of c_prey across chemically analysed prey species in the phylum. Each measured species contributes equally regardless of its number of chemical pools. Numeric weights count all retained stomach prey taxa; biomass weights use only taxa with an available positive individual mass. Weights are renormalised across chemically represented phyla for each compound when constructing c_diet.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_phylum"

#' Compound concentrations and ratios for chemically analysed prey species
#'
#' Chemical means by prey species, paired with the corresponding mean sole concentration and operational concentration ratio.
#'
#' @format A data frame with 1792 rows and 11 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{species}}{Scientific name/identifier of the chemically analysed organism (not necessarily an observed stomach-prey identification).}
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{c_prey}}{Mean chemically analysed prey-species concentration (compound-specific units).}
#'   \item{\code{n_prey_samples}}{Number of prey chemical records contributing finite concentrations to the species mean.}
#'   \item{\code{c_sole}}{Mean sole-pool concentration in the compound-specific unit field.}
#'   \item{\code{n_sole}}{Number of sole chemical records with finite concentrations used for the mean.}
#'   \item{\code{bmf_prey}}{Operational sole-to-species concentration ratio, c_sole / c_prey for a positive denominator.}
#' }
#'
#' @details Concentrations are sample-level analytical means within the same chemical scenario and family-specific normalisation; individual prey species are not necessarily identified among the actual stomach contents.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_prey"

#' BMF sensitivity to bivalve siphon mass correction
#'
#' Recomputed reference-scenario biomass-weighted BMFs for alternative multipliers applied to bivalve individual mass estimates.
#'
#' @format A data frame with 256 rows and 17 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{n_phyla}}{Number of prey phyla contributing chemically determined diet concentrations for the compound (count).}
#'   \item{\code{n_species}}{Number of species or references represented in the dataset-specific grouping (count; see Details).}
#'   \item{\code{coverage_numeric}}{Sum of the original numerical diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{coverage_biomass}}{Sum of the original biomass diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{c_sole}}{Mean sole-pool concentration in the compound-specific unit field.}
#'   \item{\code{n_sole}}{Number of sole chemical records with finite concentrations used for the mean.}
#'   \item{\code{diet_method}}{Phylum weight model: equal, numeric (prey counts) or biomass (prey count times estimated individual mass).}
#'   \item{\code{c_diet}}{Chemically reconstructed phylum-weighted diet concentration in the compound-specific normalisation unit (not a measured meal concentration).}
#'   \item{\code{bmf}}{Dimensionless operational predator-to-diet normalised concentration ratio, c_sole / c_diet for a positive denominator.}
#'   \item{\code{threshold_above_one}}{TRUE when the computed operational BMF ratio exceeds one.}
#'   \item{\code{factor_siphon}}{Dimensionless multiplier applied to the individual AFDW mass of Bivalvia; reference value is 0.0375.}
#'   \item{\code{bmf_reference}}{BMF under the reference bivalve siphon factor 0.0375 (dimensionless).}
#'   \item{\code{delta_pct}}{100 * (bmf / bmf_reference - 1), percent relative to the reference factor.}
#' }
#'
#' @details The factor_siphon range covers 0.02 to 0.05, with a 0.0375 reference. The factor is applied to bivalves (Class Bivalvia), not indiscriminately to all molluscs. delta_pct = 100 * (bmf / bmf_reference - 1) for positive reference BMF. These are deterministic sensitivity outcomes rather than confidence intervals.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_siphon"

#' BMF classification stability across twelve methodological combinations
#'
#' For each compound, counts of BMF estimates above, below or exactly at one, with minimum, maximum and a method-dependence classification.
#'
#' @format A data frame with 32 rows and 9 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{n_above_1}}{Number of BMF method combinations with a ratio strictly above one.}
#'   \item{\code{n_below_1}}{Number of BMF method combinations with a ratio strictly below one.}
#'   \item{\code{n_equal_1}}{Number of BMF method combinations with a ratio exactly equal to one.}
#'   \item{\code{bmf_min}}{Smallest finite BMF across the scenarios represented (dimensionless).}
#'   \item{\code{bmf_max}}{Largest finite BMF across the scenarios represented (dimensionless).}
#'   \item{\code{flip}}{TRUE when at least one estimate lies below and one lies above one.}
#'   \item{\code{classification}}{Interpretation of the position of all 12 scenario BMFs relative to one and of completeness.}
#' }
#'
#' @details The 12 combinations are four treatments of left-censored concentrations times three diet-weighting methods. The output omits scenario-level rows and classifies results as Incomplete, Method-dependent, Always > 1, Always < 1, Always = 1, or Includes BMF = 1. A value above one is an operational concentration comparison, not proof of causal trophic transfer.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"bmf_threshold"

#' Current chemical compound eligibility and censoring decisions
#'
#' Selection table for PCB congeners, HBCDD isomers and PFAS compounds based on biotic coverage and family-specific censoring rules.
#'
#' @format A data frame with 51 rows and 10 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{n_compartments}}{Number of represented biological compartments (expected maximum: two).}
#'   \item{\code{high_censoring}}{TRUE when at least one compartment has a censoring percentage at or above the configured 80-percent threshold.}
#'   \item{\code{all_censored_in_compartment}}{TRUE when at least one represented compartment has no quantified observations for the compound.}
#'   \item{\code{censored_in_sole_pct}}{Percentage censored among age-zero sole chemical pools.}
#'   \item{\code{censored_in_prey_pct}}{Percentage censored among benthic prey chemical samples.}
#'   \item{\code{excluded}}{TRUE when the compound is excluded from current BMF estimates under the selection rules.}
#'   \item{\code{exclusion_reason}}{Human-readable reason for exclusion or retention, including the HBCDD exception.}
#'   \item{\code{evidence_tier}}{Evidence designation: excluded, standard, or exploratory HBCDD retention due to censoring.}
#' }
#'
#' @details Exclusion requires both biological compartments; for PCB/PFAS, any compartment with at least 80 percent censored observations triggers exclusion. The HBCDD isomers are retained irrespective of that threshold where both compartments are represented. A compound not present in this table has no current decision (not assessed); its absence must not be interpreted as retained or never quantified. The table reflects only chemicals processed in the current analysis.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"compound_decisions"

#' Coverage of stomach prey by chemically characterised taxa and phylum
#'
#' Abundance- and biomass-based coverage of stomach-content taxonomic identifications by measured chemical prey taxa.
#'
#' @format A data frame with 6 rows and 8 columns.
#' \describe{
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{match_class}}{Type of stomach-prey versus chemical-prey taxon match: Direct, Genus-level proxy or No specific match.}
#'   \item{\code{n_prey}}{Total counted stomach prey individuals across retained age-zero soles, or an aggregate of those counts.}
#'   \item{\code{biomass}}{Summed share of reconstructed dietary prey biomass assigned to this group (fraction, 0-1; not mg).}
#'   \item{\code{numeric_pct_within_phylum}}{Percent of individuals in this match class within the same phylum.}
#'   \item{\code{biomass_pct_within_phylum}}{Percent of estimated phylum biomass in this match class.}
#'   \item{\code{numeric_pct_total}}{Percent of all counted stomach prey individuals belonging to this phylum/match class.}
#'   \item{\code{biomass_pct_total}}{Percent of the total reconstructed prey biomass assigned to this phylum/match class.}
#' }
#'
#' @details Match classes are Direct, Genus-level proxy and No specific match. The biomass column is a proportion of the total reference diet, despite its short name: it is not biomass in mg. Within-phylum percentages are calculated using the sum within each phylum; total percentages are calculated across all phyla. These matches concern chemical coverage rather than uncertain identification of stomach prey.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"coverage_phylum"

#' Diet proportions under bivalve siphon mass assumptions
#'
#' Phylum-level numerical and biomass proportions repeatedly reconstructed at each bivalve siphon multiplier.
#'
#' @format A data frame with 32 rows and 8 columns.
#' \describe{
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{p_numeric}}{Numerical proportion of all retained stomach prey counts assigned to the phylum (0-1).}
#'   \item{\code{p_biomass}}{Biomass proportion of all prey taxa with valid mass estimates at this siphon factor (0-1).}
#'   \item{\code{p_biomass_q25}}{Reconstructed prey biomass proportion (0-1) using the simultaneous Q25 mass scenario.}
#'   \item{\code{p_biomass_q75}}{Reconstructed prey biomass proportion (0-1) using the simultaneous Q75 mass scenario.}
#'   \item{\code{n_prey}}{Total counted stomach prey individuals across retained age-zero soles, or an aggregate of those counts.}
#'   \item{\code{n_prey_with_mass}}{Number of counted prey individuals from taxa with usable positive mass assumptions.}
#'   \item{\code{factor_siphon}}{Dimensionless multiplier applied to the individual AFDW mass of Bivalvia; reference value is 0.0375.}
#' }
#'
#' @details The bivalve-specific correction factor varies between 0.02 and 0.05 (including reference 0.0375). All prey observations contribute to numerical proportions, but only prey taxa with all three usable positive mass assumptions contribute to biomass indices. Q25, median and Q75 masses are propagated simultaneously across taxa and are deterministic alternatives, not estimated confidence bounds.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"diet_sensitivity"

#' Family-level summary of log-based diet and censoring sensitivities
#'
#' Chemical-family medians of log2 sensitivity to dietary weighting and censored-data treatment, with counts of which component dominates.
#'
#' @format A data frame with 3 rows and 7 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{n_compounds}}{Number of compounds included in this aggregation (count).}
#'   \item{\code{median_S_diet}}{Median across eligible compounds of the log2 dietary-weighting sensitivity (dimensionless).}
#'   \item{\code{median_S_censor}}{Median across eligible compounds of the log2 censored-data-treatment sensitivity (dimensionless).}
#'   \item{\code{median_S_numeric_biomass}}{Median across eligible compounds of log2 sensitivity comparing only numeric and biomass diet weights.}
#'   \item{\code{n_diet_dominant}}{Number of compounds whose log2 dietary-weighting sensitivity exceeds censoring sensitivity.}
#'   \item{\code{n_censor_dominant}}{Number of compounds whose log2 censoring sensitivity exceeds dietary-weighting sensitivity.}
#' }
#'
#' @details Based only on compounds with strictly positive, finite BMF estimates in all 12 methodological combinations. The sensitivity measures are log2 fold contrasts, dimensionless; a zero-BMF scenario is ineligible for this log-based hierarchy. Sensitivity to numeric versus biomass weighting excludes the equal-weighting method.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"family_sensitivity"

#' Family comparisons of numeric and biomass weighting effects
#'
#' Within each censoring scenario and chemical family, summarises the magnitude and frequency of changes between numeric and biomass-based BMF estimates.
#'
#' @format A data frame with 12 rows and 5 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{n_compounds}}{Number of compounds included in this aggregation (count).}
#'   \item{\code{median_delta}}{Median log2(bmf_biomass / bmf_numeric) among eligible positive ratios.}
#'   \item{\code{n_flip}}{Number of compounds whose numeric versus biomass weights switch the strict BMF > 1 classification.}
#' }
#'
#' @details median_delta is the median log2(bmf_biomass / bmf_numeric) over compounds where both ratios are positive. n_flip counts compounds crossing one strictly between numeric and biomass weighting.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"family_summary"

#' Gamma MLE fitting and fallback audit by analytical group
#'
#' Fitting diagnostics for the gamma left-censoring model applied within each retained compound, compartment and biological group.
#'
#' @format A data frame with 480 rows and 21 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{compartment}}{Biological compartment, Sole or Prey.}
#'   \item{\code{group}}{Chemical fitting stratum: sole group or chemically analysed prey species.}
#'   \item{\code{distribution}}{Fitted concentration distribution (gamma in the packaged fit_audit object).}
#'   \item{\code{n}}{Number of analytical sample-compound records used for a compound-by-compartment-by-group gamma fit.}
#'   \item{\code{n_censored}}{Number of censored analytical records.}
#'   \item{\code{n_uncensored}}{Number of uncensored analytical records.}
#'   \item{\code{pct_censored}}{Percentage of censored analytical records, 100 times n_censored divided by n.}
#'   \item{\code{fit_status}}{Status of likelihood optimisation: fitted, insufficient_data, failed or excluded, as appropriate.}
#'   \item{\code{fallback_method}}{Alternative imputation used when maximum-likelihood fitting is not available (LOQ/2).}
#'   \item{\code{fit_error}}{Diagnostic message explaining a failed/infeasible fit when applicable.}
#'   \item{\code{parameter_1}}{First distribution parameter: gamma shape for the gamma audit (dimensionless).}
#'   \item{\code{parameter_2}}{Second distribution parameter: gamma scale for the gamma audit (same units as concentration).}
#'   \item{\code{n_imputed}}{Number of censored records requiring imputation in the fitting group.}
#'   \item{\code{n_imputed_mle}}{Number of censored records imputed using successful gamma maximum-likelihood fitting.}
#'   \item{\code{n_imputed_fallback}}{Number of censored records assigned a fallback value (LOQ/2).}
#'   \item{\code{n_unresolved}}{Number of records with non-finite imputed concentrations.}
#'   \item{\code{imputed_min}}{Minimum imputed censored concentration (same concentration units as the compound).}
#'   \item{\code{imputed_median}}{Median imputed censored concentration (same concentration units as the compound).}
#'   \item{\code{imputed_max}}{Maximum imputed censored concentration (same concentration units as the compound).}
#' }
#'
#' @details A gamma distribution is fitted by censored-data maximum likelihood where there are sufficient observations and uncensored values. Shape and scale are gamma parameters, not regression coefficients. For unsuccessful fits, censored observations receive LOQ/2. Parameter_1 is gamma shape (dimensionless); parameter_2 is gamma scale with the same concentration units as the compound. imputed_min, imputed_median and imputed_max refer to imputed censored concentrations, not uncertainty limits.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"fit_audit"

#' Individual prey-mass reference table after siphon adjustment
#'
#' Best available individual prey ash-free dry mass estimates matched taxonomically to stomach-content prey taxa.
#'
#' @format A data frame with 26 rows and 16 columns.
#' \describe{
#'   \item{\code{observed_species}}{Scientific taxon identified in the juvenile-sole stomach-content dataset.}
#'   \item{\code{Class}}{Taxonomic class of the stomach prey; used particularly to identify Bivalvia.}
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{n_species}}{Number of eligible mass-reference species matched to this stomach taxon.}
#'   \item{\code{mass_median}}{Reference-individual median mass (mg AFDW), multiplied by factor_siphon ONLY for Bivalvia.}
#'   \item{\code{mass_minq25}}{Lowest Q25 mass among matched reference species (mg AFDW), multiplied by factor_siphon ONLY for Bivalvia.}
#'   \item{\code{mass_maxq75}}{Highest Q75 mass among matched reference species (mg AFDW), multiplied by factor_siphon ONLY for Bivalvia.}
#'   \item{\code{source}}{Name(s) of the literature or package source(s) used for individual prey masses.}
#'   \item{\code{mass_iqr}}{mass_maxq75 - mass_minq25 BEFORE the bivalve correction (mg AFDW), because these quality fields were joined from mass_stats_matches.}
#'   \item{\code{mass_iqr_rel}}{mass_iqr / mass_median for taxonomic matches based on at least two candidate reference species (dimensionless).}
#'   \item{\code{mass_match_quality}}{Categorical descriptor of mass-evidence strength and between-candidate spread.}
#'   \item{\code{shared_rank}}{Deepest shared taxonomic rank between stomach prey and mass reference candidates.}
#'   \item{\code{priority}}{Ordered index of that shared taxonomic rank, larger indices representing more specific levels.}
#'   \item{\code{n_candidates}}{Number of candidate reference taxa at the selected shared taxonomic rank.}
#'   \item{\code{taxonomic_matching_quality}}{Matching-quality category based on shared taxonomic rank (Perfect, High, Moderate, Low, or missing).}
#'   \item{\code{factor_siphon}}{Dimensionless multiplier applied to the individual AFDW mass of Bivalvia; reference value is 0.0375.}
#' }
#'
#' @details Individual masses are expressed as mg ash-free dry weight (AFDW). Reference estimates originate from Btrait and the Gogina et al. source workbook and may use taxonomic substitutes. mass_median is the median of reference-species medians; mass_minq25 is the minimum reference-species Q25 and mass_maxq75 is the maximum reference-species Q75, rather than a statistical confidence interval. For Class Bivalvia only, all three mass estimates are multiplied by the factor_siphon reference value of 0.0375. The mass spread and matching-quality fields are descriptive diagnostics, not statistical standard errors.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"mass_chopin"

#' Stomach prey counts, imputed masses and reconstructed contributions
#'
#' Prey-taxon-level summary of occurrence counts, literature-derived individual masses, inferred prey biomass and proportions of the reconstructed juvenile-sole diet.
#'
#' @format A data frame with 26 rows and 19 columns.
#' \describe{
#'   \item{\code{observed_species}}{Scientific taxon identified in the juvenile-sole stomach-content dataset.}
#'   \item{\code{n_prey}}{Total counted stomach prey individuals across retained age-zero soles, or an aggregate of those counts.}
#'   \item{\code{n_fish}}{Number of distinct sole stomach samples in which the prey taxon was recorded.}
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{mass_median}}{Reference individual prey mass estimate (mg ash-free dry weight, AFDW): median candidate median, corrected for bivalve siphon when applicable.}
#'   \item{\code{mass_q25}}{Lower alternative individual mass estimate (mg AFDW), based on the candidate-reference Q25 envelope.}
#'   \item{\code{mass_q75}}{Upper alternative individual mass estimate (mg AFDW), based on the candidate-reference Q75 envelope.}
#'   \item{\code{mass_available}}{TRUE when median, Q25 and Q75 individual masses are all finite and positive.}
#'   \item{\code{biomass_index}}{n_prey * mass_median (mg AFDW), if all three reference masses are available.}
#'   \item{\code{biomass_index_q25}}{n_prey * mass_q25 (mg AFDW), if mass_available.}
#'   \item{\code{biomass_index_q75}}{n_prey * mass_q75 (mg AFDW), if mass_available.}
#'   \item{\code{p_numeric}}{Fraction of ALL stomach prey counts belonging to this prey taxon (0-1).}
#'   \item{\code{p_biomass}}{Taxon biomass fraction using median mass among valid-mass prey taxa (0-1).}
#'   \item{\code{p_biomass_q25}}{Taxon biomass fraction using Q25 masses simultaneously (0-1).}
#'   \item{\code{p_biomass_q75}}{Taxon biomass fraction using Q75 masses simultaneously (0-1).}
#'   \item{\code{factor_siphon}}{Dimensionless multiplier applied to the individual AFDW mass of Bivalvia; reference value is 0.0375.}
#'   \item{\code{source}}{Name(s) of the literature or package source(s) used for individual prey masses.}
#'   \item{\code{mass_match_quality}}{Categorical descriptor of mass-evidence strength and between-candidate spread.}
#'   \item{\code{taxonomic_matching_quality}}{Matching-quality category based on shared taxonomic rank (Perfect, High, Moderate, Low, or missing).}
#' }
#'
#' @details n_prey is the sum of stomach N_tract prey counts across retained autumn G0 fish; n_fish is the number of distinct sampled soles containing the taxon. Individual mass assumptions are mg AFDW; the biomass indices multiply mass by prey count (mg AFDW). The associated p_biomass variables are proportions, not physical concentrations. Biomass denominators include only taxa with positive finite mass values for all three scenarios, whereas numeric proportions include all retained stomach taxa.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"mass_review"

#' Aggregated candidate prey-mass references by taxonomic match
#'
#' Reference-mass lookup summary for each stomach-content taxon, before applying the bivalve siphon correction.
#'
#' @format A data frame with 62 rows and 13 columns.
#' \describe{
#'   \item{\code{observed_species}}{Scientific taxon identified in the juvenile-sole stomach-content dataset.}
#'   \item{\code{n_species}}{Number of candidate reference species selected at the most specific available shared taxonomic rank.}
#'   \item{\code{mass_median}}{Median of the candidate reference-species median individual masses (mg AFDW, before siphon correction).}
#'   \item{\code{mass_minq25}}{Minimum of the candidate reference-species Q25 masses (mg AFDW, before siphon correction).}
#'   \item{\code{mass_maxq75}}{Maximum of the candidate reference-species Q75 masses (mg AFDW, before siphon correction).}
#'   \item{\code{source}}{Name(s) of the literature or package source(s) used for individual prey masses.}
#'   \item{\code{mass_iqr}}{Difference mass_maxq75 - mass_minq25 (mg AFDW); cross-reference spread, not necessarily an empirical IQR.}
#'   \item{\code{mass_iqr_rel}}{mass_iqr / mass_median for taxonomic matches based on at least two candidate reference species (dimensionless).}
#'   \item{\code{mass_match_quality}}{Categorical descriptor of mass-evidence strength and between-candidate spread.}
#'   \item{\code{shared_rank}}{Deepest shared taxonomic rank between stomach prey and mass reference candidates.}
#'   \item{\code{priority}}{Ordered index of that shared taxonomic rank, larger indices representing more specific levels.}
#'   \item{\code{n_candidates}}{Number of candidate reference taxa at the selected shared taxonomic rank.}
#'   \item{\code{taxonomic_matching_quality}}{Matching-quality category based on shared taxonomic rank (Perfect, High, Moderate, Low, or missing).}
#' }
#'
#' @details The best matching taxonomic rank is chosen first. The median is formed from candidate reference-species medians, the lower mass from the minimum candidate Q25, and the upper mass from the maximum candidate Q75. mass_iqr is defined operationally as mass_maxq75 - mass_minq25 and mass_iqr_rel as that difference divided by mass_median when at least two reference species exist. Despite the name, this need not be a classical interquartile range of a single distribution. The quality classifications describe reference matching and spread, not uncertainty intervals or validation against measured stomach-prey masses.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"mass_stats_matches"

#' Prepared benthic prey chemical observations
#'
#' Long-format benthic prey analytical concentrations, measurement censoring indicators and detection limits from the biotic source dataset.
#'
#' @format A data frame with 2225 rows and 15 columns.
#' \describe{
#'   \item{\code{sample_id}}{Biological sample or pool identifier from the source chemical dataset.}
#'   \item{\code{compartment}}{Biological compartment, Sole or Prey.}
#'   \item{\code{group}}{Chemical fitting stratum: sole group or chemically analysed prey species.}
#'   \item{\code{species}}{Scientific name/identifier of the chemically analysed organism (not necessarily an observed stomach-prey identification).}
#'   \item{\code{taxa}}{Broad source taxonomy grouping, used to assign prey phyla (Polychaetes, Crustaceans or Bivalves).}
#'   \item{\code{feeding_mode}}{Source ecological feeding-mode category when available; may be missing for soles.}
#'   \item{\code{dry_fraction}}{Dry matter percentage from the source biometrics (source field named *_percent); recorded as provided, not assumed a 0-to-1 fraction.}
#'   \item{\code{lipid_fraction}}{Total lipid content in percent dry weight from the source field lip_dw_percent; despite the name, not assumed a 0-to-1 fraction.}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound}}{Compound identifier combining chemical name and analytical mass-normalisation suffix.}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{concentration}}{Measured or source-prepared concentration in the unit field; censored records follow the family-specific source encoding.}
#'   \item{\code{censored}}{Logical flag: TRUE if the analytical observation is left-censored at its LOQ.}
#'   \item{\code{loq}}{Limit of quantification (same units as concentration), recorded as 0 for uncensored PCB and HBCDD in the preparation.}
#' }
#'
#' @details Each row describes a sample-compound observation for a chemically analysed prey organism. Original family-specific mass normalisations are retained: lipid weight for PCB/HBCDD and dry weight for PFAS. For HBCDD, an original value is treated as an LOQ when its censored companion is zero; for PFAS, a zero censored companion indicates non-quantification and the supplied reference LOQ is used; PCBs are treated as uncensored in this pipeline. This object is prepared before current BMF compound exclusions and does not imply that every chemical will be retained.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"prey"

#' Numerical and biomass diet proportions by prey phylum
#'
#' Aggregated stomach-derived prey abundance and estimated biomass weights for phylum-level dietary reconstruction.
#'
#' @format A data frame with 4 rows and 7 columns.
#' \describe{
#'   \item{\code{Phylum}}{Prey phylum (typically Annelida, Arthropoda or Mollusca).}
#'   \item{\code{p_biomass_phylum}}{Biomass proportion of the total reconstructed diet assigned to the phylum using median masses (0-1).}
#'   \item{\code{p_biomass_q25_phylum}}{Biomass proportion using all taxa Q25 mass assumptions simultaneously (0-1).}
#'   \item{\code{p_biomass_q75_phylum}}{Biomass proportion using all taxa Q75 mass assumptions simultaneously (0-1).}
#'   \item{\code{p_numeric_phylum}}{Numerical proportion of all stomach prey counts assigned to the phylum (0-1).}
#'   \item{\code{n_prey}}{Total counted stomach prey individuals across retained age-zero soles, or an aggregate of those counts.}
#'   \item{\code{n_prey_with_mass}}{Number of counted prey individuals from taxa with usable positive mass assumptions.}
#' }
#'
#' @details Numerical proportions use all retained stomach-count observations: number in phylum divided by total prey count. Biomass proportions use summed n_prey times estimated individual mass (mg AFDW), divided by the summed valid biomass indices across prey taxa. Q25 and Q75 use the corresponding mass assumptions for all taxa simultaneously. These alternatives do not describe probabilities or confidence intervals, and their effects on BMF need not be ordered.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"prop_diet"

#' Prepared concentrations under four censored-data scenarios
#'
#' Sample-by-compound analytical data after compound selection and imputation under zero, LOQ/2, LOQ and gamma-MLE scenarios.
#'
#' @format A data frame with 1690 rows and 20 columns.
#' \describe{
#'   \item{\code{sample_id}}{Biological sample or pool identifier from the source chemical dataset.}
#'   \item{\code{compartment}}{Biological compartment, Sole or Prey.}
#'   \item{\code{group}}{Chemical fitting stratum: sole group or chemically analysed prey species.}
#'   \item{\code{species}}{Scientific name/identifier of the chemically analysed organism (not necessarily an observed stomach-prey identification).}
#'   \item{\code{taxa}}{Broad source taxonomy grouping, used to assign prey phyla (Polychaetes, Crustaceans or Bivalves).}
#'   \item{\code{feeding_mode}}{Source ecological feeding-mode category when available; may be missing for soles.}
#'   \item{\code{dry_fraction}}{Dry matter percentage from the source biometrics (source field named *_percent); recorded as provided, not assumed a 0-to-1 fraction.}
#'   \item{\code{lipid_fraction}}{Total lipid content in percent dry weight from the source field lip_dw_percent; despite the name, not assumed a 0-to-1 fraction.}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound}}{Compound identifier combining chemical name and analytical mass-normalisation suffix.}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{concentration}}{Measured or source-prepared concentration in the unit field; censored records follow the family-specific source encoding.}
#'   \item{\code{censored}}{Logical flag: TRUE if the analytical observation is left-censored at its LOQ.}
#'   \item{\code{loq}}{Limit of quantification (same units as concentration), recorded as 0 for uncensored PCB and HBCDD in the preparation.}
#'   \item{\code{imputation_method}}{Treatment of censored records under gamma fitting: observed, mle or loq_half_fallback.}
#'   \item{\code{value_zero}}{Concentration after replacing censored observations by zero (same units as concentration).}
#'   \item{\code{value_loq_half}}{Concentration after replacing censored observations by LOQ/2 (same units as concentration).}
#'   \item{\code{value_loq}}{Concentration after replacing censored observations by LOQ (same units as concentration).}
#'   \item{\code{value_mle}}{Concentration after gamma conditional-expectation imputation, using LOQ/2 when fitting fails (same units as concentration).}
#' }
#'
#' @details Each original record is retained with its censor indicator and original/LOQ concentration. value_zero substitutes zero for censored observations; value_loq_half substitutes LOQ/2; value_loq substitutes LOQ; value_mle uses the gamma censored-data conditional expectation or LOQ/2 after failed or infeasible fitting. All four value columns use the same unit as the original concentration. For noncensored records all scenario values equal the observed concentration. This is the principal chemistry input to the BMF pipeline.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"scenarios"

#' Joined data for the diet-weighting versus censoring sensitivity analysis
#'
#' Per-scenario compound-level BMF weighting effects combined with chemical censoring percentages in sole and prey samples.
#'
#' @format A data frame with 128 rows and 19 columns.
#' \describe{
#'   \item{\code{scenario}}{Treatment of censored concentrations: zero, loq_half, loq or mle (gamma maximum likelihood with possible LOQ/2 fallback).}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{bmf_equal}}{BMF with equal weights assigned to chemically represented phyla (dimensionless).}
#'   \item{\code{bmf_numeric}}{BMF using relative stomach-prey numerical abundances at phylum level (dimensionless).}
#'   \item{\code{bmf_biomass}}{BMF using reconstructed stomach-prey biomass proportions at phylum level (dimensionless).}
#'   \item{\code{delta_numeric_equal}}{log2(bmf_numeric / bmf_equal), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_equal}}{log2(bmf_biomass / bmf_equal), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_numeric}}{log2(bmf_biomass / bmf_numeric), defined for positive BMF values (dimensionless).}
#'   \item{\code{delta_biomass_numeric_pct}}{100 * (bmf_biomass / bmf_numeric - 1), percent relative to numeric weighting.}
#'   \item{\code{flip_biomass_numeric}}{TRUE when numeric and biomass BMFs lie strictly on opposite sides of one.}
#'   \item{\code{n_phyla}}{Number of prey phyla contributing chemically determined diet concentrations for the compound (count).}
#'   \item{\code{n_species}}{Number of species or references represented in the dataset-specific grouping (count; see Details).}
#'   \item{\code{coverage_numeric}}{Sum of the original numerical diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{coverage_biomass}}{Sum of the original biomass diet-weight proportions of chemically represented phyla (0-1).}
#'   \item{\code{prey}}{Percentage censored in the chemically analysed prey compartment.}
#'   \item{\code{sole}}{Percentage censored in the analysed sole compartment.}
#'   \item{\code{censor_max}}{Largest censored-observation percentage among sole and prey compartments (percentage points).}
#' }
#'
#' @details The three weighting BMFs and their log2 contrasts are derived from bmf_compare. The prey and sole columns are percentages of censored observations; censor_max is their maximum. The table includes scenarios that may produce zero BMFs; those results are subsequently excluded from strictly positive log2-based sensitivity summaries.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"sensitivity_data"

#' Hierarchical comparison of diet- and censoring-driven BMF sensitivity
#'
#' Per-compound dimensionless log2 sensitivity summaries that compare variation across dietary weighting against variation across left-censoring treatments.
#'
#' @format A data frame with 31 rows and 13 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{S_diet}}{Median across censoring treatments of log2(max/min BMF) across diet weighting methods, for positive BMFs.}
#'   \item{\code{S_diet_max}}{Maximum across censoring treatments of the within-treatment log2 dietary sensitivity.}
#'   \item{\code{S_censor}}{Median across diet weightings of log2(max/min BMF) across censoring treatments, for positive BMFs.}
#'   \item{\code{S_censor_max}}{Maximum across diet weightings of the within-weighting log2 censoring sensitivity.}
#'   \item{\code{censor_max}}{Largest censored-observation percentage among sole and prey compartments (percentage points).}
#'   \item{\code{dominant_source}}{Whether diet weighting, censoring treatment, or neither produced the larger median log2 sensitivity.}
#'   \item{\code{fold_diet}}{2 raised to S_diet, a multiplicative fold range (dimensionless).}
#'   \item{\code{fold_censor}}{2 raised to S_censor, a multiplicative fold range (dimensionless).}
#'   \item{\code{delta_S}}{S_diet - S_censor (log2 units).}
#'   \item{\code{S_numeric_biomass}}{Median across censoring scenarios of log2(max/min BMF) comparing numeric versus biomass weighting only.}
#' }
#'
#' @details Only compounds with 12 finite strictly positive BMF estimates across four censoring treatments and three dietary methods are included. For each censoring scenario, S_diet is log2(max BMF/min BMF) across equal, numeric and biomass methods, then summarised as its median across censoring treatments; S_diet_max is the maximum across these treatments. Similarly, S_censor is the median across dietary methods of log2(max BMF/min BMF) across censoring treatments; S_censor_max is the corresponding maximum. delta_S = S_diet - S_censor. fold_diet = 2^S_diet, fold_censor = 2^S_censor. S_numeric_biomass uses only those two dietary methods.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"sensitivity_hierarchy"

#' Family-level summary of sensitivity to the bivalve siphon factor
#'
#' Medians of compound-level fold ranges, log2 ranges and largest absolute percentage changes under alternative bivalve siphon corrections.
#'
#' @format A data frame with 3 rows and 5 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{n_compounds}}{Number of compounds included in this aggregation (count).}
#'   \item{\code{median_fold}}{Median compound-level fold_range across the family.}
#'   \item{\code{median_S}}{Median compound-level S_siphon across the family (log2 units).}
#'   \item{\code{median_max_delta_pct}}{Median compound-level largest absolute percentage difference from the reference siphon correction.}
#' }
#'
#' @details All summaries describe deterministic variation across a finite grid of correction factors, not statistical variability. Fold and log2 measures are dimensionless; median_max_delta_pct is a percentage.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"siphon_family"

#' Compound-level sensitivity to bivalve siphon mass correction
#'
#' Compound-specific extrema and fold ranges of the biomass-weighted gamma-MLE BMF over alternative bivalve siphon correction factors.
#'
#' @format A data frame with 32 rows and 8 columns.
#' \describe{
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{bmf_min}}{Smallest finite BMF across the scenarios represented (dimensionless).}
#'   \item{\code{bmf_max}}{Largest finite BMF across the scenarios represented (dimensionless).}
#'   \item{\code{max_abs_delta_pct}}{Largest absolute percent difference in BMF relative to the reference siphon correction.}
#'   \item{\code{fold_range}}{Largest BMF divided by smallest strictly positive BMF in the scenario grid (dimensionless).}
#'   \item{\code{S_siphon}}{log2(fold_range), dimensionless sensitivity to the siphon correction factor.}
#' }
#'
#' @details Minimum and maximum are observed over the factor grid (0.02 to 0.05, including 0.0375), not a confidence or credible interval. fold_range = bmf_max / bmf_min where the denominator is positive; S_siphon = log2(fold_range). max_abs_delta_pct is the largest absolute percentage change from the reference factor.
#'
#' @source CHOPIN pipeline: \code{data-raw/03_CHOPIN_BMF_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"siphon_summary"

#' Prepared age-zero sole chemical observations
#'
#' Long-format concentrations in chemical pools of age-zero Solea solea selected for the configured season, before filtering chemical compounds for BMF eligibility.
#'
#' @format A data frame with 510 rows and 15 columns.
#' \describe{
#'   \item{\code{sample_id}}{Biological sample or pool identifier from the source chemical dataset.}
#'   \item{\code{compartment}}{Biological compartment, Sole or Prey.}
#'   \item{\code{group}}{Chemical fitting stratum: sole group or chemically analysed prey species.}
#'   \item{\code{species}}{Scientific name/identifier of the chemically analysed organism (not necessarily an observed stomach-prey identification).}
#'   \item{\code{taxa}}{Broad source taxonomy grouping, used to assign prey phyla (Polychaetes, Crustaceans or Bivalves).}
#'   \item{\code{feeding_mode}}{Source ecological feeding-mode category when available; may be missing for soles.}
#'   \item{\code{dry_fraction}}{Dry matter percentage from the source biometrics (source field named *_percent); recorded as provided, not assumed a 0-to-1 fraction.}
#'   \item{\code{lipid_fraction}}{Total lipid content in percent dry weight from the source field lip_dw_percent; despite the name, not assumed a 0-to-1 fraction.}
#'   \item{\code{family}}{Chemical family (PCB, HBCDD or PFAS).}
#'   \item{\code{compound}}{Compound identifier combining chemical name and analytical mass-normalisation suffix.}
#'   \item{\code{compound_name}}{Chemical identifier as used in the CHOPIN analytical pipeline.}
#'   \item{\code{unit}}{Concentration mass-normalisation: ng_glw means ng per g total lipid weight for PCB and HBCDD; ng_gdw means ng per g dry weight for PFAS. The audit table may instead spell out the full unit.}
#'   \item{\code{concentration}}{Measured or source-prepared concentration in the unit field; censored records follow the family-specific source encoding.}
#'   \item{\code{censored}}{Logical flag: TRUE if the analytical observation is left-censored at its LOQ.}
#'   \item{\code{loq}}{Limit of quantification (same units as concentration), recorded as 0 for uncensored PCB and HBCDD in the preparation.}
#' }
#'
#' @details One row per sole pool and analyte. Columns retain sample ID, organism, chemical family, mass normalisation, censoring and LOQ status. PCB and HBCDD concentrations are ng/g lipid weight; PFAS are ng/g dry weight. Dry and lipid fields originate from the percent-based source variables dry_content_percent and lip_dw_percent, despite having fraction in their column names. These are chemical sample pools, not individual fish.
#'
#' @source CHOPIN pipeline: \code{data-raw/01_CHOPIN_contamination_gamma.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"sole"

#' Summary of autumn age-zero sole-pool biometrics
#'
#' Single-row descriptive summary of total length and wet mass in autumn G0 chemical sole pools.
#'
#' @format A data frame with 1 rows and 6 columns.
#' \describe{
#'   \item{\code{median_length}}{Median total length across autumn age-zero sole samples (cm).}
#'   \item{\code{min_length_68}}{Minimum over pools of length_TL_cm - length_TL_cm_sd (cm); descriptive source-value envelope, not a CI.}
#'   \item{\code{max_length_68}}{Maximum over pools of length_TL_cm - length_TL_cm_sd (cm).}
#'   \item{\code{median_ww}}{Median total wet mass across autumn age-zero sole pools (g wet weight).}
#'   \item{\code{min_ww_68}}{Minimum over pools of mass_tot_gww - mass_tot_gww_sd (g wet weight).}
#'   \item{\code{max_ww_68}}{Maximum over pools of mass_tot_gww + mass_tot_gww_sd (g wet weight).}
#' }
#'
#' @details Length is in cm, wet mass in g. The middle columns are medians across pools, not averages across fish. Fields ending in 68 are extrema of source mean minus or plus a reported SD, and are not empirical quantiles or 68-percent confidence intervals.
#'
#' @source CHOPIN pipeline: \code{data-raw/02_CHOPIN_diet.R}; input biotic records in the \code{chopin.data} R package, plus literature-derived mass references where relevant.
#'
#' @docType data
#' @keywords datasets
"soles_contam_biometrics"
