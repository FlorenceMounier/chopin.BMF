
# Run from the chopin.BMF project root

dir.create("results/provenance",
           recursive = TRUE,
           showWarnings = FALSE)

# R version and installed package versions
writeLines(
  capture.output(sessionInfo()),
  "results/provenance/sessionInfo.txt"
)

# Exact source file fingerprints
files_to_hash <- c(
  "inst/extdata/Gogina_data-2021-0002-01.xlsx",
  "DESCRIPTION",
  "data-raw/00_RUN_ALL_FROM_SCRATCH.R",
  "manuscript/Article.Rmd",
  "manuscript/Article_SI.Rmd"
)

stopifnot(all(file.exists(files_to_hash)))

hashes <- tibble::tibble(
  file = files_to_hash,
  sha256 = vapply(
    files_to_hash,
    function(path) {
      digest::digest(file = path, algo = "sha256")
    },
    character(1)
  )
)

readr::write_csv(
  hashes,
  "results/provenance/source_hashes.csv"
)

# Installed source-data package metadata
source_package <- packageDescription("chopin.data")

source_package_info <- data.frame(
  field = c("Version", "RemoteSha", "RemoteRef"),
  value = vapply(
    c("Version", "RemoteSha", "RemoteRef"),
    function(x) {
      value <- source_package[[x]]
      if (is.null(value)) NA_character_ else value
    },
    character(1)
  )
)

readr::write_csv(
  source_package_info,
  "results/provenance/chopin_data_version.csv"
)


# ============================================================
# Git provenance: record current commit and working tree status
# ============================================================

git_output <- function(args) {
  tryCatch(
    suppressWarnings(
      system2("git", args = args, stdout = TRUE, stderr = TRUE)
    ),
    error = function(e) character(0)
  )
}

# Check that we are inside a Git repository
git_root <- git_output(c("rev-parse", "--show-toplevel"))

if (length(git_root) == 1L &&
    is.null(attr(git_root, "status"))) {

  # Current commit SHA
  git_commit <- git_output(c("rev-parse", "HEAD"))

  # Current branch
  git_branch <- git_output(c("branch", "--show-current"))

  # Check whether tracked or untracked files have changed
  git_status <- git_output(
    c("status", "--porcelain", "--untracked-files=normal")
  )

  git_dirty <- length(git_status) > 0L

  git_info <- tibble::tibble(
    git_commit = git_commit[1],
    git_branch = if (length(git_branch)) git_branch[1] else NA_character_,
    git_dirty = git_dirty,
    generated_at = format(
      Sys.time(), "%Y-%m-%d %H:%M:%S %z"
    )
  )

  readr::write_csv(
    git_info,
    "results/provenance/git_provenance.csv"
  )

  message("Git commit: ", git_commit[1])

  if (git_dirty) {
    warning(
      "Git working tree contains uncommitted or untracked files. ",
      "The commit SHA alone does not fully identify the analysis code."
    )
  }

} else {
  warning(
    "No Git repository detected. Git provenance was not recorded."
  )
}
