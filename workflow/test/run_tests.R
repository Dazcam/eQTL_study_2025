library(testthat)

# Resolve the directory this script lives in
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  getwd()  # fallback, e.g. when run interactively
}

script_dir <- get_script_dir()

source(file.path(script_dir, "..", "scripts", "functions.R"))
test_dir(script_dir, reporter = "summary")
