# Run the regression tests.  Returns a list with the counts and ok = TRUE if all pass.
qsp2_run_tests <- function() {
  library(testthat)
  env <- new.env(parent = globalenv())
  sys.source(file.path(QSP2_DIR, "R", "qsp2_model.R"), envir = env)
  sys.source(file.path(QSP2_DIR, "R", "qsp2_rhs_reference.R"), envir = env)
  # The observation model is as much a part of what has to keep working as the
  # differential equations are.
  sys.source(file.path(QSP2_DIR, "R", "lib_vpop.R"), envir = env)
  # the population engine, which every block from B6 on uses
  sys.source(file.path(QSP2_DIR, "R", "lib_pop.R"), envir = env)
  # the calibration definitions and the samplers (C1, the targets, R1, the search, end cells)
  sys.source(file.path(QSP2_DIR, "R", "lib_lhs.R"), envir = env)
  sys.source(file.path(QSP2_DIR, "R", "lib_cal.R"), envir = env)
  res <- test_file(file.path(QSP2_DIR, "tests", "test_qsp2_model.R"),
                   env = env, reporter = "silent")
  df <- as.data.frame(res)
  list(n = sum(df$nb), failed = sum(df$failed), errors = sum(df$error),
       ok = sum(df$failed) == 0 && sum(df$error) == 0)
}
