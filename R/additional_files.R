# Additional files 2, 3 and 5 of the article, built from the data and output files.
# Run from the repository folder after R/posthoc.R:  Rscript R/additional_files.R [output folder]
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
args <- commandArgs(trailingOnly = TRUE)
dest <- if (length(args)) args[1] else file.path("output", "additional_files")
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
source("R/lib_cal.R", local = (e <- new.env()))  # for QSP2_PRIOR and the population grid only

note <- function(x) data.table(note = x)
wb_write <- function(sheets, file) {
  wb <- createWorkbook()
  for (nm in names(sheets)) { addWorksheet(wb, nm); writeData(wb, nm, sheets[[nm]]); setColWidths(wb, nm, 1:ncol(sheets[[nm]]), "auto") }
  saveWorkbook(wb, file, overwrite = TRUE)
}

# ---- Additional file 2: targets, extraction, fixed values, priors, sources ----------------
tg  <- fread("output/B08_targets.csv")[use == "calibration"]
tg[, `:=`(sd_obs = sqrt(v_obs), sd_model_discrepancy = sqrt(v_md))]
pri <- rbindlist(lapply(names(e$QSP2_PRIOR), function(n) {
  v <- e$QSP2_PRIOR[[n]]; data.table(parameter = n, lower = as.numeric(v[1]), upper = as.numeric(v[2]),
                                     scale = c(unif = "uniform", log = "log-uniform", int = "integer uniform")[[v[3]]]) }))
grid <- data.table(parameter = c("mu", "sigma", "sd_tau"),
                   values = c("-12 to 4 in steps of 0.5 (33 values)", "0.4 to 4.4 in steps of 0.25 (17 values)",
                              paste(sort(unique(e$QSP2_TRI_GRID$sd_tau)), collapse = ", ")),
                   role = c("log-scale mean of precursor supply S", "log-scale spread of S", "log-scale spread of tempo tau"))
wb_write(list(
  README = note(c("Additional file 2. Calibration targets, extraction, fixed values and priors.",
                  "Targets: the 26 calibration statistics with observed value, kind (point or lower bound), observation SD and model-discrepancy SD (10% of the observed value).",
                  "Target definitions: cohort-level definitions of C1 to C6.",
                  "Extraction: every value used for C1 to C6, fixed physiology and priors, with the verbatim sentence and its location in the source.",
                  "Priors: the 36 fitted model parameters; Population grid: the three population parameters evaluated on a grid.",
                  "Sources: bibliographic records of the sources used.")),
  Targets = tg[, .(id, observed = obs, kind, sd_obs, sd_model_discrepancy)],
  `Target definitions` = fread("data/targets_calibration.csv"),
  Extraction = fread("data/extraction.csv"),
  Priors = pri, `Population grid` = grid,
  Sources = fread("data/sources.csv")
), file.path(dest, "Additional file 2 - calibration targets and priors.xlsx"))

# ---- Additional file 3: calibration results ------------------------------------------------
wb_write(list(
  README = note(c("Additional file 3. Results of the calibration (block B8).",
                  "Waves: per-wave counts and best largest implausibility.",
                  "All sets: each of the 4,800 parameter sets with its best population on the 25-node S grid (parameters, statistics, implausibilities, requirements).",
                  "Fine grid: the 440 candidate sets re-evaluated on 81 S nodes; acceptance on 25 and on 81 nodes.",
                  "Ensemble: the 4 accepted members (81-node grid). Acceptance: ensemble median against each target. Intervals: parameter ranges over the members.",
                  "Every kept population (B08_pairs.csv, 29 MB) is in the code and data archive.")),
  Waves = fread("output/B08_waves.csv"), `All sets` = fread("output/B08_all.csv"),
  `Fine grid` = fread("output/B08_fine_grid.csv"), Ensemble = fread("output/B08_ensemble.csv"),
  Acceptance = fread("output/B08_acceptance.csv"), Intervals = fread("output/B08_intervals.csv")
), file.path(dest, "Additional file 3 - calibration results.xlsx"))

# ---- Additional file 5: post hoc analyses --------------------------------------------------
wb_write(list(
  README = note(c("Additional file 5. Post hoc, exploratory analyses of the calibration outputs (R/posthoc.R).",
                  "No model was run; every value is computed from Additional file 3 and the archived B08_pairs.csv.",
                  readLines("output/posthoc_b08_summary.txt"))),
  `Worst target` = fread("output/posthoc_b08_worst_target.csv"),
  Direction = fread("output/posthoc_b08_direction.csv"),
  `Timing and troponin` = fread("output/posthoc_b08_tradeoff.csv"),
  `Leave one out` = fread("output/posthoc_b08_leave_one_out.csv")
), file.path(dest, "Additional file 5 - post hoc analyses.xlsx"))
cat("written to", dest, "\n")
