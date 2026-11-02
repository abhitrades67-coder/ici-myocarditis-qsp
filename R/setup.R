# ---------------------------------------------------------------------------
# setup.R  -  paths, packages, logging, gates and shared helpers.
#
# This file holds DEFINITIONS ONLY.  It is sourced by R/run_block.R and by parallel
# workers, and its only side effects are creating the output folders and loading
# the compiled library.  Blocks exchange results through files in output/, never
# through the R session, and each block runs in its own environment.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(deSolve)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

# --- locations -------------------------------------------------------------
# QSP2_DIR is the repository folder.  In RStudio: setwd() to it first.  The marker file
# guards against writing anywhere else.
QSP2_DIR <- normalizePath(Sys.getenv("QSP2_DIR", unset = getwd()), mustWork = TRUE)
if (!file.exists(file.path(QSP2_DIR, ".qsp2_root")))
  stop("not the repository folder (no .qsp2_root marker): ", QSP2_DIR,
       "\n  In RStudio use setwd() to the repository folder first.")
OUT_DIR  <- file.path(QSP2_DIR, "output")
LOG_DIR  <- file.path(QSP2_DIR, "logs")
SRC_DIR  <- file.path(QSP2_DIR, "src")
DATA_DIR <- file.path(QSP2_DIR, "data")
DOC_DIR  <- file.path(QSP2_DIR, "docs")
for (.d in c(OUT_DIR, LOG_DIR)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)
rm(.d)

# --- compiled right-hand side -----------------------------------------------
# Built from src/qsp2_rhs.c with R CMD SHLIB (Rtools on Windows).  Rebuilt when the source
# is newer than the library, so an edited model cannot run against a stale build.
QSP2_DLL <- file.path(SRC_DIR, paste0("qsp2_rhs", .Platform$dynlib.ext))
build_qsp2_dll <- function(force = FALSE) {
  src <- file.path(SRC_DIR, "qsp2_rhs.c")
  if (!file.exists(src)) stop("missing ", src)
  stale <- !file.exists(QSP2_DLL) || file.mtime(src) > file.mtime(QSP2_DLL)
  if (!(force || stale)) return(invisible(FALSE))
  if ("qsp2_rhs" %in% names(getLoadedDLLs())) dyn.unload(QSP2_DLL)
  owd <- setwd(SRC_DIR); on.exit(setwd(owd), add = TRUE)
  out <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", "qsp2_rhs.c"),
                 stdout = TRUE, stderr = TRUE)
  if (!file.exists(QSP2_DLL) || file.mtime(QSP2_DLL) < file.mtime(src))
    stop("compilation of qsp2_rhs.c failed:\n", paste(out, collapse = "\n"))
  invisible(TRUE)
}
load_qsp2_dll <- function() {
  build_qsp2_dll()
  if (!("qsp2_rhs" %in% names(getLoadedDLLs()))) dyn.load(QSP2_DLL)
  invisible(TRUE)
}
load_qsp2_dll()

# --- logging ----------------------------------------------------------------
# Each block writes a plain-text transcript to logs/<block>.log AND to the console.
start_log <- function(name) {
  path <- file.path(LOG_DIR, paste0(name, ".log"))
  con <- file(path, open = "wt")
  sink(con, split = TRUE)
  assign(".qsp_log_con", con, envir = globalenv())
  assign(".qsp_log_t0", Sys.time(), envir = globalenv())
  cat(strrep("=", 78), "\n", sep = "")
  cat(sprintf("%s   R %s   %s\n", name, getRversion(), format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(strrep("=", 78), "\n", sep = "")
  invisible(path)
}
end_log <- function() {
  t0 <- get(".qsp_log_t0", envir = globalenv())
  cat(sprintf("\n[block finished in %.1f s]\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  sink()
  close(get(".qsp_log_con", envir = globalenv()))
  invisible(NULL)
}
#' cat() that reaches the log file immediately (the log connection is buffered).
say <- function(...) {
  cat(...)
  con <- get0(".qsp_log_con", envir = globalenv())
  if (!is.null(con) && isOpen(con)) try(flush(con), silent = TRUE)
  flush.console()
}
#' The session details every log carries, so a result can be tied to its software.
log_session_info <- function() {
  si <- sessionInfo()
  cat("\n--- session ---\n")
  cat(sprintf("R %s (%s), platform %s\n", getRversion(), si$R.version$nickname, si$platform))
  pk <- c(si$otherPkgs, si$loadedOnly)
  v <- vapply(pk, function(p) p$Version, character(1))
  cat("packages:", paste(sprintf("%s %s", names(v), v)[order(names(v))], collapse = "; "), "\n")
  cat("---------------\n\n")
  invisible(si)
}

# --- small helpers ----------------------------------------------------------------
fnum <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")
save_json <- function(obj, filename) {
  path <- file.path(OUT_DIR, filename)
  write(toJSON(obj, auto_unbox = TRUE, digits = 12, null = "null", na = "null", pretty = TRUE), path)
  invisible(path)
}
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

# --- calibration data only --------------------------------------------------------
# Calibration code never sees validation data.  Two guards: read_targets() serves the
# calibration targets only, and run_block() refuses to start a calibration block whose
# source even names a validation targets file.
CALIBRATION_BLOCKS <- c("B04", "B05", "B06", "B07", "B08")
read_targets <- function(kind = "calibration") {
  if (!identical(kind, "calibration")) stop("only the calibration targets can be read")
  path <- file.path(DATA_DIR, "targets_calibration.csv")
  if (!file.exists(path)) stop("missing ", path)
  fread(path, encoding = "UTF-8")
}
#' TRUE when a block's source text names a validation targets file (static guard).
names_validation_file <- function(lines) any(grepl("targets_validation", lines, fixed = TRUE))
