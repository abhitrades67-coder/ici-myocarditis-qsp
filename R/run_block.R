# ---------------------------------------------------------------------------
# run_block.R  -  run ONE block of the analysis, with its log and its gate.
#
#   Rscript:  from the repository folder,  Rscript R/run_block.R B04
#   RStudio:  setwd() to the repository folder; source("R/run_block.R"); run_block("B04")
#
# The blocks run in the order B04, B05, B06, B07, B08.  What the runner guarantees:
#  * the block runs in its OWN environment, so it cannot silently rebind a name the
#    setup defined.  Any rebinding of a protected name in the global environment fails
#    the gate;
#  * a calibration block whose source names a validation targets file is refused before
#    it runs;
#  * a block runs only when the latest run of the block before it passed (B04 has no
#    predecessor), so the blocks cannot be run out of order or past a failure;
#  * the block must register at least one check through check(name, ok, detail),
#    and passes only if every check passes and it ran without error;
#  * the result is appended to output/00_block_status.csv and written in full to
#    output/<block>_gate.json;  a failed gate stops with an error (non-zero exit
#    under Rscript), so a chain of blocks cannot run past a failure.
# ---------------------------------------------------------------------------

if (!exists("QSP2_DIR", envir = globalenv(), inherits = FALSE))
  source(file.path("R", "setup.R"), local = globalenv())

BLOCK_ORDER <- c("B04", "B05", "B06", "B07", "B08")

run_block <- function(id) {
  i_block <- match(id, BLOCK_ORDER)
  if (is.na(i_block)) stop("unknown block ", id, "; the blocks are ", paste(BLOCK_ORDER, collapse = ", "))
  rdir <- file.path(QSP2_DIR, "R")
  cand <- list.files(rdir, pattern = paste0("^", id, "_.*\\.R$"))
  if (length(cand) != 1L)
    stop(sprintf("expected exactly one R/%s_*.R block file, found %d: %s",
                 id, length(cand), paste(cand, collapse = ", ")))
  file <- file.path(rdir, cand)
  src_lines <- readLines(file, warn = FALSE)
  if (id %in% CALIBRATION_BLOCKS && names_validation_file(src_lines))
    stop(sprintf("%s is a calibration block and its source names validation targets; refused", id))

  # The chain: a block may run only when the latest run of the block before it passed.
  if (i_block > 1L) {
    prev <- BLOCK_ORDER[i_block - 1L]
    sfile0 <- file.path(OUT_DIR, "00_block_status.csv")
    st <- if (file.exists(sfile0)) fread(sfile0, colClasses = "character") else data.table()
    last_pass <- if (nrow(st) && all(c("block", "pass") %in% names(st)))
      tail(st$pass[st$block == prev], 1L) else character(0)
    if (!length(last_pass) || !identical(toupper(last_pass), "TRUE"))
      stop(sprintf("%s cannot run: the latest run of %s has not passed (see output/00_block_status.csv)",
                   id, prev))
  }

  assign("CURRENT_BLOCK", id, envir = globalenv())
  protected <- setdiff(ls(globalenv(), all.names = TRUE), c(".qsp_log_con", ".qsp_log_t0", ".Random.seed"))
  protected_val <- mget(protected, envir = globalenv())

  start_log(id)
  cat(sprintf("block file: R/%s   sha256 %s\n", cand, sha256_file(file)))
  cat(sprintf("runner:     R/run_block.R   sha256 %s\n", sha256_file(file.path(rdir, "run_block.R"))))
  log_session_info()

  benv <- new.env(parent = globalenv())
  benv$.checks <- data.table(check = character(), pass = logical(), detail = character())
  benv$check <- function(name, ok, detail = "") {
    ok <- isTRUE(ok)
    cat(sprintf("  [%s]  %s\n", if (ok) "PASS" else "FAIL", name))
    if (nzchar(detail))
      for (ln in strsplit(sub("\n$", "", detail), "\n")[[1]]) cat("           ", ln, "\n", sep = "")
    benv$.checks <- rbind(benv$.checks, data.table(check = name, pass = ok, detail = detail))
    invisible(ok)
  }

  ran_ok <- TRUE; err <- ""
  tryCatch(sys.source(file, envir = benv, keep.source = FALSE),
           error = function(e) { ran_ok <<- FALSE; err <<- conditionMessage(e)
             cat("\n!! BLOCK STOPPED WITH AN ERROR: ", err, "\n", sep = "") })

  rebound <- protected[!vapply(protected, function(nm)
    identical(get0(nm, envir = globalenv(), inherits = FALSE), protected_val[[nm]]), logical(1))]
  # Put back whatever the block rebound (or removed) BEFORE anything below uses it:
  # the status record must be written with the real paths, not the block's values.
  for (nm in rebound) assign(nm, protected_val[[nm]], envir = globalenv())
  checks <- benv$.checks
  pass <- ran_ok && !length(rebound) && nrow(checks) > 0L && all(checks$pass)

  cat("\n", strrep("-", 78), "\n", sep = "")
  cat(sprintf("GATE %s: %s   (%d checks, %d failed%s%s)\n", id, if (pass) "PASSED" else "FAILED",
              nrow(checks), sum(!checks$pass),
              if (!ran_ok) "; stopped with an error" else "",
              if (length(rebound)) paste0("; rebound protected names: ", paste(rebound, collapse = ", ")) else ""))
  if (nrow(checks) == 0L) cat("  the block registered no checks, so it cannot pass\n")
  cat(strrep("-", 78), "\n", sep = "")
  end_log()

  status <- data.table(block = id, file = cand, checks = nrow(checks), failed = sum(!checks$pass),
                       error = err, rebound = paste(rebound, collapse = ";"), pass = pass)
  sfile <- file.path(OUT_DIR, "00_block_status.csv")
  fwrite(status, sfile, append = file.exists(sfile))
  save_json(list(block = id, file = cand, file_sha256 = sha256_file(file), pass = pass,
                 ran_without_error = ran_ok, error = err, rebound = rebound,
                 checks = checks), paste0(id, "_gate.json"))
  if (!pass) stop(sprintf("GATE FAILED for %s - see logs/%s.log; nothing downstream may run", id, id))
  invisible(TRUE)
}

if (!interactive()) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (length(.args) == 1L) run_block(.args[1]) else if (length(.args) > 1L)
    stop("give one block id, e.g. Rscript R/run_block.R B04")
}
