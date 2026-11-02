# ---------------------------------------------------------------------------
# lib_parallel.R  -  a small PSOCK pool for the population blocks (B6 onwards).
#
# A deliberately plain design.  Windows has no fork(), so workers are fresh R processes that
# build the model themselves from the library files, which hold definitions only.
#
# QSP2_WORKERS caps the pool (default 10, and never more than the machine's cores less 2).
# Set QSP2_WORKERS=1 to run everything serially, which gives the same numbers.
#
# Work is sent by function NAME, so what a worker runs is the library code it sourced, and
# nothing from the calling block's environment travels with it.  Anything the function
# needs besides its item goes in `ctx`, which is sent to each worker once.
#
# QSP2_CACHE_DIR, when set, makes a long run resumable.  Each worker writes every finished job
# to that folder as it completes, under a key built from a fingerprint of the code (every
# library file, the C source, the data files, R and deSolve versions), the function name, the
# context and the job itself.  A rerun after a crash or a reboot reads back every job whose
# key matches and computes the rest, so only the jobs that were in flight are lost.  A job that
# returns an error status is not kept, so it is tried again.  Any change to the code or the
# data changes every key, so nothing stale can be read.  Unset, nothing is read or written.
# ---------------------------------------------------------------------------

QSP2_WORKERS <- as.integer(Sys.getenv("QSP2_WORKERS", unset = "10"))
QSP2_WORKER_FILES <- c("setup.R", "qsp2_model.R", "lib_lhs.R", "lib_pop.R", "lib_cal.R")

#' `files` are the library files each worker sources.
qsp2_cluster <- function(n = QSP2_WORKERS, files = QSP2_WORKER_FILES) {
  n <- min(n, max(1L, parallel::detectCores() - 2L))
  if (n <= 1L) return(NULL)
  cl <- parallel::makePSOCKcluster(n)
  init <- function(dir, files) {
    Sys.setenv(QSP2_DIR = dir)
    setwd(dir)
    for (f in files) source(file.path(dir, "R", f))
    TRUE
  }
  environment(init) <- globalenv()
  ok <- tryCatch(parallel::clusterCall(cl, init, QSP2_DIR, files),
                 error = function(e) list(FALSE))
  if (!all(unlist(ok))) {
    parallel::stopCluster(cl)
    stop("workers failed to initialise")
  }
  cat(sprintf("  parallel pool: %d workers\n", n))
  cl
}

qsp2_cluster_stop <- function(cl) if (!is.null(cl)) parallel::stopCluster(cl)

#' The fingerprint of everything that decides what a job returns: the library files, the C
#' source, the data files, and the versions of R and deSolve.  Over-inclusive on purpose, so a
#' change anywhere invalidates the cache rather than risking a stale result.
qsp2_code_fingerprint <- function() {
  files <- sort(c(list.files(file.path(QSP2_DIR, "R"), pattern = "^(lib_.*|qsp2_model|setup)[.]R$",
                             full.names = TRUE),
                  list.files(file.path(QSP2_DIR, "src"), pattern = "[.]c$", full.names = TRUE),
                  list.files(file.path(QSP2_DIR, "data"), pattern = "[.]csv$", full.names = TRUE)))
  h <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
  digest::digest(list(basename(files), unname(h), R.version.string,
                      as.character(utils::packageVersion("deSolve"))), algo = "sha256")
}

#' Apply the library function named `fname` to every element of X, as fname(x, ctx = ctx).
#' Load-balanced across the pool when there is one, serial when there is not.  With
#' QSP2_CACHE_DIR set, finished jobs are kept and read back on a rerun (see the header).
qsp2_papply <- function(cl, X, fname, ctx = NULL) {
  stopifnot(is.character(fname), length(fname) == 1L)
  cache <- Sys.getenv("QSP2_CACHE_DIR", unset = "")
  paths <- rep(NA_character_, length(X))
  if (nzchar(cache)) {
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    base <- digest::digest(list(qsp2_code_fingerprint(), fname,
                                digest::digest(ctx, algo = "sha256")), algo = "sha256")
    paths <- vapply(X, function(x) file.path(cache, paste0(
      digest::digest(list(base, x), algo = "sha256"), ".rds")), "")
    n_have <- sum(file.exists(paths))
    if (n_have) cat(sprintf("  cache: %d of %d jobs already done, read back\n", n_have, length(X)))
  }
  items <- Map(function(x, path) list(x = x, path = path), X, paths)
  names(items) <- names(X)
  # a job is computed only if it is not already in the cache; a result is written to a
  # temporary name and renamed, so a reboot mid-write cannot leave a damaged file behind
  run_one <- function(it, f, ctx) {
    if (!is.na(it$path) && file.exists(it$path)) {
      r <- tryCatch(readRDS(it$path), error = function(e) NULL)
      if (!is.null(r)) return(r)
    }
    r <- f(it$x, ctx = ctx)
    failed <- is.list(r) && is.character(r$status) && length(r$status) == 1L &&
      grepl("^error", r$status)
    if (!is.na(it$path) && !failed) {
      tmp <- paste0(it$path, ".tmp", Sys.getpid())
      saveRDS(r, tmp)
      file.rename(tmp, it$path)
    }
    r
  }
  if (is.null(cl)) {
    f <- get(fname, mode = "function")
    return(lapply(items, run_one, f = f, ctx = ctx))
  }
  setter <- function(v) { assign(".qsp2_ctx", v, envir = globalenv()); TRUE }
  environment(setter) <- globalenv()
  parallel::clusterCall(cl, setter, ctx)
  environment(run_one) <- globalenv()
  g <- eval(bquote(function(it) .(run_one)(it, f = .(as.name(fname)),
                                           ctx = get(".qsp2_ctx", envir = globalenv()))),
            envir = globalenv())
  parallel::parLapplyLB(cl, items, g)
}
