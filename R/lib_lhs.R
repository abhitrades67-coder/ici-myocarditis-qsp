# ---------------------------------------------------------------------------
# lib_lhs.R  -  Latin hypercube sampling and history matching, in base R.
#
# Written here rather than taken from a package, so that every step can be read.  Both
# are short.
# ---------------------------------------------------------------------------

#' A Latin hypercube over a named list of ranges.
#' Each element is c(low, high, kind) with kind one of "unif", "log", "int".
qsp2_lhs <- function(n, prior) {
  stopifnot(n >= 1, length(prior) >= 1)
  out <- vector("list", length(prior))
  names(out) <- names(prior)
  for (j in seq_along(prior)) {
    spec <- prior[[j]]
    lo <- as.numeric(spec[1]); hi <- as.numeric(spec[2]); kind <- spec[3]
    # one stratified draw per interval, then shuffled, which is what makes it Latin
    u <- (sample.int(n) - 1 + runif(n)) / n
    out[[j]] <- switch(kind,
      unif = lo + u * (hi - lo),
      log  = exp(log(lo) + u * (log(hi) - log(lo))),
      int  = floor(lo + u * (hi - lo + 1 - 1e-12)),
      stop("unknown kind: ", kind))
  }
  as.data.table(out)
}

#' A Latin hypercube over ONE box: for every parameter, from the smallest to the largest value
#' among the points in `keep`, widened on each side by `inflate` times that span (on the log scale
#' for "log" parameters) and clipped to the prior.  No draw is filtered by its distance to the
#' points, so with points far apart most of the box lies away from all of them.  Block B8 uses
#' qsp2_lhs_near below instead.
qsp2_lhs_in_region <- function(n, prior, keep, inflate = 0.1) {
  nm <- names(prior)
  box <- prior
  for (j in nm) {
    v <- keep[[j]]
    spec <- prior[[j]]
    lo0 <- as.numeric(spec[1]); hi0 <- as.numeric(spec[2]); kind <- spec[3]
    if (kind == "log") {
      lv <- log(v); pad <- inflate * (max(lv) - min(lv) + 1e-9)
      box[[j]] <- c(max(lo0, exp(min(lv) - pad)), min(hi0, exp(max(lv) + pad)), kind)
    } else {
      pad <- inflate * (max(v) - min(v) + 1e-9)
      box[[j]] <- c(max(lo0, min(v) - pad), min(hi0, max(v) + pad), kind)
    }
  }
  qsp2_lhs(n, box)
}

#' Draws around EACH point of `keep` separately: every point gets its own small
#' box, `half` times the prior's range on either side of it (on the log scale for "log"
#' parameters), clipped to the prior, and an equal share of the n draws, each share a Latin
#' hypercube within its box.  Unlike one box spanning all the points, this keeps the draws near
#' the points that came closest, whichever corners of the prior they sit in.
qsp2_lhs_near <- function(n, prior, keep, half = 0.1) {
  k <- nrow(keep)
  stopifnot(k >= 1, n >= 1)
  per <- rep(n %/% k, k); per[seq_len(n %% k)] <- per[seq_len(n %% k)] + 1L
  out <- lapply(seq_len(k), function(i) {
    if (per[i] == 0) return(NULL)
    box <- prior
    for (j in names(prior)) {
      spec <- prior[[j]]; lo0 <- as.numeric(spec[1]); hi0 <- as.numeric(spec[2]); kind <- spec[3]
      v <- as.numeric(keep[[j]][i])
      if (kind == "log") {
        w <- half * (log(hi0) - log(lo0))
        box[[j]] <- c(max(lo0, exp(log(v) - w)), min(hi0, exp(log(v) + w)), kind)
      } else {
        w <- half * (hi0 - lo0)
        box[[j]] <- c(max(lo0, v - w), min(hi0, v + w), kind)
      }
    }
    qsp2_lhs(per[i], box)
  })
  d <- rbindlist(out)
  d[sample.int(nrow(d))]
}

#' Implausibility for one target.
#' A target may be a point value or a lower bound.
qsp2_implausibility <- function(model, obs, v_obs, v_md, v_sim = 0, bound = c("point", "lower")) {
  bound <- match.arg(bound)
  d <- model - obs
  if (bound == "lower") d <- pmin(d, 0)       # only a shortfall counts
  abs(d) / sqrt(v_obs + v_md + v_sim)
}

#' The rule: largest implausibility at most 3, second largest at most 2.5.
qsp2_not_implausible <- function(I, cut1 = 3, cut2 = 2.5) {
  if (is.null(dim(I))) I <- matrix(I, nrow = 1)
  apply(I, 1, function(r) {
    r <- sort(r[is.finite(r)], decreasing = TRUE)
    if (!length(r)) return(FALSE)
    r[1] <= cut1 && (length(r) < 2 || r[2] <= cut2)
  })
}
