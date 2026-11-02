# Post hoc, exploratory summaries of block B8's saved outputs (Additional file 5 of the article).
# Reads only output/B08_all.csv (each set's best population on the search grid),
# output/B08_pairs.csv (every kept population, with the 81-node results for candidates) and
# output/B08_targets.csv.  Runs no model.  Writes output/posthoc_b08_summary.txt and
# output/posthoc_b08_*.csv.  Run from the repository folder:  Rscript R/posthoc.R
suppressPackageStartupMessages(library(data.table))
out <- "output"
sink(file.path(out, "posthoc_b08_summary.txt"), split = TRUE)
cat("Post hoc summaries of B08 outputs (exploratory; not part of any gate)\n\n")
A <- fread(file.path(out, "B08_all.csv"))
P <- fread(file.path(out, "B08_pairs.csv"))
tg <- fread(file.path(out, "B08_targets.csv"))
cal <- tg[use == "calibration"]$id
cal <- setdiff(cal, "C5_p_presented")
ic <- paste0("I_", cal)
stopifnot(all(ic %in% names(A)), length(cal) == 26)
cohort <- sub("_.*", "", cal)

req <- function(d) with(d, is.finite(R1_p7_mono) & R1_p7_mono >= 0.95 & is.finite(R3_t3_share_mono) &
                          R3_t3_share_mono < 0.5 & lost <= 1e-3 & edge_ok == 1)
Im <- as.matrix(A[, ..ic]); Im[is.na(Im)] <- Inf
A[, worst := sub("^I_", "", ic[apply(Im, 1, which.max)])]
A[, req_ok := req(A)]

cat("1. Sets and their best population on the 25-node search grid\n")
cat(sprintf("   sets %d; meeting R1, R3, no loss and the end-cell rule: %d (%.1f%%)\n",
            nrow(A), sum(A$req_ok), 100 * mean(A$req_ok)))
for (cut in c(3, 3.5, 4)) cat(sprintf("   largest I <= %.1f: %d sets; also meeting the requirements: %d\n",
                                      cut, sum(A$Imax <= cut), sum(A$Imax <= cut & A$req_ok)))
cat("\n   target with the largest I, all sets (top 10):\n")
w_all <- A[, .N, by = worst][order(-N)][1:10]; print(w_all)
cat("\n   target with the largest I, sets meeting the requirements with largest I <= 4:\n")
w_near <- A[req_ok & Imax <= 4, .N, by = worst][order(-N)]; print(w_near)
fwrite(A[, .N, by = .(worst, near = req_ok & Imax <= 4)], file.path(out, "posthoc_b08_worst_target.csv"))

cat("\n2. Direction of the misses, sets meeting the requirements with largest I <= 4\n")
obs <- setNames(tg$obs, tg$id)
near <- A[req_ok & Imax <= 4]
dirt <- rbindlist(lapply(cal, function(t) {
  v <- near[[t]]; v <- v[is.finite(v)]
  data.table(target = t, observed = obs[[t]], n = length(v), model_median = median(v),
             q25 = quantile(v, 0.25), q75 = quantile(v, 0.75), share_above_obs = mean(v > obs[[t]]),
             median_I = median(near[[paste0("I_", t)]], na.rm = TRUE))
}))
print(dirt[order(-median_I)], digits = 3)
fwrite(dirt, file.path(out, "posthoc_b08_direction.csv"))
cat(sprintf("\n   C6_trop_dx over all 4,800 best populations: maximum %.1f ng/L against %g observed; above 182 (half): %d sets\n",
            max(A$C6_trop_dx, na.rm = TRUE), obs[["C6_trop_dx"]], sum(A$C6_trop_dx > 182, na.rm = TRUE)))

cat("\n3. Trade-off between timing and troponin, best populations of all sets\n")
tim <- c("C4_t_mean", "C5_t_myo", "C5_q50", "C3_t_median")
for (t in tim) cat(sprintf("   Spearman correlation of model C6_trop_dx with model %-11s %6.2f\n", t,
                           cor(A$C6_trop_dx, A[[t]], method = "spearman", use = "complete.obs")))
ok_t <- A$I_C4_t_mean <= 2 & A$I_C5_t_myo <= 2
ok_c <- A$I_C6_trop_dx <= 2
cat(sprintf("   I <= 2 on both C4_t_mean and C5_t_myo: %d sets; I <= 2 on C6_trop_dx: %d sets; on all three: %d sets\n",
            sum(ok_t, na.rm = TRUE), sum(ok_c, na.rm = TRUE), sum(ok_t & ok_c, na.rm = TRUE)))
cat(sprintf("   smallest I on C6_trop_dx among sets with both timing I <= 2: %.2f\n",
            if (any(ok_t, na.rm = TRUE)) min(A$I_C6_trop_dx[which(ok_t)]) else NA))
cat(sprintf("   smallest largest-of-two timing I among sets with C6_trop_dx I <= 2: %.2f\n",
            if (any(ok_c, na.rm = TRUE)) min(pmax(A$I_C4_t_mean, A$I_C5_t_myo)[which(ok_c)]) else NA))
fwrite(A[, .(set, wave, Imax, I2nd, req_ok, C6_trop_dx, C4_t_mean, C5_t_myo, C5_q50, C3_t_median,
             I_C6_trop_dx, I_C4_t_mean, I_C5_t_myo)], file.path(out, "posthoc_b08_tradeoff.csv"))

cat("\n4. Leaving out one cohort at a time (B8's own per-population flags).  Each set contributes its\n")
cat("   population with the smallest largest I over the remaining targets; the ensemble-median rule\n")
cat("   is then applied.\n")
loo <- rbindlist(lapply(c("C1", "C2", "C3", "C4", "C5", "C6"), function(k) {
  f <- paste0("pass_no_", k)
  keep <- P[get(f) == 1]
  rest <- ic[cohort != k]
  if (!nrow(keep)) return(data.table(left_out = k, sets = 0L, targets_median_I_over_2 = NA_character_,
                                     left_out_median_I = NA_character_))
  M <- as.matrix(keep[, ..rest]); M[is.na(M)] <- Inf
  keep[, imax_rest := apply(M, 1, max)]
  one <- keep[keep[, .I[which.min(imax_rest)], by = set]$V1]
  med <- sapply(rest, function(x) median(one[[x]]))
  medk <- sapply(ic[cohort == k], function(x) median(one[[x]]))
  data.table(left_out = k, sets = nrow(one),
             targets_median_I_over_2 = if (any(med > 2)) paste(sprintf("%s %.2f", sub("^I_", "", rest[med > 2]), med[med > 2]), collapse = "; ") else "none",
             left_out_median_I = paste(sprintf("%s %.2f", sub("^I_", "", names(medk)), medk), collapse = "; "))
}))
print(loo)
fwrite(loo, file.path(out, "posthoc_b08_leave_one_out.csv"))
sink()
