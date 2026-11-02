# Figures 2 to 5 of the article, from saved outputs only (no model is run).
# Inputs: output/B08_all.csv, B08_fine_grid.csv, B08_waves.csv, B08_ensemble.csv,
#         B08_targets.csv, posthoc_b08_direction.csv, B07_profiles.csv
# Outputs: output/figures/Figure_2..5.pdf (vector) and .png (300 dpi)
# Run from the repository folder after R/posthoc.R:  Rscript R/figures.R
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
out <- "output"; fig <- file.path(out, "figures"); dir.create(fig, showWarnings = FALSE)
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREY <- "#9a9893"; INK <- "#0b0b0b"; INK2 <- "#52514e"
theme_set(theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.25),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
        plot.tag = element_text(face = "bold", size = 9), legend.position = "bottom",
        legend.title = element_blank(), legend.key.size = unit(3, "mm"), strip.text = element_text(colour = INK)))
save2 <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}
panel2 <- function(a, b) {  # two panels side by side without extra packages
  gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
  if (!gridExtra_ok) stop("install.packages('gridExtra') is needed to place two panels side by side")
  gridExtra::arrangeGrob(a, b, ncol = 2)
}

A  <- fread(file.path(out, "B08_all.csv"))
FG <- fread(file.path(out, "B08_fine_grid.csv"))
W  <- fread(file.path(out, "B08_waves.csv"))
E  <- fread(file.path(out, "B08_ensemble.csv"))
TG <- fread(file.path(out, "B08_targets.csv"))
DR <- fread(file.path(out, "posthoc_b08_direction.csv"))
PR <- fread(file.path(out, "B07_profiles.csv"))
req <- function(d) with(d, is.finite(R1_p7_mono) & R1_p7_mono >= 0.95 & is.finite(R3_t3_share_mono) &
                          R3_t3_share_mono < 0.5 & lost <= 1e-3 & edge_ok == 1)
A[, req_ok := req(A)]

# ---- Figure 2 ---------------------------------------------------------------------------
set.seed(1)
A[, xj := wave + runif(.N, -0.32, 0.32)]
ylim2 <- c(2.5, 30)
p2a <- ggplot(A, aes(xj, pmin(pmax(Imax, ylim2[1]), ylim2[2]))) +
  geom_point(colour = GREY, size = 0.35, alpha = 0.45, stroke = 0) +
  geom_hline(yintercept = 3, linetype = "dashed", colour = INK2, linewidth = 0.3) +
  geom_line(data = W, aes(wave, best_so_far), inherit.aes = FALSE, colour = BLUE, linewidth = 0.6) +
  geom_point(data = W, aes(wave, best_so_far), inherit.aes = FALSE, colour = BLUE, size = 1.4) +
  geom_text(data = W, aes(wave, ylim2[2] * 1.12, label = accepted_sets), inherit.aes = FALSE,
            size = 2.4, colour = INK) +
  scale_x_continuous(breaks = 1:8, limits = c(0.5, 8.5)) +
  scale_y_log10(breaks = c(2.5, 3, 4, 6, 10, 20, 30), limits = c(ylim2[1], ylim2[2] * 1.2)) +
  labs(x = "Wave", y = "Largest implausibility (25-node grid)", tag = "a",
       subtitle = "Numbers above each wave: sets accepted on the search grid")
FG[, grp := fifelse(accepted_81, "Accepted on 81 nodes", fifelse(accepted_25, "Accepted on 25 nodes only", "Not accepted"))]
FG[, grp := factor(grp, levels = c("Not accepted", "Accepted on 25 nodes only", "Accepted on 81 nodes"))]
lim_b <- range(c(FG$Imax_25, FG$Imax_81)); lim_b <- c(floor(lim_b[1] * 10) / 10, min(ceiling(lim_b[2] * 10) / 10, 5.5))
p2b <- ggplot(FG[order(grp)], aes(Imax_25, pmin(Imax_81, lim_b[2]), shape = grp, colour = grp, fill = grp)) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.3) +
  geom_hline(yintercept = 3, linetype = "dashed", colour = INK2, linewidth = 0.3) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = INK2, linewidth = 0.3) +
  geom_point(size = 1.3, stroke = 0.5) +
  scale_shape_manual(values = c(1, 4, 21)) +
  scale_colour_manual(values = c(GREY, INK, ORANGE)) +
  scale_fill_manual(values = c(NA, NA, ORANGE)) +
  coord_equal(xlim = lim_b, ylim = lim_b) +
  guides(shape = guide_legend(nrow = 2), colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(x = "Largest implausibility, 25-node grid", y = "Largest implausibility, 81-node grid", tag = "b")
save2(panel2(p2a, p2b), "Figure_2", 170, 90)

# ---- Figure 3 ---------------------------------------------------------------------------
lab <- c(C1_q25_combo = "C1 onset 25th pct, combination", C1_q50_combo = "C1 onset median, combination",
         C1_q75_combo = "C1 onset 75th pct, combination", C1_q25_mono = "C1 onset 25th pct, monotherapy",
         C1_q50_mono = "C1 onset median, monotherapy", C1_q75_mono = "C1 onset 75th pct, monotherapy",
         C1_ratio = "C1 ratio of medians", C2_p_combo = "C2 incidence, combination (lower bound)",
         C2_p_mono = "C2 incidence, monotherapy (lower bound)", C2_ratio = "C2 ratio of incidences",
         C3_p_damage = "C3 damage", C3_t_median = "C3 median time", C4_p_myo = "C4 myocarditis",
         C4_t_mean = "C4 mean time to diagnosis", C4_p_screen = "C4 found by troponin first (lower bound)",
         C4_ppv = "C4 myocarditis among elevations", C4_p_baseline = "C4 raised baseline",
         C5_p_pos = "C5 positive, any cause", C5_q25 = "C5 time to positive 25th pct",
         C5_q50 = "C5 time to positive median", C5_q75 = "C5 time to positive 75th pct",
         C5_ppv = "C5 myocarditis among positives", C5_t_myo = "C5 time to positive, myocarditis",
         C6_p_rise_nodx = "C6 rise without myocarditis", C6_trop_dx = "C6 troponin T at diagnosis",
         C6_m_bg = "C6 size of background rises")
near <- A[req_ok & Imax <= 4]
obs <- setNames(TG$obs, TG$id)
tl <- rbindlist(lapply(names(lab), function(t) data.table(target = t, r = near[[t]] / obs[[t]])))
ts <- tl[is.finite(r), .(lo = quantile(pmax(r, 1e-3), 0.05), q1 = quantile(pmax(r, 1e-3), 0.25),
                         md = median(pmax(r, 1e-3)), q3 = quantile(pmax(r, 1e-3), 0.75),
                         hi = quantile(pmax(r, 1e-3), 0.95)), by = target]
ord <- DR[order(median_I)]$target
ts[, f := factor(lab[target], levels = lab[ord])]
em <- rbindlist(lapply(names(lab), function(t) data.table(target = t, r = E[[t]] / obs[[t]])))
em[, f := factor(lab[target], levels = lab[ord])]
p3 <- ggplot(ts, aes(y = f)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.35) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = GREY, linewidth = 0.35) +
  geom_crossbar(aes(x = md, xmin = q1, xmax = q3), width = 0.55, fill = "#dce9f8", colour = BLUE,
                linewidth = 0.3, middle.linewidth = 0.6) +
  geom_point(data = em, aes(x = pmax(r, 1e-3), y = f), shape = 21, size = 1.4, fill = ORANGE,
             colour = "white", stroke = 0.3) +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10), labels = c("0.01", "0.03", "0.1", "0.3", "1", "3", "10")) +
  coord_cartesian(xlim = c(0.01, 12)) +
  labs(x = "Model value / observed value (log scale)", y = NULL)
save2(p3, "Figure_3", 170, 140)

# ---- Figure 4 ---------------------------------------------------------------------------
tsd <- function(t) sqrt(TG[id == t]$v_obs + TG[id == t]$v_md)
band <- function(t) obs[[t]] + c(-2, 2) * tsd(t)
yb <- band("C6_trop_dx")
A[, grp := factor(fifelse(req_ok, "Meets both structural requirements", "Other sets"),
                  levels = c("Other sets", "Meets both structural requirements"))]
mk4 <- function(xv, xl, tag) {
  xb <- band(xv)
  ggplot(A[order(grp)], aes(.data[[xv]], C6_trop_dx, colour = grp)) +
    annotate("rect", xmin = max(xb[1], 1), xmax = xb[2], ymin = max(yb[1], 1), ymax = yb[2],
             fill = "#f0efec", colour = NA) +
    geom_point(size = 0.45, alpha = 0.6, stroke = 0) +
    annotate("point", x = obs[[xv]], y = obs[["C6_trop_dx"]], shape = 4, size = 2.6, stroke = 0.9, colour = INK) +
    scale_colour_manual(values = c(GREY, BLUE)) +
    scale_y_log10(breaks = c(3, 10, 30, 100, 300, 1000, 3000)) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
    labs(x = xl, y = "C6 troponin T at diagnosis (ng/L, log scale)", tag = tag)
}
save2(panel2(mk4("C5_t_myo", "C5 time to first positive in myocarditis (days)", "a"),
             mk4("C4_t_mean", "C4 mean time to diagnosis (days)", "b")), "Figure_4", 170, 90)

# ---- Figure 5 ---------------------------------------------------------------------------
PR[, parameter := factor(parameter, levels = PR[order(-span_frac)]$parameter)]
PR[, grp := factor(fifelse(identified, "Narrowed to about half of the prior range", "Not narrowed"),
                   levels = c("Not narrowed", "Narrowed to about half of the prior range"))]
p5 <- ggplot(PR, aes(span_frac, parameter, fill = grp)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.3) +
  scale_fill_manual(values = c(GREY, BLUE)) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Span among the 75 feasible sets / prior range", y = NULL)
save2(p5, "Figure_5", 85, 130)
cat("figures written to", fig, "\n")
