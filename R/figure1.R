# Figure 1 of the article: model structure, calibration sources and the gated workflow.
# Drawn with grid only; 170 x 128 mm.  Run from the repository folder:  Rscript R/figure1.R
suppressPackageStartupMessages(library(grid))
fig <- file.path("output", "figures"); dir.create(fig, showWarnings = FALSE, recursive = TRUE)
W <- 170; H <- 128
INK <- "#0b0b0b"; INK2 <- "#52514e"; FILL <- "#f0efec"; BLUEF <- "#dce9f8"; BLUE <- "#2a78d6"
ORANGE <- "#eb6834"; GREY <- "#9a9893"

mm <- function(v) unit(v, "mm")
box <- function(x0, y0, w, h, title, body = NULL, fill = FILL, col = INK2, lty = 1, tcol = INK,
                lwd = 0.6, tsize = 7, bsize = 6.1) {
  grid.roundrect(mm(x0), mm(y0), mm(w), mm(h), r = unit(1.2, "mm"), just = c("left", "bottom"),
                 gp = gpar(fill = fill, col = col, lty = lty, lwd = lwd))
  nt <- length(strsplit(title, "\n")[[1]])
  ty <- if (is.null(body)) y0 + h / 2 else y0 + h - 1.6
  grid.text(title, mm(x0 + w / 2), mm(ty), just = c("centre", if (is.null(body)) "centre" else "top"),
            gp = gpar(fontsize = tsize, fontface = "bold", col = tcol, lineheight = 0.95))
  if (!is.null(body))
    grid.text(body, mm(x0 + w / 2), mm(y0 + h - 2.4 - nt * tsize * 0.38), just = c("centre", "top"),
              gp = gpar(fontsize = bsize, col = tcol, lineheight = 1.05))
}
arr <- function(x0, y0, x1, y1, col = INK2, lty = 1)
  grid.lines(mm(c(x0, x1)), mm(c(y0, y1)), gp = gpar(col = col, lwd = 0.7, lty = lty, fill = col),
             arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
tag <- function(label, x, y, text) {
  grid.text(label, mm(x), mm(y), just = c("left", "top"), gp = gpar(fontsize = 9, fontface = "bold", col = INK))
  grid.text(text, mm(x + 4), mm(y - 0.3), just = c("left", "top"), gp = gpar(fontsize = 7, col = INK2))
}

draw <- function() {
  grid.newpage()
  pushViewport(viewport(width = mm(W), height = mm(H)))
  # ---- (a) model -------------------------------------------------------------------------
  tag("a", 1, 127, "Model modules and patient variation")
  box(2, 101, 32, 22, "Pharmacokinetics and\nreceptor occupancy",
      "nivolumab, ipilimumab,\npembrolizumab;\nPD-1 and CTLA-4\nblockade")
  box(40, 101, 34, 22, "Priming, expansion\nand trafficking",
      "T cells primed only\nunder blockade;\neffector pool with a\nceiling; transit to heart")
  box(80, 101, 42, 22, "Myocardium",
      "effector and licensed resident\nT cells, macrophages, Tregs,\nexhaustion; killing, repair,\nscar; ejection fraction")
  box(128, 101, 40, 22, "Troponin",
      "fast and slow release;\nclearance at the measured\nhalf-lives of troponin I\nand troponin T")
  arr(34, 112, 40, 112); arr(74, 112, 80, 112); arr(122, 112, 128, 112)
  box(80, 76, 42, 19, "Tiers and presentation",
      "T1 injury, T2 myocarditis,\nT3 severe; presentation when\ntroponin or ejection fraction\ncrosses its threshold")
  box(128, 76, 40, 19, "Diagnosis and treatment",
      "at a positive test after T2,\nor at presentation; the ICI\nis stopped and killing\nhalts after a lag")
  arr(101, 101, 101, 95); arr(148, 101, 148, 95); arr(122, 85.5, 128, 85.5)
  box(2, 76, 72, 19, "Patient variation, integrated on grids",
      "precursor supply S sets the size of the response;\ntempo τ sets its pace; baseline troponin,\nsymptom threshold and baseline ejection fraction\nvary by quadrature", fill = BLUEF, col = BLUE, lty = 2)
  arr(57, 95, 57, 101, col = BLUE); arr(74, 85.5, 80, 85.5, col = BLUE)
  # ---- (b) calibration sources -------------------------------------------------------------
  tag("b", 1, 70, "Calibration sources, each simulated with its own assay, threshold, testing calendar and follow-up")
  cards <- list(c("C1  FAERS", "onset quartiles by\nregimen; no screening"),
                c("C2  Safety databases", "incidence by regimen;\nno screening"),
                c("C3  Monthly hs-cTnI", "6 months; damage\nand its timing"),
                c("C4  cTnI before\neach infusion", "> 34 ng/L; myocarditis,\ntiming, raised baselines"),
                c("C5  cTnI before\neach dose", ">= 55 ng/L; positives,\ntiming, myocarditis"),
                c("C6  hs-cTnT", "background rises;\ntroponin at diagnosis"))
  for (i in seq_along(cards)) {
    x0 <- 2 + (i - 1) * 27.8
    box(x0, 49, 26.4, 15.5, cards[[i]][1], gsub(">= ", "≥ ", cards[[i]][2]), tsize = 6.6, bsize = 6)
  }
  # ---- (c) workflow -------------------------------------------------------------------------
  tag("c", 1, 43, "Gated workflow")
  box(2, 14, 27, 22, "B4  Verification", "72 checks; engine\nagrees with an\nindividual-patient\nMonte Carlo")
  box(34, 14, 27, 22, "B6  Feasibility", "2,000 sets;\n75 meet the\nstructural\nrequirements")
  box(66, 14, 27, 22, "B7  Sensitivity", "partial rank\ncorrelation;\n3 of 39 parameters\nnarrowed")
  box(98, 14, 32, 22, "B8  History matching", "8 waves of 600 sets\non a 25-node grid;\n440 re-evaluated\non 81 nodes")
  box(135, 14, 33, 22, "Calibration gate",
      "needs ≥ 500 sets and an\nensemble median\nI ≤ 2 on all 26 targets;\nresult: 4 sets, 6 missed", col = ORANGE, lwd = 1.1)
  arr(29, 25, 34, 25); arr(61, 25, 66, 25); arr(93, 25, 98, 25); arr(130, 25, 135, 25)
  box(98, 3, 70, 7, "Validation (B9) and schedule comparison (B10): not run", NULL,
      fill = "white", col = GREY, lty = 2, tcol = GREY, tsize = 6.6)
  arr(151.5, 14, 151.5, 10, col = GREY, lty = 2)
  popViewport()
}
cairo_pdf(file.path(fig, "Figure_1.pdf"), width = W / 25.4, height = H / 25.4); draw(); invisible(dev.off())
png(file.path(fig, "Figure_1.png"), width = W, height = H, units = "mm", res = 300, type = "cairo"); draw(); invisible(dev.off())
cat("Figure 1 written to", fig, "\n")
