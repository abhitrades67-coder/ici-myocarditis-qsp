# ---------------------------------------------------------------------------
# faers_c1_cases.R  -  the C1 case file, built from the FAERS quarterly ASCII files.
#
# Usage, from the repository folder:
#   Rscript R/faers_c1_cases.R <folder with the FAERS .txt files> [output file]
#
# Input: the DEMO, DRUG, REAC and THER files of the 46 quarters from the first quarter of
# 2015 to the second quarter of 2026, in the FDA's dollar-delimited ASCII format, extracted
# into one folder (the files are found by name, e.g. DEMO15Q1.txt, DRUG26Q2.txt).
# Output (default data/c1_cases_precision.csv): one row per myocarditis case with a
# checkpoint inhibitor as primary or secondary suspect drug, with the onset in days from the
# earliest start of a suspect checkpoint inhibitor to the event, the two dates and the
# precision each was reported to ("day", "month" or "year"), the regimen and the country of
# occurrence.  Onsets from 0 to 730 days are kept; block B8 applies the C1 rules (reports from
# Japan removed, onset within 365 days, both dates complete to the day).
#
# The case definition:
#   * a report carries a myocarditis preferred term (MYOCARDITIS_PT) and no infective
#     myocarditis term (EXCLUDE_PT);
#   * the latest version of each case (caseid) among the reports that carry a myocarditis or
#     infective myocarditis term and name a checkpoint inhibitor;
#   * the checkpoint inhibitor is identified in the active-ingredient or drug-name field, first
#     matching pattern in the order PD-1, PD-L1, CTLA-4; only primary (PS) and secondary (SS)
#     suspect drugs are used;
#   * the regimen is combination when the suspect drugs include anti-CTLA-4 and anti-PD-(L)1;
#   * a date given only as a month or a year is completed to its first day, and its precision
#     is recorded so that such dates can be left out.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("give the folder holding the FAERS .txt files")
FAERS_DIR <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
OUT_FILE <- if (length(args) >= 2) args[2] else file.path("data", "c1_cases_precision.csv")
fnum <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")

# --- drug identification ---------------------------------------------------
# ORDER MATTERS: identify_ici assigns the FIRST matching pattern, so this list must stay in the
# order PD-1, PD-L1, CTLA-4.
PAT_PD1 <- c(
  nivolumab     = "\\bnivolumab\\b|\\bopdivo\\b",
  pembrolizumab = "\\bpembrolizumab\\b|\\bkeytruda\\b",
  cemiplimab    = "\\bcemiplimab\\b|\\blibtayo\\b",
  dostarlimab   = "\\bdostarlimab\\b|\\bjemperli\\b",
  toripalimab   = "\\btoripalimab\\b",
  tislelizumab  = "\\btislelizumab\\b")
PAT_PDL1 <- c(
  atezolizumab  = "\\batezolizumab\\b|\\btecentriq\\b",
  durvalumab    = "\\bdurvalumab\\b|\\bimfinzi\\b",
  avelumab      = "\\bavelumab\\b|\\bbavencio\\b")
PAT_CTLA4 <- c(
  ipilimumab    = "\\bipilimumab\\b|\\byervoy\\b",
  tremelimumab  = "\\btremelimumab\\b|\\bimjudo\\b")
ICI_PATTERNS <- c(PAT_PD1, PAT_PDL1, PAT_CTLA4)
CLASS_OF <- setNames(
  c(rep("PD1", length(PAT_PD1)), rep("PDL1", length(PAT_PDL1)), rep("CTLA4", length(PAT_CTLA4))),
  names(ICI_PATTERNS))
ICI_ANY <- paste0("(", paste(ICI_PATTERNS, collapse = "|"), ")")

# --- outcome terms (MedDRA preferred terms as they appear in REAC.PT) -------
MYOCARDITIS_PT <- c("myocarditis", "immune-mediated myocarditis", "autoimmune myocarditis",
                    "myocarditis post infection", "eosinophilic myocarditis", "giant cell myocarditis",
                    "hypersensitivity myocarditis")
# a report carrying any of these infective myocarditis terms is excluded
EXCLUDE_PT <- c("myocarditis infectious", "viral myocarditis", "myocarditis septic",
                "myocarditis bacterial", "enterovirus myocarditis", "cytomegalovirus myocarditis",
                "myocarditis toxoplasmal", "myocarditis mycotic", "coxsackie myocarditis",
                "influenza myocarditis")

# --- reading -----------------------------------------------------------------
quarter_files <- function(table) {
  f <- list.files(FAERS_DIR, pattern = paste0("^", table, ".*\\.txt$"),
                  full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  sort(unique(normalizePath(f)))
}

#' Read one FAERS quarter, lower-casing the header and keeping `want`.
read_quarter <- function(f, want) {
  hdr <- tryCatch(names(fread(f, sep = "$", nrows = 0L, header = TRUE, showProgress = FALSE)),
                  error = function(e) character(0))
  if (!length(hdr)) return(NULL)
  sel <- hdr[tolower(trimws(hdr)) %in% want]
  if (!length(sel)) return(NULL)
  dt <- tryCatch(fread(f, sep = "$", header = TRUE, select = sel, colClasses = "character",
                       quote = "", fill = TRUE, encoding = "Latin-1", showProgress = FALSE,
                       na.strings = c("", "NA")),
                 error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  setnames(dt, tolower(trimws(names(dt))))
  for (m in setdiff(want, names(dt))) dt[, (m) := NA_character_]
  dt
}

#' Concatenate all quarters of one table, filtering each quarter before binding.
scan_table <- function(table, want, filt) {
  files <- quarter_files(table)
  if (!length(files)) stop("no ", table, "*.txt under ", FAERS_DIR)
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    dt <- read_quarter(files[i], want)
    if (is.null(dt)) { cat("  unreadable:", basename(files[i]), "\n"); next }
    out[[i]] <- filt(dt)
    rm(dt); gc(verbose = FALSE, full = FALSE)
  }
  cat(sprintf("  %s: %d quarterly files\n", table, length(files)))
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

as_i64 <- function(x) suppressWarnings(as.numeric(x))

#' TRUE where either name field names a checkpoint inhibitor (tested on distinct strings).
ici_row_hit <- function(dt) {
  hit_field <- function(v) {
    u <- unique(v[!is.na(v)])
    if (!length(u)) return(rep(FALSE, length(v)))
    keep <- u[grepl(ICI_ANY, tolower(u), perl = TRUE)]
    if (!length(keep)) return(rep(FALSE, length(v)))
    v %chin% keep
  }
  h <- hit_field(dt$prod_ai) | hit_field(dt$drugname)
  h[is.na(h)] <- FALSE
  h
}

#' Tag rows with a checkpoint inhibitor and its target class, first pattern wins.
identify_ici <- function(dt) {
  txt <- tolower(paste0(ifelse(is.na(dt$prod_ai), "", dt$prod_ai), " | ",
                        ifelse(is.na(dt$drugname), "", dt$drugname)))
  u <- unique(txt)
  assign_u <- rep(NA_character_, length(u))
  for (nm in names(ICI_PATTERNS)) {
    todo <- is.na(assign_u)
    if (!any(todo)) break
    hit <- grepl(ICI_PATTERNS[[nm]], u[todo], perl = TRUE)
    assign_u[which(todo)[hit]] <- nm
  }
  ici <- assign_u[match(txt, u)]
  dt[, `:=`(text = txt, ici = ici, ici_class = unname(CLASS_OF[ici]))]
  dt[!is.na(ici)]
}

#' FAERS dates are YYYYMMDD, YYYYMM or YYYY; partial dates are completed to the first day of
#' the period, and the precision of each is returned.
parse_faers_date <- function(s) {
  x <- gsub("[^0-9]", "", ifelse(is.na(s), "", trimws(s)))
  n <- nchar(x)
  out <- rep(as.Date(NA), length(x))
  prec <- rep("", length(x))
  i8 <- which(n == 8L)
  if (length(i8)) { out[i8] <- as.Date(x[i8], "%Y%m%d"); prec[i8] <- "day" }
  i6 <- which(n == 6L)
  if (length(i6)) { out[i6] <- as.Date(paste0(x[i6], "01"), "%Y%m%d"); prec[i6] <- "month" }
  i4 <- which(n == 4L)
  if (length(i4)) { out[i4] <- as.Date(paste0(x[i4], "0101"), "%Y%m%d"); prec[i4] <- "year" }
  prec[is.na(out)] <- ""
  list(date = out, prec = prec)
}

classify_regimen <- function(classes) {
  has_ctla4 <- "CTLA4" %in% classes
  has_pd <- any(c("PD1", "PDL1") %in% classes)
  if (has_ctla4 && has_pd) return("COMBO (PD-[L]1 + CTLA-4)")
  if (has_ctla4) return("CTLA-4 monotherapy")
  if ("PD1" %in% classes) return("PD-1 monotherapy")
  if ("PDL1" %in% classes) return("PD-L1 monotherapy")
  "unclassified"
}

# --- the pass ------------------------------------------------------------------
cat("FAERS folder:", FAERS_DIR, "\n")
reac <- scan_table("REAC", c("primaryid", "pt"),
                   function(d) d[tolower(trimws(pt)) %chin% c(MYOCARDITIS_PT, EXCLUDE_PT)])
reac[, `:=`(pt_l = tolower(trimws(pt)), primaryid = as_i64(primaryid))]
drug <- scan_table("DRUG", c("primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai"),
                   function(d) d[which(ici_row_hit(d))])
drug <- identify_ici(drug)
drug[, `:=`(primaryid = as_i64(primaryid), caseid = as_i64(caseid), role_cod = toupper(trimws(role_cod)))]
keep_ids <- intersect(unique(reac$primaryid), unique(drug$primaryid))
cat(sprintf("  %s reports name a checkpoint inhibitor and a myocarditis or infective myocarditis term\n",
            fnum(length(keep_ids))))
demo <- scan_table("DEMO", c("primaryid", "caseid", "caseversion", "event_dt", "occr_country"),
                   function(d) d[as_i64(primaryid) %in% keep_ids])
demo[, `:=`(primaryid = as_i64(primaryid), caseid = as_i64(caseid),
            caseversion_n = suppressWarnings(as.numeric(caseversion)))]
demo[is.na(caseversion_n), caseversion_n := 0]
ther <- scan_table("THER", c("primaryid", "dsg_drug_seq", "start_dt"),
                   function(d) d[as_i64(primaryid) %in% keep_ids])
ther[, primaryid := as_i64(primaryid)]
setnames(ther, "dsg_drug_seq", "drug_seq")

# the latest version of each case
setorder(demo, caseid, caseversion_n)
demo <- demo[!duplicated(caseid, fromLast = TRUE)]
myo_ids <- setdiff(unique(reac[pt_l %chin% MYOCARDITIS_PT, primaryid]),
                   unique(reac[pt_l %chin% EXCLUDE_PT, primaryid]))
drug <- drug[role_cod %chin% c("PS", "SS")]
merged <- merge(drug[primaryid %in% myo_ids], ther, by = c("primaryid", "drug_seq"),
                all.x = TRUE, allow.cartesian = TRUE)
ev <- parse_faers_date(demo$event_dt[match(merged$primaryid, demo$primaryid)])
st <- parse_faers_date(merged$start_dt)
merged[, `:=`(event_date = ev$date, start_date = st$date, event_prec = ev$prec, start_prec = st$prec)]
merged[, onset_days := as.numeric(event_date - start_date)]
merged <- merged[primaryid %in% demo$primaryid]
regimen <- drug[primaryid %in% myo_ids, .(regimen = classify_regimen(unique(na.omit(ici_class)))),
                by = primaryid]
# onset from the earliest start of a suspect checkpoint inhibitor
SDC <- c("ici", "onset_days", "start_date", "event_date", "start_prec", "event_prec")
dated <- merged[!is.na(start_date) & !is.na(event_date)]
setorder(dated, start_date, drug_seq)
D3 <- dated[, .SD[1L], by = primaryid, .SDcols = SDC][onset_days >= 0 & onset_days <= 730]
D3 <- merge(D3, regimen, by = "primaryid", all.x = TRUE)
D3[, occr_country := demo$occr_country[match(primaryid, demo$primaryid)]]
setorder(D3, primaryid)
fwrite(D3[, .(primaryid, ici, onset_days, start_date, event_date, start_prec, event_prec,
              regimen, occr_country)], OUT_FILE)
cat(sprintf("written: %s (%s cases)\n", OUT_FILE, fnum(nrow(D3))))
print(D3[, .N, by = regimen][order(-N)])
