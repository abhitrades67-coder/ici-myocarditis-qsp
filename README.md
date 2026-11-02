# ici-myocarditis-qsp

Code, data and outputs for the article:

> Anil A, Saravanan A, Roy SS, Srinivasan A. A quantitative systems pharmacology model of immune
> checkpoint inhibitor myocarditis does not reproduce clinical onset, incidence and troponin data:
> a history-matching calibration study.

The repository holds a population quantitative systems pharmacology (QSP) model of myocarditis
caused by immune checkpoint inhibitors, the virtual population engine and observation model, the
verification checks, and the history-matching calibration against 26 statistics from six clinical
data sources (C1 to C6). Every output reported in the article is included in `output/` and can be
regenerated with the commands below.

## Contents

```
R/          model, population engine, calibration and analysis code (see "Code" below)
src/        qsp2_rhs.c, the model's right-hand side in C (compiled on first use)
tests/      regression tests (testthat), run by blocks B4 to B8
data/       calibration targets, extracted source values, source records, and the C1 case file
docs/       model_specification.md, the full model specification (Additional file 1)
output/     every output of blocks B4 to B8, the post hoc summaries and the figures
```

### Code

| File | Role |
|---|---|
| `R/setup.R` | paths, packages, logging and shared helpers; compiles and loads `src/qsp2_rhs.c` |
| `R/run_block.R` | runs one block with its log and gate (see "Blocks and gates") |
| `R/qsp2_model.R` | parameters, dosing calendars and the ODE simulator |
| `R/qsp2_rhs_reference.R` | an independent R implementation of the right-hand side, used only by block B4 |
| `R/lib_pop.R` | the virtual population engine and the observation model |
| `R/lib_cal.R` | calibration cohorts, targets, priors, statistics and the implausibility rule |
| `R/lib_vpop.R`, `R/lib_mc.R` | a simpler population summary and a patient-by-patient Monte Carlo, used by the checks |
| `R/lib_lhs.R` | Latin hypercube sampling |
| `R/lib_parallel.R` | the parallel worker pool and the resumable cache |
| `R/run_tests.R` | runs the regression tests in `tests/` |
| `R/B04_implementation.R` | B4: implementation and verification |
| `R/B05_magnitudes.R` | B5: magnitude audit |
| `R/B06_feasibility.R` | B6: structural feasibility (requirements R1 to R5) |
| `R/B07_sensitivity.R` | B7: sensitivity (PRCC) and how much of each prior the feasible region keeps |
| `R/B08_calibration.R` | B8: history matching and the acceptance gate |
| `R/posthoc.R` | post hoc summaries of the B8 outputs (Additional file 5) |
| `R/figure1.R`, `R/figures.R` | Figure 1, and Figures 2 to 5 |
| `R/additional_files.R` | Additional files 2, 3 and 5 |
| `R/faers_c1_cases.R` | builds the C1 case file from the FAERS quarterly files |

### Data

| File | Content |
|---|---|
| `data/targets_calibration.csv` | definitions of the six calibration data sources C1 to C6 |
| `data/extraction.csv` | every value used for C1 to C6, the fixed physiology and the priors, with the verbatim sentence and its location in the source |
| `data/sources.csv` | bibliographic records of the sources (source PDFs are not redistributed) |
| `data/c1_cases_precision.csv` | the FAERS myocarditis cases behind C1, one row per case, with the precision of each date |

## Requirements

- R 4.5.2 and a C compiler (Rtools 4.5 on Windows; the standard toolchain on Linux or macOS).
- R packages, with the versions used: deSolve 1.41, data.table 1.18.2.1, matrixStats 1.5.0,
  jsonlite 2.0.0, digest 0.6.39, testthat 3.3.2 and ggplot2 4.0.3; for the figures and additional
  files also gridExtra 2.3 and openxlsx 4.2.8.1.

```r
install.packages(c("deSolve", "data.table", "matrixStats", "jsonlite", "digest", "testthat",
                   "ggplot2", "gridExtra", "openxlsx"))
```

Run every command from the repository folder (the folder that holds `.qsp2_root`). In RStudio,
`setwd()` to that folder first.

## Reproducing the results

### From the saved outputs (minutes)

The post hoc summaries, the figures and the additional files read only `output/`; no model is run.

```sh
Rscript R/posthoc.R            # output/posthoc_b08_summary.txt and output/posthoc_b08_*.csv
Rscript R/figure1.R            # output/figures/Figure_1.pdf and .png
Rscript R/figures.R            # output/figures/Figure_2 to Figure_5 (.pdf and .png)
Rscript R/additional_files.R   # output/additional_files/ (Additional files 2, 3 and 5)
```

### Running the model (blocks B4 to B8)

```sh
Rscript R/run_block.R B04   # implementation and verification, about 10 minutes
Rscript R/run_block.R B05   # magnitude audit, a few minutes
Rscript R/run_block.R B06   # structural feasibility, 2,000 parameter sets, about 22 hours on 10 workers
Rscript R/run_block.R B07   # sensitivity, about a minute
Rscript R/run_block.R B08   # calibration: 8 waves of 600 sets and an 81-node re-evaluation of
                            # 440 sets, about 70 hours on 10 workers; it exits with an error
                            # because the acceptance gate fails, as reported in the article
```

Each block writes `logs/<block>.log`, its outputs in `output/`, a full record of its checks in
`output/<block>_gate.json` and one line in `output/00_block_status.csv`.

- **Workers.** Blocks B6 and B8 run on `QSP2_WORKERS` parallel R processes (default 10, at most
  the number of cores less 2). `QSP2_WORKERS=1` runs serially and gives the same numbers.
- **Resumable runs.** Set `QSP2_CACHE_DIR` to a folder before B6 or B8. Every finished parameter
  set is then saved under a key that hashes the code (`R/lib_*.R`, `R/qsp2_model.R`, `R/setup.R`,
  `src/*.c`), the data files (`data/*.csv`), the R and deSolve versions, and the job itself. An
  interrupted run started again repeats only the sets that were in progress.
- **Seeds.** B6 uses the seed 20260928 and B8 the seed 20260929; the C1 bootstrap uses 20260928.
  With the same R and package versions the outputs are identical to those in `output/`.

### Blocks and gates

The blocks run in the order B4, B5, B6, B7, B8. A block runs only when the latest run of the block
before it passed (the record is in `output/00_block_status.csv`, included here, so any single block
can be rerun against the saved outputs). Each block registers its checks, and its gate passes only
if every check passes. Calibration blocks read only `data/targets_calibration.csv` through
`read_targets()`.

### The C1 case file from FAERS

`data/c1_cases_precision.csv` is built from the FDA Adverse Event Reporting System quarterly data
extract files, first quarter of 2015 to second quarter of 2026 (46 quarters), available at
https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html. Download the ASCII zip file of
each quarter (about 16 GB of text files once extracted), extract the DEMO, DRUG, REAC and THER
`.txt` files of all quarters into one folder, and run:

```sh
Rscript R/faers_c1_cases.R <folder with the FAERS .txt files> c1_rebuilt.csv
```

This takes about 10 minutes and gives a file identical to `data/c1_cases_precision.csv` (2,237 cases). The case
definition is in the header of `R/faers_c1_cases.R`; block B8 then removes reports from Japan,
keeps onsets within 365 days and keeps cases whose start and event dates are both complete to the
day.

## Outputs

| File | Content |
|---|---|
| `output/B04_verification.csv`, `B04_engine_vs_mc.csv`, `B04_implementation.json` | the verification checks, and the population engine against the patient Monte Carlo |
| `output/B05_magnitudes.csv`, `B05_nominal_pars.csv`, `B05_magnitudes.json` | the magnitude audit |
| `output/B06_lhs.csv`, `B06_feasible.csv`, `B06_waves.csv`, `B06_feasibility.json` | every set of the feasibility screen and the feasible region |
| `output/B07_sensitivity.csv`, `B07_profiles.csv`, `B07_sensitivity.json` | PRCC and surviving prior spans |
| `output/B08_targets.csv` | the 26 calibration targets with observational and model-discrepancy variances |
| `output/B08_all.csv`, `B08_waves.csv` | every one of the 4,800 parameter sets with its best population, and the waves |
| `output/B08_fine_grid.csv` | the 440 candidate sets re-evaluated on 81 S nodes |
| `output/B08_pairs.csv` | every kept population, with the flags for each alternative reading |
| `output/B08_ensemble.csv`, `B08_acceptance.csv`, `B08_intervals.csv`, `B08_c6_mix.csv`, `B08_calibration.json` | the accepted members, the acceptance test against each target, parameter ranges and the C6 regimen-mix check |
| `output/*_gate.json`, `output/00_block_status.csv` | every check of every block, and the gate results |
| `output/posthoc_b08_*` | post hoc summaries (Additional file 5) |
| `output/figures/` | Figures 1 to 5 of the article |

## Licence

The code is released under the MIT licence (`LICENSE`). The data and outputs are released under
the Creative Commons Attribution 4.0 International licence (CC BY 4.0). Quotations in
`data/extraction.csv` remain the copyright of their publishers and are reproduced for the
purpose of documenting the calibration targets.

## Contact

Anand Srinivasan, Department of Pharmacology, All India Institute of Medical Sciences (AIIMS)
Bhubaneswar, Odisha 751019, India; anandsrinivasan@aiimsbhubaneswar.edu.in
