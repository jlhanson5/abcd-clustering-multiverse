# =============================================================================
# ABCD 5.1: clustering multiverse of CBCL total problems x baseline income
# 00_setup.R: paths, constants, shared helpers. Sourced by every other script.
#
# Reads tables directly from the 5.1 tabulated release (no local copies).
# Derived analysis files go to data_derived/ (inside Dropbox, approved storage).
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Override on another machine with the ABCD_REL / ABCD_DICT environment variables
# (e.g. in ~/.Renviron) instead of editing these paths.
REL  <- Sys.getenv("ABCD_REL", "~/Library/CloudStorage/Dropbox/_Research_In_Progress/abcd-data-release-5.1/core")
DICT <- Sys.getenv("ABCD_DICT", "~/Library/CloudStorage/Dropbox/_Research_In_Progress/ABCD_datadict51.xlsx")

DER  <- "data_derived"
OUT  <- "output"
UNIV <- file.path(OUT, "universes")
FIG  <- file.path(OUT, "figures")
LOGS <- file.path(OUT, "logs")
RES  <- file.path(OUT, "results")     # results_*, agreement_*, spec_*, gmm_refit_summary_*, cbma_poverty (since 2026-11-09)
SUMM <- file.path(OUT, "summaries")   # summary_*.txt
for (d in c(RES, SUMM)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# Figure folders, named by what they show. Other keys (pilot, smoke) keep their own name.
FIG_DIR <- c(full_refit = "1_spec_curve", full = "2_main_run", stability = "3_stability",
             cbma = "4_clusterbma", extra = "5_single_vs_consensus")
fig_dir <- function(key) {
  d <- file.path(FIG, if (key %in% names(FIG_DIR)) FIG_DIR[[key]] else key)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

SEED <- 20261023

EV_BL <- "baseline_year_1_arm_1"
EV_Y1 <- "1_year_follow_up_y_arm_1"
EV_Y2 <- "2_year_follow_up_y_arm_1"
EV_Y3 <- "3_year_follow_up_y_arm_1"
EVENTS <- c(BL = EV_BL, Y1 = EV_Y1, Y2 = EV_Y2, Y3 = EV_Y3)   # BL-Y3 window
WAVES  <- names(EVENTS)

# Baseline combined household income bins (demo_comb_income_v2) -> midpoint $.
# Edge bins: <$5k -> 2,500; $200k+ -> 250,000. 777/999 (refuse/don't know) excluded.
INCOME_MID <- c(`1` = 2500, `2` = 8500, `3` = 14000, `4` = 20500, `5` = 30000,
                `6` = 42500, `7` = 62500, `8` = 87500, `9` = 150000, `10` = 250000)

AGE_CENTER <- 9          # years; age time metric is (age - 9)
CLASS_FLOOR <- 0.05      # solutions with any class < 5% are flagged (kept)
N_STARTS <- 20           # random starts per universe (revisit later)
PILOT_N <- 2000          # approximate pilot size (whole families sampled)

# Poverty model (decided 2026-11-09): the race/ethnicity-adjusted model is primary;
# the age/sex/site model is the sensitivity analysis ("_unadj" in outputs). Race is
# ABCD's 5-level race_ethnicity (01b_race.R -> data_derived/race.rds), reference White.
# Children with no race code are left out of the adjusted model only.
POV_F <- list(adj   = class ~ income_halving + age_bl + sex + site + race,
              unadj = class ~ income_halving + age_bl + sex + site)
RACE_LEV <- c(`1` = "White", `2` = "Black", `3` = "Hispanic", `4` = "Asian", `5` = "Other")
add_race <- function(d) {
  r <- readRDS(file.path(DER, "race.rds"))
  d$race <- factor(r$race[match(d$src_subject_id, r$src_subject_id)], levels = RACE_LEV)
  d
}
# Highest-vs-lowest RRR is not estimable (set to NA, both models) when the highest class
# has fewer than MIN_TOP_N children: with 21 site and 4 race dummies the model separates
# (decided 2026-11-09; 2 GMM k = 6 universes per run, classes of 1-2 children).
MIN_TOP_N <- 10
# Also not estimable: a fit that separates (SE zero, non-finite or above MAX_SE). Added
# 2026-11-09 after 6 stability solutions (classes of 12-22 children) gave RRRs of 1e37+.
# No spec-curve or clusterBMA estimate is affected.
MAX_SE <- 5
# Rows used by a poverty formula: drops missing race for the adjusted model.
pov_rows <- function(dat, f) if ("race" %in% all.vars(f)) dat[!is.na(dat$race), , drop = FALSE] else dat

# Read one release table, all columns as character, then keep `cols`.
read_abcd <- function(path, cols = NULL) {
  x <- read_csv(file.path(REL, path),
                col_types = cols(.default = col_character()),
                show_col_types = FALSE, progress = FALSE)
  if (!is.null(cols)) x <- select(x, src_subject_id, any_of("eventname"), all_of(cols))
  x
}

num <- function(x, na = c(777, 999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% na] <- NA
  x
}

# Append a line to the exclusion / decision log.
log_line <- function(..., file = file.path(LOGS, "sample_flow.txt")) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = file, append = TRUE)
}
