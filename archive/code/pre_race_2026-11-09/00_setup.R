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
