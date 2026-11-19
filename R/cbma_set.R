# R/cbma_set.R: which clusterBMA input set the downstream scripts read (12b, 13, 14, 15).
# CBMA_SET = rule (default; primary, output/cbma) or four (sensitivity, output/cbma_four).
# Design (decided 2026-11-19, docs/DECISIONS_memory.md): S_Dbw is undefined (NaN or Inf)
# for some k = 8 inputs on the raw and log scales in both sets, so raw and log have no S_Dbw
# and no averaged fit (as in the cognition project's precedent); t keeps all four
# weightings. 9 fits per set. Fits in CBMA_NOT_RUN are by design, not missing: any other
# absent fit stops the script.

CBMA_SET <- Sys.getenv("CBMA_SET", "rule")
stopifnot(CBMA_SET %in% c("rule", "four"))
CBMA_DIR <- file.path(OUT, if (CBMA_SET == "rule") "cbma" else "cbma_four")
CBMA_SFX <- if (CBMA_SET == "rule") "" else "_four"
CBMA_FITS <- c("raw_ch", "t_ch", "log_ch", "t_sdbw", "raw_equal", "t_equal", "log_equal",
               "all_equal", "t_avg")
CBMA_NOT_RUN <- c("raw_sdbw", "raw_avg", "log_sdbw", "log_avg")
CBMA_NOT_RUN_WHY <- "S_Dbw undefined for some k = 8 inputs, so no S_Dbw or averaged fit"
CBMA_SET_LAB <- if (CBMA_SET == "rule") "primary set: k-means and kml, k 2-8" else
  "sensitivity set: k-means, kml, PAM and GBTM, k 2-8"

cbma_files <- function() {
  f <- sort(list.files(CBMA_DIR, "\\.rds$", full.names = TRUE))
  got <- sub("\\.rds$", "", basename(f))
  miss <- setdiff(CBMA_FITS, got)
  if (length(miss)) stop(sprintf("clusterBMA fits missing from %s: %s", CBMA_DIR, paste(miss, collapse = ", ")))
  extra <- setdiff(got, CBMA_FITS)
  if (length(extra)) stop(sprintf("unexpected clusterBMA fits in %s: %s", CBMA_DIR, paste(extra, collapse = ", ")))
  f
}
cbma_not_run_note <- sprintf("Not run by design: %s (%s).",
                             paste(sub("_", " ", CBMA_NOT_RUN), collapse = ", "), CBMA_NOT_RUN_WHY)
