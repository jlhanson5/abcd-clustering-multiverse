# One-off (2026-10-23): delete full-run task files whose fits all failed because
# one child had no interview_age ("time" contains NAs), so 03_run.R refits them.
f <- list.files("output/universes/full", pattern = "_age_.*\\.rds$", full.names = TRUE)
bad <- f[vapply(f, function(x) {
  e <- vapply(readRDS(x)$fits, function(r) as.character(r$fit$error), "")
  all(!is.na(e) & grepl("time column \"time\" contains NAs", e, fixed = TRUE))
}, logical(1))]
cat("Removing", length(bad), "failed age-task files\n"); print(basename(bad))
file.remove(bad)
