# =============================================================================
# R/llpa_reseed.R: refit the Gaussian-mixture (llpa) tasks with fit_llpa's seed.
#   Rscript --vanilla R/llpa_reseed.R [mode] [cores]      # mode default: full_refit
# The main run fitted llpa unseeded, and at the full sample mclust initialises on a
# random 2,000-child subset, so those fits cannot be reproduced. fit_llpa seeds
# now (seed = SEED + task position, as in 03_run.R). Writes the llpa task files in
# output/universes/<mode>/, replacing the copies R/gmm_refit_merge.R made of the
# main run's. output/universes/full keeps the unseeded originals for comparison.
# Decided 2026-11-07. Afterwards run 04_poverty.R, 05_agreement.R, 06_figures.R.
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
suppressPackageStartupMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "full_refit"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
stopifnot(mode != "full", dir.exists(file.path(UNIV, mode)))

wide  <- readRDS(file.path(DER, "wide.rds"))
long  <- readRDS(file.path(DER, "long.rds"))
tasks <- readRDS(file.path(OUT, "spec_tasks.rds")) |> mutate(all_ids = list(task_id))
todo  <- filter(tasks, method == "llpa")
dest  <- file.path(UNIV, mode)

done <- mclapply(seq_len(nrow(todo)), function(i) {
  task <- todo[i, ]
  t0 <- Sys.time()
  out <- fit_task(task, wide, long, pilot = FALSE)
  f <- file.path(dest, paste0(task$task_id, ".rds"))
  saveRDS(list(task = select(task, -all_ids), fits = out, reseeded = Sys.Date()), paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  sprintf("%-32s %5.1f min", task$task_id, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}, mc.cores = cores, mc.preschedule = FALSE)

errs <- vapply(done, inherits, logical(1), "try-error")
cat(sprintf("llpa tasks refitted: %d ok, %d worker errors\n", sum(!errs), sum(errs)))
if (any(errs)) print(done[errs])
