# =============================================================================
# R/gbtm_k78.R: GBTM at k = 7-8 for the clusterBMA input set (decided 2026-11-15,
# docs/DECISIONS_memory.md). Six tasks: gbtm x raw/t/log x one/full family, age time,
# ge3_imp (the anchor configuration), quadratic (fit_gbtm, the main-run fitter).
#   Rscript R/gbtm_k78.R fit [cores]   # one job per task x k, into output/universes/gbtm_k78/
#   Rscript R/gbtm_k78.R merge         # append to the full_refit task files
# Seeds: the task's own seed (SEED + its row in spec_tasks), as 03_run.R gives it for
# every k, so these fits are what 03_run.R would produce with ks = 2:8 for these tasks
# (lcmm_fit sets the seed per k). k-means, kml and PAM already had k = 7-8.
# fit is restartable (finished task x k files are skipped; temp-and-rename saves).
# merge moves each original full_refit task file to output/archive/pre_k78_2026-11-18/
# before writing the version with k = 7-8 appended; then rerun 02_spec.R (GBTM k = 7-8
# supplement universes for these tasks) and 04_poverty.R full_refit.
# Note: R/gmm_refit_merge.R copies non-GMM task files from universes/full, so rerunning
# it would drop these fits; rerun this merge after it.
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
suppressPackageStartupMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
step  <- args[1]; stopifnot(step %in% c("fit", "merge"))
cores <- if (length(args) >= 2) as.integer(args[2]) else 3
K78 <- 7:8
dest <- file.path(UNIV, "gbtm_k78")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

tasks <- readRDS(file.path(RES, "spec_tasks.rds")) |> mutate(all_ids = list(task_id)) |>
  filter(method == "gbtm", time == "age", config == "ge3_imp")
stopifnot(nrow(tasks) == 6)
jobs <- expand_grid(row = seq_len(nrow(tasks)), k = K78) |>
  mutate(task_id = tasks$task_id[row], file = file.path(dest, sprintf("%s_k%d.rds", task_id, k)))

if (step == "fit") {
  wide <- readRDS(file.path(DER, "wide.rds"))
  long <- readRDS(file.path(DER, "long.rds"))
  todo <- jobs[!file.exists(jobs$file), ]
  # Full-family, k = 8 first (slowest).
  todo <- todo[order(!grepl("_full_", todo$task_id), -todo$k), ]
  cat(sprintf("%d of %d GBTM k = 7-8 fits to run on %d cores\n", nrow(todo), nrow(jobs), cores))
  done <- mclapply(seq_len(nrow(todo)), function(i) {
    j <- todo[i, ]
    task <- tasks[j$row, ]; task$ks <- list(j$k)
    t0 <- Sys.time()
    out <- fit_task(task, wide, long, pilot = FALSE)
    saveRDS(out, paste0(j$file, ".tmp")); file.rename(paste0(j$file, ".tmp"), j$file)
    msg <- sprintf("%s  %-32s k=%d %6.1f min  %s", format(Sys.time(), "%m-%d %H:%M"), j$task_id, j$k,
                   as.numeric(difftime(Sys.time(), t0, units = "mins")),
                   if (!is.na(out[[1]]$fit$error)) out[[1]]$fit$error else
                     sprintf("conv %s, min class %.3f", out[[1]]$fit$converged, out[[1]]$fit$min_prop))
    cat(msg, "\n", file = file.path(LOGS, "gbtm_k78.log"), append = TRUE)
    msg
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- vapply(done, inherits, logical(1), "try-error")
  cat(sprintf("Finished: %d ok, %d worker errors\n", sum(!errs), sum(errs)))
  if (any(errs)) print(done[errs])
}

if (step == "merge") {
  stopifnot(all(file.exists(jobs$file)))
  arch <- file.path(OUT, "archive", "pre_k78_2026-11-18", "universes_full_refit")
  dir.create(arch, recursive = TRUE, showWarnings = FALSE)
  for (tid in tasks$task_id) {
    f <- file.path(UNIV, "full_refit", paste0(tid, ".rds"))
    x <- readRDS(f)
    have <- vapply(x$fits, \(r) r$fit$k, numeric(1))
    if (all(K78 %in% have)) { cat("already merged:", tid, "\n"); next }
    stopifnot(!any(K78 %in% have))
    new <- unlist(lapply(jobs$file[jobs$task_id == tid], readRDS), recursive = FALSE)
    x$fits <- c(x$fits, new)
    x$task$ks <- list(sort(c(have, K78)))
    x$k78_added <- Sys.Date()
    stopifnot(file.rename(f, file.path(arch, basename(f))))   # archive by moving
    saveRDS(x, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    cat(sprintf("merged %s: k = %s\n", tid, paste(vapply(x$fits, \(r) r$fit$k, numeric(1)), collapse = ",")))
  }
}
