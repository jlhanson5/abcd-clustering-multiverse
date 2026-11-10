# =============================================================================
# 03_run.R: fit every task in the spec table, in parallel, one .rds per task.
#   Rscript 03_run.R pilot [cores] [method,method,...]
#   Rscript 03_run.R full  [cores] [method,method,...]
# Tasks whose .rds already exists are skipped, so the run can be restarted.
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
suppressPackageStartupMessages(library(parallel))

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "pilot"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
only  <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else NULL
stopifnot(mode %in% c("pilot", "full"))

wide  <- readRDS(file.path(DER, "wide.rds"))
long  <- readRDS(file.path(DER, "long.rds"))
tasks <- readRDS(file.path(RES, "spec_tasks.rds")) |> mutate(all_ids = list(task_id))
if (!is.null(only)) tasks <- filter(tasks, method %in% only)

dest <- file.path(UNIV, mode)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
todo <- tasks[!file.exists(file.path(dest, paste0(tasks$task_id, ".rds"))), ]
# Slow methods first so the long fits don't end up running alone at the end.
todo <- todo[order(match(todo$method, c("gmm", "gbtm", "llpa", "kml", "pam", "kmeans", "ward"))), ]
cat(sprintf("[%s] %d of %d tasks to fit on %d cores\n", mode, nrow(todo), nrow(tasks), cores))

done <- mclapply(seq_len(nrow(todo)), function(i) {
  task <- todo[i, ]
  t0 <- Sys.time()
  out <- fit_task(task, wide, long, pilot = mode == "pilot")
  # Write to a temp name, then rename, so a crash mid-save never leaves a
  # truncated .rds that a restart would treat as done.
  f <- file.path(dest, paste0(task$task_id, ".rds"))
  saveRDS(list(task = select(task, -all_ids), fits = out), paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  msg <- sprintf("%s  %-40s %6.1f min", format(Sys.time(), "%H:%M:%S"), task$task_id,
                 as.numeric(difftime(Sys.time(), t0, units = "mins")))
  cat(msg, "\n", file = file.path(LOGS, paste0("run_", mode, ".log")), append = TRUE)
  msg
}, mc.cores = cores, mc.preschedule = FALSE)

errs <- vapply(done, inherits, logical(1), "try-error")
cat(sprintf("Finished: %d ok, %d worker errors\n", sum(!errs), sum(errs)))
if (any(errs)) print(done[errs])
writeLines(capture.output(sessionInfo()), file.path(LOGS, paste0("sessionInfo_", mode, ".txt")))
