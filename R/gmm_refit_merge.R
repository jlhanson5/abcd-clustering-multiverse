# =============================================================================
# R/gmm_refit_merge.R: build a universes folder with the refitted GMM solutions.
#   Rscript R/gmm_refit_merge.R <mode>       # mode: e.g. full_refit
# Then run 04_poverty.R, 05_agreement.R and 06_figures.R with that mode. The main
# run's output/universes/full and results_full.rds are not changed.
#
# For each GMM task x k = 2-4 with starts in output/gmm_refit/: the best converged
# solution (conv == 1) replaces the main-run fit. Runs within TOL log-likelihood
# units count as one optimum. If no start converged, the main-run fit is kept and
# flagged. Every other task file is copied as it is. gmm_flag marks GMM universes:
#   "k5-6 main run"   not refitted, not reliably estimated (README)
#   "best found once" only one start reached the best converged solution
#   "better run not converged"  a non-converged run beat the best converged solution
#   "no start converged"        main-run fit kept
#   "replicated"      at least two starts reached it, and no such run
# Output: output/universes/<mode>/, output/gmm_refit_summary_<mode>.rds
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")

mode <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(mode), mode != "full")
TOL <- 0.5
src <- file.path(UNIV, "full")
dest <- file.path(UNIV, mode)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

wide <- readRDS(file.path(DER, "wide.rds"))
long <- readRDS(file.path(DER, "long.rds"))
tasks <- readRDS(file.path(RES, "spec_tasks.rds"))
files <- list.files(file.path(OUT, "gmm_refit"), "\\.rds$", full.names = TRUE)
recs <- lapply(files, readRDS)
key <- tibble(task_id = map_chr(recs, "task_id"), k = map_int(recs, "k"), row = seq_along(recs))
summ <- list()

for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  f <- file.path(src, paste0(task$task_id, ".rds"))
  # Non-GMM tasks come from the main run, except llpa tasks already refitted with a
  # seed by R/llpa_reseed.R (kept, so rerunning this merge does not undo them).
  if (task$method != "gmm") {
    out <- file.path(dest, basename(f))
    if (!(task$method == "llpa" && file.exists(out) && !is.null(readRDS(out)$reseeded)))
      file.copy(f, out, overwrite = TRUE)
    next
  }
  old <- readRDS(f)
  d <- NULL
  old$fits <- lapply(old$fits, function(rec) {
    k <- rec$fit$k
    rows <- key$row[key$task_id == task$task_id & key$k == k]
    if (!length(rows)) {
      flag <- if (k >= 5) "k5-6 main run" else "not refitted"
      rec$fit <- mutate(rec$fit, gmm_flag = flag, refit = FALSE)
      return(rec)
    }
    st <- bind_rows(lapply(recs[rows], \(r) as_tibble(r[c("start", "seed", "ll_short", "ll_full", "conv",
                                                          "bic", "aic", "entropy", "mins", "error")])))
    okc <- which(st$conv %in% 1)
    row0 <- tibble(task_id = task$task_id, k = k, starts = nrow(st), errors = sum(!is.na(st$error)),
                   converged = length(okc))
    if (!length(okc)) {
      summ[[length(summ) + 1]] <<- mutate(row0, at_best = 0L, nonconv_above = NA_integer_, ll_best = NA_real_,
                                          flag = "no start converged", cpu_mins = sum(st$mins))
      rec$fit <- mutate(rec$fit, gmm_flag = "no start converged", refit = FALSE)
      return(rec)
    }
    best <- okc[which.max(st$ll_full[okc])]
    ll <- st$ll_full[best]
    at_best <- sum(st$ll_full[okc] > ll - TOL)
    above <- sum(!(st$conv %in% 1) & !is.na(st$ll_full) & st$ll_full > -1e8 & st$ll_full > ll + TOL)
    flag <- if (above > 0) "better run not converged" else if (at_best < 2) "best found once" else "replicated"
    if (is.null(d)) d <<- task_data(task, wide, long, pilot = FALSE)
    cls <- recs[[rows[match(st$start[best], st$start)]]]$cls
    props <- tabulate(cls, k) / length(cls)
    summ[[length(summ) + 1]] <<- mutate(row0, at_best = at_best, nonconv_above = above, ll_best = ll,
                                        flag = flag, cpu_mins = sum(st$mins))
    list(fit = tibble(task_id = task$task_id, k = k, n = length(d$ids), seed = st$seed[best],
                      secs = sum(st$mins) * 60, error = NA_character_, converged = TRUE,
                      bic = st$bic[best], aic = st$aic[best], entropy = st$entropy[best],
                      asw = asw(d$X, cls, st$seed[best]), k_found = sum(props > 0),
                      min_prop = min(props), small_class = min(props) < CLASS_FLOOR,
                      gmm_flag = flag, refit = TRUE),
         assign = tibble(src_subject_id = d$ids, class = cls))
  })
  saveRDS(old, file.path(dest, basename(f)))
}

summ <- bind_rows(summ)
saveRDS(summ, file.path(RES, paste0("gmm_refit_summary_", mode, ".rds")))
cat(sprintf("Task files written: %d. GMM universes refitted: %d.\n", length(list.files(dest)), nrow(summ)))
print(dplyr::count(summ, k, flag), n = 30)
