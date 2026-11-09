# =============================================================================
# R/gmm_refit.R: refit the main-run GMM tasks at k = 2-4, converging every start.
#   Rscript R/gmm_refit.R <n_starts> <cores> [task_regex] [ks]
#   e.g. Rscript R/gmm_refit.R 50 8
#        Rscript R/gmm_refit.R 2 2 gmm_log_wave_one_ge3_native 2
# Decided 2026-11-03 after archive/code/gmm_optima.R showed gridsearch's 10-iteration screen
# does not find the best optimum (README "GMM optima experiment").
#
# Each start: 10 hlme iterations from random(minit), as lcmm::gridsearch does, then
# convergence from there. No screening. Data come from task_data() in R/fit_task.R,
# so each task uses the same children and scale as the main run. One .rds per start
# in output/gmm_refit/ (restartable, atomic saves; seeded by task, k and start, so
# n_starts can be raised later). R/gmm_refit_merge.R picks the best per task x k.
# =============================================================================

source("00_setup.R")
source("R/fit_task.R")
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
n_starts <- as.integer(args[1])
cores <- as.integer(args[2])
regex <- if (length(args) >= 3) args[3] else "."
ks <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else 2:4
dest <- file.path(OUT, "gmm_refit")
dir.create(dest, showWarnings = FALSE)

wide <- readRDS(file.path(DER, "wide.rds"))
long <- readRDS(file.path(DER, "long.rds"))
all_tasks <- readRDS(file.path(OUT, "spec_tasks.rds")) |> filter(method == "gmm")
tasks <- filter(all_tasks, grepl(regex, task_id))
f <- lcmm_formulas(gmm = TRUE)
logf <- file.path(LOGS, "gmm_refit.log")

for (i in seq_len(nrow(tasks))) {
  task <- tasks[i, ]
  jobs <- expand_grid(k = ks, s = seq_len(n_starts)) |>
    mutate(out = file.path(dest, sprintf("%s_k%d_s%03d.rds", task$task_id, k, s))) |>
    filter(!file.exists(out))
  cat(sprintf("%s  %s: %d starts to run\n", format(Sys.time(), "%m-%d %H:%M"), task$task_id, nrow(jobs)))
  if (nrow(jobs) == 0) next

  d <- task_data(task, wide, long, pilot = FALSE)
  dl <- d$long; dl$id <- match(dl$id, d$ids)
  b <- list(fixed = f$fixed, random = f$random, subject = "id", data = dl, verbose = FALSE)
  # Bound globally before forking: hlme looks up random()'s argument from its
  # namespace (see lcmm_record in R/fit_task.R).
  assign(".lcmm_minit", do.call(lcmm::hlme, c(b, list(ng = 1))), envir = globalenv())
  task_idx <- match(task$task_id, all_tasks$task_id)

  res <- mclapply(seq_len(nrow(jobs)), function(j) {
    k <- jobs$k[j]; s <- jobs$s[j]
    seed <- as.integer(SEED + 2e6 + 1e4 * task_idx + 100 * k + s)
    set.seed(seed)
    t0 <- Sys.time()
    rec <- list(task_id = task$task_id, k = k, start = s, seed = seed, error = NA_character_,
                ll_short = NA_real_, ll_full = NA_real_, conv = NA_real_, bic = NA_real_,
                aic = NA_real_, entropy = NA_real_, niter = NA_real_, cls = NULL)
    r <- tryCatch({
      short <- do.call(lcmm::hlme, c(b, list(mixture = f$mixture, ng = k, maxiter = 10,
                                            B = quote(random(.lcmm_minit)))))
      full <- do.call(lcmm::hlme, c(b, list(mixture = f$mixture, ng = k, B = short$best)))
      p <- full$pprob
      list(ll_short = short$loglik, ll_full = full$loglik, conv = full$conv, bic = full$BIC,
           aic = full$AIC, entropy = rel_entropy(as.matrix(p[, -(1:2)])), niter = full$niter,
           cls = as.integer(p$class[match(seq_along(d$ids), p$id)]))
    }, error = function(e) list(error = conditionMessage(e)))
    rec[names(r)] <- r
    rec$mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    saveRDS(rec, paste0(jobs$out[j], ".tmp")); file.rename(paste0(jobs$out[j], ".tmp"), jobs$out[j])
    sprintf("k%d s%03d %.1f min", k, s, rec$mins)
  }, mc.cores = cores, mc.preschedule = FALSE)

  errs <- vapply(res, inherits, logical(1), "try-error")
  cat(sprintf("%s  %-40s %d starts, %d worker errors\n", format(Sys.time(), "%m-%d %H:%M"),
              task$task_id, nrow(jobs), sum(errs)), file = logf, append = TRUE)
  if (any(errs)) print(res[errs])
}
