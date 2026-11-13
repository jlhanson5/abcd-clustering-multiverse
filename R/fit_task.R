# =============================================================================
# R/fit_task.R: fit one task (method x data config) for every k in the task.
# fit_task() returns one record per k: hard (modal) class assignments plus fit
# indices. Sourced by 03_run.R.
# =============================================================================

suppressPackageStartupMessages({
  library(latrend)
  library(cluster)
  library(mclust)
})

# ---- Data for one task ---------------------------------------------------------
# X     : n x 4 wide matrix on the task's scale (imputed where needed). Used by the
#         vector methods, and by every method for the silhouette.
# long  : long data for kml / GBTM / GMM (id, time, y).
task_data <- function(task, wide, long_obs, pilot) {
  w <- wide
  if (pilot) w <- filter(w, pilot)
  if (task$family == "one") w <- filter(w, one_per_family)
  if (task$config == "all4") w <- filter(w, all4)

  raw <- as.matrix(w[paste0("imp_r_", WAVES)])
  tsc <- as.matrix(w[paste0("imp_t_", WAVES)])
  X <- switch(task$scale, raw = raw, t = tsc, log = log1p(raw))
  dimnames(X) <- list(w$src_subject_id, WAVES)

  tvar <- function(wave_num, age) if (task$time == "age") age - AGE_CENTER else wave_num
  if (task$config == "ge3_native") {
    lg <- long_obs |>
      filter(src_subject_id %in% w$src_subject_id, !is.na(tot_r)) |>
      transmute(id = src_subject_id, time = tvar(wave_num, age),
                y = switch(task$scale, raw = tot_r, t = tot_t, log = log1p(tot_r)))
  } else {
    ages <- as.matrix(w[paste0("age_", WAVES)])
    lg <- tibble(id = rep(w$src_subject_id, times = 4),
                 wave_num = rep(0:3, each = nrow(w)),
                 age = as.vector(ages), y = as.vector(X)) |>
      transmute(id, time = tvar(wave_num, age), y)
  }
  list(X = X, long = as.data.frame(arrange(lg, id, time)), ids = w$src_subject_id)
}

# ---- Fit indices ----------------------------------------------------------------
rel_entropy <- function(z) {               # 1 = perfectly separated posteriors
  k <- ncol(z); if (k < 2) return(NA_real_)
  z <- pmax(z, 1e-12)
  1 - sum(-z * log(z)) / (nrow(z) * log(k))
}

asw <- function(X, cls, seed, max_n = 3000) {
  if (length(unique(cls)) < 2) return(NA_real_)
  i <- seq_len(nrow(X))
  if (nrow(X) > max_n) { set.seed(seed); i <- sort(sample(i, max_n)) }
  # The subsample can hold a single class when the others are tiny (e.g. average
  # linkage splitting off singletons); silhouette() then returns NA, not a matrix.
  if (length(unique(cls[i])) < 2) return(NA_real_)
  mean(silhouette(cls[i], dist(X[i, , drop = FALSE]))[, "sil_width"])
}

# ---- Methods ----------------------------------------------------------------------
# Each returns list(cls = integer vector aligned to rownames(X)/ids, bic, aic,
# entropy, converged, model). `model` is what R/predict_class.R needs to assign
# new children; 03_run.R ignores it. The settings arguments (metric, model_name,
# shape, ri_only) are for the stability analyses; their defaults are the main
# multiverse settings.

fit_kmeans <- function(d, k, seed) {
  set.seed(seed)
  f <- kmeans(d$X, k, nstart = N_STARTS, iter.max = 100)
  list(cls = f$cluster, bic = NA, aic = NA, entropy = NA, converged = f$ifault == 0,
       model = f$centers)
}

# PAM via CLARA: N_STARTS subsamples of up to 2,000 children, each PAM'd, best
# medoids kept on the full data. With n <= 2,000 this is exact PAM.
fit_pam <- function(d, k, seed, metric = "euclidean") {
  set.seed(seed)
  f <- clara(d$X, k, metric = metric, samples = N_STARTS, sampsize = min(nrow(d$X), 2000),
             pamLike = TRUE, rngR = TRUE)
  list(cls = f$clustering, bic = NA, aic = NA, entropy = NA, converged = TRUE,
       model = f$medoids)
}

fit_ward <- function(d, k, seed, tree) {
  list(cls = cutree(tree, k), bic = NA, aic = NA, entropy = NA, converged = TRUE,
       model = NULL)
}

# Gaussian mixture on the 4 wave scores, full unconstrained covariance (VVV).
# mclust initialises from hierarchical clustering. That is deterministic only
# for n <= 2,000: above that (mclust.options("subset")) it runs on a random
# subset of 2,000, so the fit depends on the RNG. Found 2026-10-24; the seed was
# added then, so the main run's full-sample fits (unseeded) are not exactly
# reproducible (refits with this seed agree at ARI .28-.84).
fit_llpa <- function(d, k, seed, model_name = "VVV") {
  set.seed(seed)
  f <- Mclust(d$X, G = k, modelNames = model_name, verbose = FALSE)
  if (is.null(f)) stop(sprintf("Mclust returned NULL (no %s solution)", model_name))
  list(cls = f$classification, bic = -f$bic, aic = -2 * f$loglik + 2 * f$df,
       entropy = rel_entropy(f$z), converged = TRUE, model = f)
}

latrend_assign <- function(m, ids) {
  a <- trajectoryAssignments(m)
  names(a) <- ids(m)
  as.integer(a[ids])
}

fit_kml <- function(d, k, seed) {
  set.seed(seed)
  m <- latrend(lcMethodKML(response = "y", id = "id", time = "time", nClusters = k,
                           nbRedrawing = N_STARTS, seed = seed),
               data = d$long, verbose = FALSE)
  list(cls = latrend_assign(m, d$ids), bic = NA, aic = NA, entropy = NA, converged = TRUE,
       model = NULL)
}

# Class-specific quadratic in time. GMM adds a random intercept and slope; GBTM
# has no random effects. latrend's lcMethodLcmmGBTM ignores init = "gridsearch"
# (one start only), so GBTM is fitted through lcMethodLcmmGMM with random = ~ -1,
# which is the same hlme model but runs the N_STARTS gridsearch.
# Stability settings: shape = "linear" drops the quadratic term; ri_only = TRUE
# keeps only the random intercept (GMM).
lcmm_formulas <- function(gmm, shape = "quadratic", ri_only = FALSE) {
  quad <- shape == "quadratic"
  list(fixed   = if (quad) y ~ time + I(time^2) else y ~ time,
       mixture = if (quad) ~ time + I(time^2) else ~ time,
       random  = if (!gmm) ~ -1 else if (ri_only) ~ 1 else ~ time)
}

lcmm_fit <- function(d, k, seed, gmm, shape = "quadratic", ri_only = FALSE) {
  set.seed(seed)
  f <- lcmm_formulas(gmm, shape, ri_only)
  meth <- lcMethodLcmmGMM(fixed = f$fixed, mixture = f$mixture, random = f$random,
                          id = "id", time = "time", nClusters = k,
                          init = "gridsearch", gridsearch.rep = N_STARTS,
                          gridsearch.maxiter = 10, gridsearch.parallel = FALSE,
                          seed = seed)
  m <- latrend(meth, data = d$long, verbose = FALSE)
  list(cls = latrend_assign(m, d$ids), bic = BIC(m), aic = AIC(m),
       entropy = rel_entropy(postprob(m)), converged = isTRUE(converged(m) == 1),
       model = m)
}
fit_gbtm <- function(d, k, seed, ...) lcmm_fit(d, k, seed, gmm = FALSE, ...)
fit_gmm  <- function(d, k, seed, ...) lcmm_fit(d, k, seed, gmm = TRUE, ...)

# ---- Start recording (stability) -----------------------------------------------------
# The same procedure as lcmm::gridsearch (which latrend calls): a 1-class model
# for initial values, N_STARTS random starts of 10 iterations each, then the best
# start run to convergence. Here every start's log-likelihood and modal classes
# are kept, and the top `n_converge` starts (by short-run log-likelihood) are all
# run to convergence, so we can see whether the best log-likelihood replicates.
# hlme is called directly (numeric subject ids), not through latrend. hlme
# evaluates the argument of random() from its namespace, which falls through to
# the global environment, so the 1-class model is bound there for the duration
# (each forked worker has its own global environment).
lcmm_record <- function(d, k, seed, gmm, n_converge = 5, shape = "quadratic", ri_only = FALSE) {
  f <- lcmm_formulas(gmm, shape, ri_only)
  key <- d$ids
  dl <- d$long; dl$id <- match(dl$id, key)
  set.seed(seed)
  base <- list(fixed = f$fixed, random = f$random, subject = "id", data = dl, verbose = FALSE)
  minit <- do.call(lcmm::hlme, c(base, list(ng = 1)))
  if (minit$conv != 1) stop("1-class initial model did not converge")
  assign(".lcmm_minit", minit, envir = globalenv())
  on.exit(rm(".lcmm_minit", envir = globalenv()), add = TRUE)
  cls_of <- function(m) { p <- m$pprob; as.integer(p$class[match(seq_along(key), p$id)]) }

  short <- lapply(seq_len(N_STARTS), function(r) {
    m <- do.call(lcmm::hlme, c(base, list(mixture = f$mixture, ng = k, maxiter = 10,
                                         B = quote(random(.lcmm_minit)))))
    list(ll = m$loglik, best = m$best, cls = cls_of(m))
  })
  ll_short <- vapply(short, `[[`, numeric(1), "ll")
  top <- head(order(ll_short, decreasing = TRUE, na.last = NA), n_converge)
  full <- lapply(top, function(r) {
    m <- do.call(lcmm::hlme, c(base, list(mixture = f$mixture, ng = k, B = short[[r]]$best)))
    list(m = m, ll = m$loglik, conv = m$conv, cls = cls_of(m))
  })
  ll_full <- vapply(full, `[[`, numeric(1), "ll")
  conv_full <- vapply(full, `[[`, numeric(1), "conv")
  okc <- which(conv_full == 1)
  b <- if (length(okc)) okc[which.max(ll_full[okc])] else which.max(ll_full)
  m <- full[[b]]$m
  list(cls = full[[b]]$cls, bic = m$BIC, aic = m$AIC, entropy = rel_entropy(as.matrix(m$pprob[, -(1:2)])),
       converged = m$conv == 1, model = NULL,
       starts = tibble(start = seq_len(N_STARTS), ll_short = ll_short,
                       converged_run = match(seq_len(N_STARTS), top),
                       ll_full = ll_full[match(seq_len(N_STARTS), top)],
                       conv_full = conv_full[match(seq_len(N_STARTS), top)]),
       start_cls = sapply(short, `[[`, "cls"),
       full_cls = sapply(full, `[[`, "cls"))
}

# ---- GMM with every start converged (stability refit) -----------------------------------
# Decided 2026-11-07 after the GMM refit (README "GMM refit"): gridsearch's 10-iteration
# screen does not find the best optimum, so each of n_starts random starts gets 10
# iterations from random(minit), as gridsearch does, and is then run to convergence.
# The best converged start is kept; if none converges, the best start is returned with
# converged = FALSE. at_best counts converged starts within 0.5 log-likelihood of it.
# `model` is the hlme object, which predict_lcmm() accepts directly.
fit_gmm_conv <- function(d, k, seed, n_starts = 10, shape = "quadratic", ri_only = FALSE) {
  f <- lcmm_formulas(TRUE, shape, ri_only)
  key <- d$ids
  dl <- d$long; dl$id <- match(dl$id, key)
  set.seed(seed)
  base <- list(fixed = f$fixed, random = f$random, subject = "id", data = dl, verbose = FALSE)
  minit <- do.call(lcmm::hlme, c(base, list(ng = 1)))
  if (minit$conv != 1) stop("1-class initial model did not converge")
  assign(".lcmm_minit", minit, envir = globalenv())
  on.exit(rm(".lcmm_minit", envir = globalenv()), add = TRUE)
  runs <- lapply(seq_len(n_starts), function(s) tryCatch({
    short <- do.call(lcmm::hlme, c(base, list(mixture = f$mixture, ng = k, maxiter = 10,
                                             B = quote(random(.lcmm_minit)))))
    m <- do.call(lcmm::hlme, c(base, list(mixture = f$mixture, ng = k, B = short$best)))
    list(m = m, ll = m$loglik, conv = m$conv)
  }, error = function(e) list(m = NULL, ll = NA_real_, conv = NA_real_)))
  ll <- vapply(runs, `[[`, numeric(1), "ll"); conv <- vapply(runs, `[[`, numeric(1), "conv")
  ll[!is.finite(ll) | ll < -1e8] <- NA
  okc <- which(conv %in% 1)
  if (!length(okc) && all(is.na(ll))) stop("no GMM start produced a fit")
  b <- if (length(okc)) okc[which.max(ll[okc])] else which.max(ll)
  m <- runs[[b]]$m
  p <- m$pprob
  list(cls = as.integer(p$class[match(seq_along(key), p$id)]), bic = m$BIC, aic = m$AIC,
       entropy = rel_entropy(as.matrix(p[, -(1:2)])), converged = m$conv == 1, model = m,
       n_conv = length(okc), at_best = if (length(okc)) sum(ll[okc] > ll[b] - 0.5) else 0L)
}

# ---- Skew-aware GBTM (added 2026-11-13) ------------------------------------------------
# lcmm::lcmm with an estimated link function ("5-quant-splines"): the observed score is a
# monotone transformation of a latent Gaussian outcome, so the model estimates the
# transformation instead of assuming raw, T or log scores. No random effects (GBTM);
# class-specific quadratic in time. Starts as in fit_gmm_conv: n_starts random starts of
# 10 iterations from random(minit), each run to convergence, best converged kept. Tested
# in R/lcmm_link_pilot.R and R/lcmm_link_stability.R. Fitted on raw scores only (spec).
fit_gbtm_link <- function(d, k, seed, n_starts = 5) {
  key <- d$ids
  dl <- d$long; dl$id <- match(dl$id, key)
  set.seed(seed)
  base <- list(fixed = y ~ time + I(time^2), random = ~ -1, subject = "id", data = dl,
               link = "5-quant-splines", verbose = FALSE)
  minit <- do.call(lcmm::lcmm, base)
  if (minit$conv != 1) stop("1-class initial model did not converge")
  assign(".lcmm_minit", minit, envir = globalenv())
  on.exit(rm(".lcmm_minit", envir = globalenv()), add = TRUE)
  runs <- lapply(seq_len(n_starts), function(s) tryCatch({
    short <- do.call(lcmm::lcmm, c(base, list(mixture = ~ time + I(time^2), ng = k, maxiter = 10,
                                             B = quote(random(.lcmm_minit)))))
    m <- do.call(lcmm::lcmm, c(base, list(mixture = ~ time + I(time^2), ng = k, B = short$best)))
    list(m = m, ll = m$loglik, conv = m$conv)
  }, error = function(e) list(m = NULL, ll = NA_real_, conv = NA_real_)))
  ll <- vapply(runs, `[[`, numeric(1), "ll"); conv <- vapply(runs, `[[`, numeric(1), "conv")
  ll[!is.finite(ll) | ll < -1e8] <- NA
  okc <- which(conv %in% 1)
  if (!length(okc) && all(is.na(ll))) stop("no start produced a fit")
  b <- if (length(okc)) okc[which.max(ll[okc])] else which.max(ll)
  m <- runs[[b]]$m
  p <- m$pprob
  list(cls = as.integer(p$class[match(seq_along(key), p$id)]), bic = m$BIC, aic = m$AIC,
       entropy = rel_entropy(as.matrix(p[, -(1:2)])), converged = m$conv == 1, model = m,
       n_conv = length(okc), at_best = if (length(okc)) sum(ll[okc] > ll[b] - 0.5) else 0L)
}

# ---- One task ----------------------------------------------------------------------
fit_task <- function(task, wide, long_obs, pilot) {
  d <- task_data(task, wide, long_obs, pilot)
  fitter <- get(paste0("fit_", task$method))
  tree <- if (task$method == "ward") hclust(dist(d$X), method = "ward.D2") else NULL
  seed <- as.integer(SEED + match(task$task_id, task$all_ids[[1]]))

  purrr::map(task$ks[[1]], function(k) {
    t0 <- Sys.time()
    res <- tryCatch(
      if (task$method == "ward") fitter(d, k, seed, tree) else fitter(d, k, seed),
      error = function(e) list(error = conditionMessage(e)))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    base <- tibble(task_id = task$task_id, k = k, n = length(d$ids), seed = seed, secs = secs)
    if (!is.null(res$error)) return(list(fit = mutate(base, error = res$error), assign = NULL))
    props <- tabulate(res$cls, k) / length(res$cls)
    list(fit = mutate(base, error = NA_character_, converged = res$converged,
                      bic = res$bic, aic = res$aic, entropy = res$entropy,
                      asw = asw(d$X, res$cls, seed),
                      k_found = sum(props > 0), min_prop = min(props),
                      small_class = min(props) < CLASS_FLOOR,
                      n_conv = res$n_conv %||% NA_integer_, at_best = res$at_best %||% NA_integer_),
         assign = tibble(src_subject_id = d$ids, class = res$cls))
  })
}
