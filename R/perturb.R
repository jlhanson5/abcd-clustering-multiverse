# =============================================================================
# R/perturb.R: data for the stability analyses (see docs/STABILITY_PLAN.md).
# Anchor data, subsample / site-split / noise / null generators, settings
# variants, and the per-replicate data builder. Sourced by 07, 08 and 09.
#
# Paired design: every data perturbation is seeded by (perturbation, replicate)
# only, so all methods and both anchors see the same children, the same site
# halves, the same jitter and the same null datasets.
# =============================================================================

STAB_ANCHORS <- c("raw", "log")
CHEAP <- c("kmeans", "pam", "ward", "llpa", "kml")
TRAJ  <- c("gbtm", "gmm")
# Seed-dependent methods. Ward is deterministic. The Gaussian mixture is too at
# n <= 2,000, but at the anchor n mclust initialises on a random 2,000-child
# subset (see fit_llpa), so it gets seed replicates.
STOCH <- c("kmeans", "pam", "kml", "llpa")
K_STAB <- 2:6
SUB_FRAC <- 0.8
# Noise reliabilities. r = .90 replaced .94 on 2026-11-07: the school-age CBCL/6-18
# one-week test-retest is reported only as a range across scales, .80-.94 (Naar-King,
# Ellis & Frey, Assessing Children's Well-Being, 2004, p. 70), so .94 was its top;
# the preschool CBCL/1.5-5 Total Problems 8-day r is .90 (Achenbach & Rescorla 2000,
# Table 8-1). The noise90 jitter uses the old noise94 draws (same P_CODE), rescaled.
NOISE_R <- c(noise90 = 0.90, noise80 = 0.80)
N_CONVERGE <- 5                          # GBTM/GMM start recording: top starts run to convergence
GMM_CONV_STARTS <- 10                    # GMM stability refit (R/gmm_stab_refit.R): starts per fit

# Replicates per perturbation: cheap methods / GBTM-GMM.
REPS <- tribble(
  ~ptype,            ~cheap, ~traj,
  "seed",                50,     0,      # stochastic cheap methods only
  "starts",               0,     1,      # one recorded-start fit per anchor
  "subsample",          100,    20,
  "site",                20,    20,
  "noise90",             20,     5,
  "noise80",             20,     5,
  "null",                20,    10,
  "null_subsample",      20,     0,      # 20 subsamples of null dataset 1
  "setting",             20,     0)      # per variant: full-data fit (rep 0) + subsamples 1..20

# Settings variants (stability only). Arguments go straight to the fitters in
# R/fit_task.R, except `linkage`, which picks the hclust method.
VARIANTS <- tribble(
  ~method, ~variant,    ~args,
  "ward",  "average",   list(linkage = "average"),
  "ward",  "complete",  list(linkage = "complete"),
  "pam",   "manhattan", list(metric = "manhattan"),
  "llpa",  "VVI",       list(model_name = "VVI"),
  "llpa",  "EEE",       list(model_name = "EEE"),
  "gbtm",  "linear",    list(shape = "linear"),
  "gmm",   "ri_only",   list(ri_only = TRUE),
  "gmm",   "linear",    list(shape = "linear"))

# Seeds. Data seeds depend only on the perturbation and replicate.
P_CODE <- c(subsample = 1, site = 2, noise90 = 3, noise80 = 4, null = 5)
data_seed <- function(p, rep) as.integer(SEED + 1e5 * P_CODE[[p]] + rep)
fit_seed  <- function(ptype, rep) as.integer(SEED + 7e5 + 1e4 * match(ptype, REPS$ptype) + rep)

# ---- Anchor data ------------------------------------------------------------------
# One child per family, >= 3 waves with the shared imputation (ge3_imp). `R` is
# the imputed raw CBCL matrix; the log anchor is log1p(R).
anchor_base <- function(wide, smoke_n = NULL) {
  w <- filter(wide, one_per_family) |> arrange(src_subject_id)
  if (!is.null(smoke_n)) { set.seed(SEED); w <- w[sort(sample(nrow(w), smoke_n)), ] }
  R <- as.matrix(w[paste0("imp_r_", WAVES)]); dimnames(R) <- list(w$src_subject_id, WAVES)
  ages <- as.matrix(w[paste0("age_", WAVES)])
  list(ids = w$src_subject_id, R = R, ages = ages, site = w$site,
       sd_w = apply(R, 2, sd))
}

# ---- Generators ---------------------------------------------------------------------
draw_subsample <- function(n, rep) {
  set.seed(data_seed("subsample", rep)); sort(sample(n, round(SUB_FRAC * n)))
}

# Random half of the sites (half A). With 22 sites, 11 vs 11.
draw_site_half <- function(site, rep) {
  s <- sort(unique(site))
  set.seed(data_seed("site", rep)); a <- sample(s, floor(length(s) / 2))
  site %in% a
}

# Measurement noise on the raw scale: SEM = SD_wave * sqrt(1 - r). Jittered raw
# scores are rounded and floored at 0 (the log anchor takes log1p of these).
jitter_raw <- function(R, sd_w, ptype, rep) {
  sem <- sd_w * sqrt(1 - NOISE_R[[ptype]])
  set.seed(data_seed(ptype, rep))
  E <- matrix(rnorm(length(R)), nrow(R)) %*% diag(sem)
  out <- pmax(round(R + E), 0); dimnames(out) <- dimnames(R); out
}

# Null generating model: one-class growth model on log(1+raw), quadratic in
# age - 9, random intercept + slope. Fitted once in 07 (aggregate parameters only).
fit_null_model <- function(base) {
  n <- length(base$ids)
  dl <- data.frame(id = rep(seq_len(n), times = 4), time = as.vector(base$ages) - AGE_CENTER,
                   y = as.vector(log1p(base$R)))
  m <- lcmm::hlme(y ~ time + I(time^2), random = ~ time, subject = "id", ng = 1,
                  data = dl, verbose = FALSE)
  if (m$conv != 1) stop("null model did not converge")
  b <- m$best
  U <- matrix(0, 2, 2); U[upper.tri(U, diag = TRUE)] <- m$cholesky
  list(beta = unname(b[c("intercept", "time", "I(time^2)")]),
       Sigma = crossprod(U), sigma = abs(unname(b[["stderr"]])), loglik = m$loglik)
}

# A null dataset: the real children's ages (so sex, site and income stay theirs),
# CBCL drawn from the 1-class model. Returns the raw matrix: expm1 of the log
# draw, rounded and floored at 0. The log anchor then takes log1p, which differs
# from the simulated log score only by that rounding.
simulate_null <- function(base, nm, rep) {
  set.seed(data_seed("null", rep))
  n <- length(base$ids); t <- base$ages - AGE_CENTER
  b <- matrix(rnorm(2 * n), n) %*% chol(nm$Sigma)
  Y <- nm$beta[1] + b[, 1] + (nm$beta[2] + b[, 2]) * t + nm$beta[3] * t^2 +
       matrix(rnorm(4 * n, 0, nm$sigma), n)
  out <- pmax(round(expm1(Y)), 0); dimnames(out) <- dimnames(base$R); out
}

# ---- Per-replicate data ----------------------------------------------------------------
# Returns a named list of parts (one, or A and B for site splits). Each part has
# idx (rows of the anchor), R (raw matrix the replicate clusters) and ages.
rep_data <- function(rep, base, nm) {
  n <- length(base$ids)
  part <- function(idx, R = base$R) list(idx = idx, R = R[idx, , drop = FALSE],
                                         ages = base$ages[idx, , drop = FALSE])
  switch(rep$ptype,
    seed = , starts = list(all = part(seq_len(n))),
    subsample = list(all = part(draw_subsample(n, rep$rep))),
    setting = if (rep$rep == 0) list(all = part(seq_len(n))) else
                list(all = part(draw_subsample(n, rep$rep))),
    site = { a <- draw_site_half(base$site, rep$rep)
             list(A = part(which(a)), B = part(which(!a))) },
    noise90 = , noise80 = list(all = part(seq_len(n), jitter_raw(base$R, base$sd_w, rep$ptype, rep$rep))),
    null = list(all = part(seq_len(n), simulate_null(base, nm, rep$rep))),
    null_subsample = list(all = part(draw_subsample(n, rep$rep), simulate_null(base, nm, 1))),
    stop("unknown ptype ", rep$ptype))
}

# The same structure task_data() builds in R/fit_task.R: X on the anchor's scale,
# and long data with time = age - 9 for GBTM/GMM, wave position for the others.
make_d <- function(p, scale, method, ids_all) {
  ids <- ids_all[p$idx]
  X <- if (scale == "log") log1p(p$R) else p$R
  dimnames(X) <- list(ids, WAVES)
  n <- length(ids)
  time <- if (method %in% TRAJ) as.vector(p$ages) - AGE_CENTER else rep(0:3, each = n)
  lg <- data.frame(id = rep(ids, times = 4), time = time, y = as.vector(X))
  list(X = X, long = lg[order(lg$id, lg$time), ], ids = ids)
}

# Class order: 1..K by mean raw CBCL across the 4 waves of the data the replicate
# clustered (the 04_poverty.R rule). Empty classes drop out. Returns the map from
# fitted label to relabelled class, so predictions can use the same map.
relabel_map <- function(cls, R) {
  m <- tapply(rowMeans(R), cls, mean)
  setNames(rank(m, ties.method = "first"), names(m))
}
