# =============================================================================
# R/gckm_pilot.R: pilot of growth-curve k-means (GCKM; latrend::lcMethodGCKM),
# added 2026-11-14. Exploratory, not part of the pipeline.
#   Rscript R/gckm_pilot.R run [cores]   # restartable, one file per run
#   Rscript R/gckm_pilot.R check         # held-out prediction reproduces in-sample BLUPs
#   Rscript R/gckm_pilot.R summary [cores]
# Method (latrend's GCKM defaults, written out so it can be checked and predicted):
#   1. lmer(y ~ time + (time | id)), time = age - AGE_CENTER, REML;
#   2. each child's coefficients = fixed + random (intercept, slope);
#   3. standardise the two columns (scale());
#   4. k-means with N_STARTS starts (latrend's default is 1), k = 2-6.
# Variant "quad": random intercept, slope and quadratic (3 coefficients), to see
# whether the curvature adds anything with 4 waves.
# Anchors: raw and log1p(raw), one child per family, >= 3 waves, imputed (N = 8,331).
# Runs per anchor and variant: the full data (reference), 20 80% subsamples (draws
# 1-20) and 20 null datasets (draws 1-20; the stability analysis uses 1-10; 11-20 added
# 2026-11-14 to check the raw-scale poverty false-positive rate).
# Stability rule for clusterBMA (proposed 2026-11-14): median pairwise subsample ARI
# >= .75 at every k on both anchors. Redundancy check: ARI of the full-data fit with
# the stability base solutions of k-means and kml (output/stability/base_assign.rds).
# Site split-half (added 2026-11-14, linear only): the same 20 splits by site as the
# stability analysis. Each half is fitted on its own; the other half's children are
# then assigned from the training half's model: BLUP coefficients from the training
# fixed effects and variance components, b_i = G Z_i' (Z_i G Z_i' + s^2 I)^-1 (y_i - X_i b),
# standardised with the training centre and scale, nearest training k-means centre.
# Metrics as in 09_stability_summary.R: ARI of the test half's own fit vs the
# predicted classes, and Tibshirani & Walther's prediction strength.
# Variant "linwave" (added 2026-11-14): linear, with wave number (0-3) as time instead of
# age, to test whether visit timing (follow-up span is slightly longer in poorer families)
# drives the raw-scale false positives in the null data.
# Poverty (added 2026-11-14): the primary race-adjusted multinomial model (POV_F$adj),
# highest vs lowest class by mean raw CBCL, on every full, subsample and null solution;
# compared with the stability base and subsample solutions of the other methods.
# Output: output/gckm/*.rds; output/summaries/summary_gckm_pilot.txt (aggregates only).
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(lme4); library(parallel); library(nnet) })

args <- commandArgs(trailingOnly = TRUE)
GDIR <- file.path(OUT, "gckm"); dir.create(GDIR, showWarnings = FALSE, recursive = TRUE)
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
N    <- length(base$ids)

RUNS <- tidyr::expand_grid(anchor = c("raw", "log"), variant = c("lin", "quad"),
                           bind_rows(tibble(ptype = "full", rep = 0L),
                                     tibble(ptype = "subsample", rep = 1:20),
                                     tibble(ptype = "null", rep = 1:20))) |>
  bind_rows(tidyr::expand_grid(anchor = c("raw", "log"), variant = c("lin", "linwave"), ptype = "site", rep = 1:20)) |>
  bind_rows(tidyr::expand_grid(anchor = c("raw", "log"), variant = "linwave",
                               bind_rows(tibble(ptype = "full", rep = 0L), tibble(ptype = "subsample", rep = 1:20),
                                         tibble(ptype = "null", rep = 1:20)))) |>
  mutate(file = file.path(GDIR, sprintf("%s_%s_%s_r%02d.rds", anchor, variant, ptype, rep)))

# One GCKM fit: lmer, child coefficients, standardised, k-means at every k. Returns the
# classes and what prediction needs.
fml_of <- function(variant) if (variant %in% c("lin", "linwave")) y ~ time + (time | id) else y ~ time + I(time^2) + (time + I(time^2) | id)
long_of <- function(X, A) data.frame(id = factor(rep(seq_len(nrow(X)), times = 4)), time = as.vector(A) - AGE_CENTER, y = as.vector(X))
gckm_fit <- function(X, A, variant) {
  n <- nrow(X)
  m <- suppressMessages(lmer(fml_of(variant), data = long_of(X, A), REML = TRUE))
  cf <- as.matrix(coef(m)$id)[as.character(seq_len(n)), , drop = FALSE]
  Z <- scale(cf)
  km <- lapply(K_STAB, \(k) { set.seed(SEED + k); kmeans(Z, k, nstart = N_STARTS, iter.max = 100) })
  cls <- sapply(km, `[[`, "cluster"); colnames(cls) <- paste0("k", K_STAB)
  list(cls = cls, m = m, cf = cf,
       pred = list(beta = fixef(m), G = as.matrix(VarCorr(m)$id), sigma = sigma(m), variant = variant,
                   center = attr(Z, "scaled:center"), scale = attr(Z, "scaled:scale"),
                   centers = lapply(km, `[[`, "centers")))
}
# BLUP coefficients for new children from a fitted model's fixed effects and variance components.
blup_coef <- function(pr, X, A) {
  t <- A - AGE_CENTER; q <- length(pr$beta)
  t(vapply(seq_len(nrow(X)), \(i) {
    Zi <- if (q == 2) cbind(1, t[i, ]) else cbind(1, t[i, ], t[i, ]^2)
    Vi <- Zi %*% pr$G %*% t(Zi) + diag(pr$sigma^2, ncol(X))
    pr$beta + drop(pr$G %*% t(Zi) %*% solve(Vi, X[i, ] - drop(Zi %*% pr$beta)))
  }, numeric(q)))
}
gckm_predict <- function(pr, X, A) {
  Z <- sweep(sweep(blup_coef(pr, X, A), 2, pr$center), 2, pr$scale, "/")
  sapply(pr$centers, \(C) apply(Z, 1, \(z) which.min(colSums((t(C) - z)^2))))
}

run_one <- function(i) {
  r <- RUNS[i, ]
  if (file.exists(r$file)) return(invisible(NULL))
  if (r$ptype == "site") return(run_site(r))
  idx <- if (r$ptype == "subsample") draw_subsample(N, r$rep) else seq_len(N)
  R <- if (r$ptype == "null") simulate_null(base, nm, r$rep) else base$R
  X <- R[idx, , drop = FALSE]; if (r$anchor == "log") X <- log1p(X)
  A <- base$ages[idx, , drop = FALSE]; n <- length(idx)
  if (r$variant == "linwave") A <- matrix(rep(0:3, each = n) + AGE_CENTER, n)   # wave number as time
  dl <- data.frame(id = factor(rep(seq_len(n), times = 4)), time = as.vector(A) - AGE_CENTER, y = as.vector(X))
  fml <- fml_of(r$variant)
  t0 <- Sys.time()
  m <- suppressMessages(lmer(fml, data = dl, REML = TRUE))
  cf <- as.matrix(coef(m)$id)[as.character(seq_len(n)), , drop = FALSE]
  Z <- scale(cf)
  cls <- sapply(K_STAB, \(k) { set.seed(SEED + k); kmeans(Z, k, nstart = N_STARTS, iter.max = 100)$cluster })
  colnames(cls) <- paste0("k", K_STAB)
  vc <- as.data.frame(VarCorr(m))
  out <- list(run = r, idx = idx, cls = cls, singular = isSingular(m), vc = vc[, c("grp", "var1", "var2", "sdcor")],
              cor_coef = cor(cf), mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, paste0(r$file, ".tmp")); file.rename(paste0(r$file, ".tmp"), r$file)
  cat(sprintf("%s %s %s %s r%02d done (%.1f min%s)\n", format(Sys.time(), "%m-%d %H:%M"), r$anchor, r$variant,
              r$ptype, r$rep, out$mins, if (out$singular) ", singular" else ""))
}

run_site <- function(r) {
  a <- draw_site_half(base$site, r$rep)
  X <- if (r$anchor == "log") log1p(base$R) else base$R
  AG <- if (r$variant == "linwave") matrix(rep(0:3, each = N) + AGE_CENTER, N) else base$ages
  t0 <- Sys.time()
  idx <- list(A = which(a), B = which(!a))
  f <- lapply(idx, \(i) gckm_fit(X[i, , drop = FALSE], AG[i, , drop = FALSE], r$variant))
  AonB <- gckm_predict(f$A$pred, X[idx$B, , drop = FALSE], AG[idx$B, , drop = FALSE])
  BonA <- gckm_predict(f$B$pred, X[idx$A, , drop = FALSE], AG[idx$A, , drop = FALSE])
  out <- list(run = r, idx = idx, A = f$A$cls, B = f$B$cls, AonB = AonB, BonA = BonA,
              mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(out, paste0(r$file, ".tmp")); file.rename(paste0(r$file, ".tmp"), r$file)
  cat(sprintf("%s %s %s site r%02d done (%.1f min)\n", format(Sys.time(), "%m-%d %H:%M"), r$anchor, r$variant, r$rep, out$mins))
}

# Check: predicting the training children from their own model must reproduce
# coef(m) (same BLUP formula) and the fitted classes.
if (length(args) >= 1 && args[1] == "check") {
  for (an in c("raw", "log")) for (va in c("lin", "quad")) {
    X <- if (an == "log") log1p(base$R) else base$R
    i <- draw_subsample(N, 1)[1:2000]
    f <- gckm_fit(X[i, ], base$ages[i, ], va)
    d <- max(abs(blup_coef(f$pred, X[i, ], base$ages[i, ]) - f$cf))
    p <- gckm_predict(f$pred, X[i, ], base$ages[i, ])
    cat(sprintf("%s %s: max |BLUP - coef(m)| = %.2e; class agreement k2-6: %s\n", an, va, d,
                paste(round(colMeans(p == f$cls), 4), collapse = " ")))
  }
  quit(save = "no")
}

if (length(args) >= 1 && args[1] == "run") {
  cores <- if (length(args) >= 2) as.integer(args[2]) else 1
  n_test <- if (length(args) >= 3) as.integer(args[3]) else nrow(RUNS)
  invisible(mclapply(seq_len(n_test), \(i) tryCatch(run_one(i), error = \(e) cat("run", i, "error:", conditionMessage(e), "\n")),
                     mc.cores = cores, mc.preschedule = FALSE))
  quit(save = "no")
}

# ---- Summary (aggregates only) ---------------------------------------------------------------
ari <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) < 2) NA_real_ else mclust::adjustedRandIndex(x[ok], y[ok]) }
fullv <- function(idx, v) { out <- rep(NA_integer_, N); out[idx] <- v; out }
cores <- if (length(args) >= 2) as.integer(args[2]) else 1
recs <- lapply(RUNS$file[file.exists(RUNS$file)], readRDS)
srecs <- Filter(\(x) x$run$ptype == "site", recs); recs <- Filter(\(x) x$run$ptype != "site", recs)
sumf <- file.path(SUMM, "summary_gckm_pilot.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(as.data.frame(x), row.names = FALSE)), collapse = "\n")
say(sprintf("Growth-curve k-means (GCKM) pilot (%s): %d of %d runs; median %.1f min per run; %d singular lmer fits.",
            format(Sys.time(), "%Y-%m-%d %H:%M"), length(recs), nrow(RUNS), median(sapply(recs, `[[`, "mins")),
            sum(sapply(recs, `[[`, "singular"))))

V <- bind_rows(lapply(recs, \(x) tibble(anchor = x$run$anchor, variant = x$run$variant, ptype = x$run$ptype, rep = x$run$rep,
                                        k = K_STAB, v = lapply(seq_along(K_STAB), \(j) fullv(x$idx, x$cls[, j])))))
wide <- \(d, col) d |> mutate(k = paste0("k", k)) |> tidyr::pivot_wider(names_from = k, values_from = all_of(col))

pairs <- V |> filter(ptype != "full") |> group_by(anchor, variant, ptype, k) |> group_modify(\(g, key) {
  if (nrow(g) < 2) return(tibble(ari = numeric()))
  ij <- t(combn(nrow(g), 2)); tibble(ari = apply(ij, 1, \(p) ari(g$v[[p[1]]], g$v[[p[2]]])))
}) |> ungroup()
say("\nMedian pairwise ARI between replicates (subsample: stability; null: same structure without classes):")
say(pt(pairs |> group_by(anchor, variant, ptype, k) |> summarise(a = round(median(ari, na.rm = TRUE), 2), .groups = "drop") |>
         wide("a") |> arrange(ptype, anchor, variant)))

ref <- V |> filter(ptype == "full") |> select(anchor, variant, k, ref = v)
say("\nMedian ARI of each subsample with the full-data fit:")
say(pt(V |> filter(ptype == "subsample") |> inner_join(ref, by = c("anchor", "variant", "k")) |>
         mutate(a = map2_dbl(v, ref, ari)) |> group_by(anchor, variant, k) |> summarise(a = round(median(a), 2), .groups = "drop") |> wide("a")))

say("\nSmallest class share (median over runs):")
say(pt(V |> mutate(m = map_dbl(v, \(x) { t <- table(x); min(t) / sum(t) })) |> group_by(anchor, variant, ptype, k) |>
         summarise(m = round(median(m), 3), .groups = "drop") |> wide("m") |> arrange(ptype, anchor, variant)))

ba <- readRDS(file.path(SOUT, "base_assign.rds"))
other <- bind_rows(lapply(c("kmeans", "kml", "gbtm", "pam"), \(meth) bind_rows(lapply(c("raw", "log"), \(an) bind_rows(lapply(K_STAB, \(k) {
  u <- ba$universes |> filter(anchor == an, method == meth, k == !!k)
  if (!nrow(u)) return(NULL)
  b <- ba$assign[[u$universe_id[1]]]
  tibble(method = meth, anchor = an, k = k, v = list(fullv(b$idx, b$cls)))
}))))))
say("\nRedundancy: ARI of the full-data GCKM fit with the stability base solutions (raw and log anchors):")
say(pt(ref |> inner_join(other, by = c("anchor", "k"), relationship = "many-to-many") |> mutate(a = round(map2_dbl(ref, v, ari), 2)) |>
         select(anchor, variant, method, k, a) |> wide("a") |> arrange(anchor, variant, method)))

say("\nRandom-effect SDs and correlation of the child coefficients, full-data fits:")
for (x in recs) if (x$run$ptype == "full") {
  say(sprintf("%s %s: %s; coefficient correlations %s", x$run$anchor, x$run$variant,
              paste(sprintf("%s %s%s %.3f", x$vc$grp, x$vc$var1, ifelse(is.na(x$vc$var2), "", paste0("/", x$vc$var2)), x$vc$sdcor), collapse = "; "),
              paste(round(x$cor_coef[upper.tri(x$cor_coef)], 2), collapse = " ")))
}

# ---- Site split-half ---------------------------------------------------------------------------
ps_tw <- function(test, pred) {           # Tibshirani & Walther (2005), as in 09_stability_summary.R
  ok <- !is.na(test) & !is.na(pred); test <- test[ok]; pred <- pred[ok]
  s <- vapply(split(pred, test), function(p) { n <- length(p); if (n < 2) return(NA_real_)
    m <- table(p); sum(m * (m - 1)) / (n * (n - 1)) }, numeric(1))
  if (all(is.na(s))) NA_real_ else min(s, na.rm = TRUE)
}
site <- bind_rows(lapply(srecs, \(x) bind_rows(lapply(seq_along(K_STAB), \(j) tibble(
  anchor = x$run$anchor, variant = x$run$variant, rep = x$run$rep, k = K_STAB[j], direction = c("A->B", "B->A"),
  ari = c(ari(x$B[, j], x$AonB[, j]), ari(x$A[, j], x$BonA[, j])),
  ps = c(ps_tw(x$B[, j], x$AonB[, j]), ps_tw(x$A[, j], x$BonA[, j])))))))
st <- readRDS(file.path(SOUT, "summary.rds"))
say(sprintf("\nSite split-half, linear GCKM, age and wave time (%d split runs): mean ARI (test fit vs predicted) and prediction strength; other methods from output/stability/summary.rds:", length(srecs)))
cmp_site <- bind_rows(site |> mutate(method = paste0("gckm_", variant)) |> select(-variant), st$site |> filter(method %in% c("kmeans", "kml", "pam", "gbtm"))) |>
  group_by(anchor, method, k) |> summarise(ari = round(mean(ari, na.rm = TRUE), 2), ps = round(mean(ps, na.rm = TRUE), 2), .groups = "drop")
say("ARI:"); say(pt(cmp_site |> select(-ps) |> wide("ari") |> arrange(anchor, method)))
say("Prediction strength:"); say(pt(cmp_site |> select(-ari) |> wide("ps") |> arrange(anchor, method)))

# ---- Poverty --------------------------------------------------------------------------------------
wd  <- readRDS(file.path(DER, "wide.rds")) |> add_race()
cov <- wd[match(base$ids, wd$src_subject_id), ] |>
  transmute(income_halving = -log2(income), age_bl, sex = factor(sex), site = factor(site), race)
poverty_one <- function(idx, cls, f = POV_F$adj) {      # as in 09_stability_summary.R
  mp <- relabel_map(cls, base$R[idx, , drop = FALSE]); cls <- unname(mp[as.character(cls)])
  dat <- cov[idx, ]; dat$class <- factor(cls)
  dat <- pov_rows(dat, f); dat$site <- droplevels(dat$site); dat$race <- droplevels(dat$race)
  K <- nlevels(dat$class); if (K < 2) return(NULL)
  top <- sum(dat$class == levels(dat$class)[K]); if (top < MIN_TOP_N) return(NULL)
  fit <- tryCatch(multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE), error = \(e) NULL)
  if (is.null(fit)) return(NULL)
  V <- tryCatch(vcov(fit), error = \(e) NULL); if (is.null(V)) return(NULL)
  term <- if (K == 2) "income_halving" else paste0(levels(dat$class)[K], ":income_halving")
  cf <- coef(fit); est <- if (K == 2) cf[["income_halving"]] else cf[levels(dat$class)[K], "income_halving"]
  se <- sqrt(max(V[term, term], 0))
  if (!is.finite(est) || !is.finite(se) || se <= 0 || se > MAX_SE) return(NULL)
  tibble(est = est, se = se, top_share = top / nrow(dat), top_mean = mean(rowMeans(base$R[idx, , drop = FALSE])[cls == K]))
}
jobs <- bind_rows(lapply(recs, \(x) tibble(anchor = x$run$anchor, variant = x$run$variant, ptype = x$run$ptype,
                                          rep = x$run$rep, k = K_STAB, idx = list(x$idx), cls = lapply(seq_along(K_STAB), \(j) x$cls[, j]))))
t0 <- Sys.time()
pv <- mclapply(seq_len(nrow(jobs)), \(i) poverty_one(jobs$idx[[i]], jobs$cls[[i]]), mc.cores = cores)
has <- !vapply(pv, is.null, logical(1))
pov <- bind_cols(jobs[has, c("anchor", "variant", "ptype", "rep", "k")], bind_rows(pv[has])) |>
  mutate(rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
say(sprintf("\nPoverty RRR (race-adjusted, highest vs lowest class per halving of income): %d of %d solutions estimable (%.1f min).",
            nrow(pov), nrow(jobs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
say("Full data: RRR (95% CI), size and mean raw CBCL of the highest class:")
say(pt(pov |> filter(ptype == "full") |> transmute(anchor, variant, k, rrr = sprintf("%.2f (%.2f-%.2f)", rrr, lo, hi),
                                                   top_pct = round(100 * top_share, 1), top_mean = round(top_mean, 1)) |> arrange(anchor, variant, k)))
say("Subsamples and null data: median RRR (share of solutions with CI above 1):")
say(pt(pov |> filter(ptype != "full") |> group_by(anchor, variant, ptype, k) |>
         summarise(v = sprintf("%.2f (%.2f)", median(rrr), mean(lo > 1)), .groups = "drop") |> wide("v") |> arrange(ptype, anchor, variant)))
say("Other methods, stability base solutions (same anchors): RRR")
say(pt(st$base |> filter(method %in% c("kmeans", "kml", "pam", "gbtm")) |> transmute(anchor, method, k, rrr = round(rrr, 2)) |>
         wide("rrr") |> arrange(anchor, method)))
