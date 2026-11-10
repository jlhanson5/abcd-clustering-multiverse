# =============================================================================
# 09_stability_summary.R: stability metrics from the replicates in 08.
#   Rscript 09_stability_summary.R full|smoke [cores]
# Output: output/stability/summary.rds (tables below) and
#         output/summaries/summary_stability_<mode>.txt
#   fits      fit indices, one row per replicate x k x part
#   pairs     ARI between replicate solutions (shared children), plus each
#             replicate vs the primary run's (full_refit) base solution (j = 0)
#   site      prediction strength per split and direction: ARI of the test half's
#             own fit vs the training half's predicted classes, and Tibshirani &
#             Walther's prediction strength (min over test classes of the share of
#             co-assigned pairs kept)
#   child     child-level consistency: per child, the share of replicates (that
#             include the child) placing them in the highest class; summarised as
#             reliably high (>= .90), never high (<= .10), boundary
#   child_hist  binned distribution of that share (Stab-2)
#   jaccard   per base class, the best-matching replicate class's Jaccard
#   profiles  class sizes and mean raw CBCL per wave, every solution
#   ksel      k picked by silhouette (max) / BIC (min), per replicate
#   poverty   RRR (highest vs lowest class, per halving of income) per solution,
#             race-adjusted (primary since 2026-11-09); poverty_unadj: age/sex/site only
#   starts    GBTM/GMM start recording: per-start log-likelihoods, replication
#   settings  each settings variant vs base: full-data ARI; subsample stability
#   base      the primary run's (full_refit) anchor fits (indices + poverty) for reference
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(nnet); library(parallel) })

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "full"
cores <- if (length(args) >= 2) as.integer(args[2]) else 8
SOUT  <- file.path(OUT, if (mode == "smoke") "stability_smoke" else "stability")
sumf  <- file.path(SUMM, paste0("summary_stability_", mode, ".txt")); unlink(sumf)
say   <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
t_all <- Sys.time()

base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
ba   <- readRDS(file.path(SOUT, "base_assign.rds"))
N    <- length(base$ids)
REAL <- c("seed", "subsample", "noise90", "noise80", "site")

# ---- Load replicates -----------------------------------------------------------------
files <- list.files(file.path(SOUT, "reps"), pattern = "\\.rds$", full.names = TRUE)
recs  <- mclapply(files, readRDS, mc.cores = cores)
meta  <- bind_rows(lapply(recs, \(r) select(r$rep, rep_id, anchor, method, ptype, variant, rep)))
fits  <- bind_rows(lapply(recs, \(r) bind_rows(lapply(r$fits, `[[`, "fit")))) |>
  left_join(meta, by = "rep_id")
sol <- bind_rows(lapply(recs, \(r) bind_rows(lapply(r$fits, \(o) {
  if (!length(o$parts)) return(NULL)
  tibble(rep_id = r$rep$rep_id, k = o$fit$k[1], part = names(o$parts),
         idx = lapply(o$parts, `[[`, "idx"), cls = lapply(o$parts, `[[`, "cls"))
})))) |> left_join(meta, by = "rep_id")
preds <- bind_rows(lapply(recs, \(r) bind_rows(lapply(r$fits, \(o) {
  if (is.null(o$pred)) return(NULL)
  tibble(rep_id = r$rep$rep_id, k = o$fit$k[1], AonB = list(o$pred$AonB), BonA = list(o$pred$BonA))
}))))
startrec <- lapply(recs[meta$ptype == "starts"], \(r) list(rep = r$rep, fits = r$fits))
rm(recs); invisible(gc())
say(sprintf("Stability summary, %s (%s). %d replicates, %d solutions, %d fit errors.",
            mode, format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(meta), nrow(sol), sum(!is.na(fits$error))))

fullv <- function(idx, cls) { v <- rep(NA_integer_, N); v[idx] <- cls; v }
base_vec <- function(anchor, method, k) {
  u <- sprintf("%s_%s_%s_one_ge3_imp_k%d", method, anchor, ifelse(method %in% TRAJ, "age", "wave"), k)
  b <- ba$assign[[u]]; if (is.null(b)) NULL else fullv(b$idx, b$cls)
}

# ARI on the children both solutions include (05_agreement.R's formula).
ari <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  n <- length(x); if (n < 2) return(c(NA, n))
  x <- match(x, sort(unique(x))); y <- match(y, sort(unique(y)))
  kx <- max(x); ky <- max(y)
  tab <- matrix(tabulate((x - 1L) * ky + y, kx * ky), kx, ky, byrow = TRUE)
  c2 <- function(v) sum(v * (v - 1) / 2)
  sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); ex <- sa * sb / (n * (n - 1) / 2)
  den <- (sa + sb) / 2 - ex
  c(if (den == 0) NA else (c2(tab) - ex) / den, n)
}

# ---- Pairwise ARI ------------------------------------------------------------------------
pair_groups <- sol |>
  filter(part == "all", ptype %in% c("seed", "subsample", "noise90", "noise80", "null_subsample", "setting"),
         !(ptype == "setting" & rep == 0)) |>
  group_by(anchor, method, variant, ptype, k) |> group_split()
pairs <- bind_rows(mclapply(pair_groups, function(g) {
  V <- do.call(cbind, mapply(fullv, g$idx, g$cls, SIMPLIFY = FALSE))
  r <- g$rep; R <- length(r)
  ij <- if (R > 1) t(combn(R, 2)) else matrix(integer(0), 0, 2)
  out <- matrix(vapply(seq_len(nrow(ij)), \(p) ari(V[, ij[p, 1]], V[, ij[p, 2]]), numeric(2)), nrow = 2)
  pr <- tibble(i = r[ij[, 1]], j = r[ij[, 2]], ari = out[1, ], n_common = out[2, ])
  b <- if (!g$ptype[1] %in% c("null_subsample")) base_vec(g$anchor[1], g$method[1], g$k[1])
  if (!is.null(b) && g$variant[1] == "base") {
    vb <- matrix(vapply(seq_len(R), \(c) ari(V[, c], b), numeric(2)), nrow = 2)
    pr <- bind_rows(pr, tibble(i = r, j = 0L, ari = vb[1, ], n_common = vb[2, ]))
  }
  bind_cols(g[1, c("anchor", "method", "variant", "ptype", "k")], pr)
}, mc.cores = cores))

# ---- Site splits: prediction strength ------------------------------------------------------
ps_tw <- function(test, pred) {           # Tibshirani & Walther (2005)
  ok <- !is.na(test) & !is.na(pred); test <- test[ok]; pred <- pred[ok]
  s <- vapply(split(pred, test), function(p) {
    n <- length(p); if (n < 2) return(NA_real_)
    m <- table(p); sum(m * (m - 1)) / (n * (n - 1))
  }, numeric(1))
  # All NA when the prediction failed (3 of 2,800 directions in the full run:
  # degenerate raw GBTM/GMM training fits at k 5-6, hlme errors on predict).
  if (all(is.na(s))) NA_real_ else min(s, na.rm = TRUE)
}
site <- sol |> filter(ptype == "site") |>
  select(rep_id, anchor, method, rep, k, part, cls) |>
  pivot_wider(names_from = part, values_from = cls) |>
  inner_join(preds, by = c("rep_id", "k"))
site <- bind_rows(lapply(seq_len(nrow(site)), function(i) {
  s <- site[i, ]
  if (is.null(s$A[[1]]) || is.null(s$B[[1]])) return(NULL)
  tibble(anchor = s$anchor, method = s$method, rep = s$rep, k = s$k, direction = c("A->B", "B->A"),
         ari = c(ari(s$B[[1]], s$AonB[[1]])[1], ari(s$A[[1]], s$BonA[[1]])[1]),
         ps = c(ps_tw(s$B[[1]], s$AonB[[1]]), ps_tw(s$A[[1]], s$BonA[[1]])))
}))

# ---- Child-level consistency -------------------------------------------------------------------
child_one <- function(g) {
  hi <- den <- numeric(N)
  for (i in seq_len(nrow(g))) {
    idx <- g$idx[[i]]; cl <- g$cls[[i]]
    den[idx] <- den[idx] + 1; hi[idx] <- hi[idx] + (cl == max(cl))
  }
  p <- ifelse(den > 0, hi / den, NA)
  b <- base_vec(g$anchor[1], g$method[1], g$k[1])
  base_top <- if (is.null(b)) rep(NA, N) else b == max(b, na.rm = TRUE)
  ok <- !is.na(p)
  list(sum = tibble(n_children = sum(ok), n_reps = nrow(g), mean_reps_per_child = mean(den[ok]),
                    reliably_high = mean(p[ok] >= .9), never_high = mean(p[ok] <= .1),
                    boundary = mean(p[ok] > .1 & p[ok] < .9), ever_high = mean(p[ok] > 0),
                    base_top_n = sum(base_top & ok, na.rm = TRUE),
                    base_top_reliable = mean(p[ok & base_top %in% TRUE] >= .9)),
       hist = tibble(bin = cut(p[ok], c(-Inf, 0, seq(.05, 1, .05)), labels = FALSE)) |> count(bin))
}
cgroups <- sol |>
  filter(variant == "base", ptype %in% REAL) |>
  bind_rows(mutate(filter(sol, variant == "base", ptype %in% REAL), ptype = "pooled")) |>
  group_by(anchor, method, ptype, k) |> group_split()
cres <- mclapply(cgroups, child_one, mc.cores = cores)
key  <- bind_rows(lapply(cgroups, \(g) g[1, c("anchor", "method", "ptype", "k")]))
child <- bind_cols(key, bind_rows(lapply(cres, `[[`, "sum")))
child_hist <- bind_rows(lapply(seq_along(cres), \(i) bind_cols(key[i, ], cres[[i]]$hist)))

# ---- Jaccard vs base, per base class -------------------------------------------------------------
jac <- function(b, r) {
  ok <- !is.na(b) & !is.na(r); b <- b[ok]; r <- r[ok]
  tab <- table(b, r)
  J <- tab / (outer(rowSums(tab), colSums(tab), "+") - tab)
  tibble(class = as.integer(rownames(tab)), jaccard = apply(J, 1, max))
}
jaccard <- bind_rows(mclapply(group_split(group_by(filter(sol, variant == "base", ptype %in% REAL),
                                                   anchor, method, k)), function(g) {
  b <- base_vec(g$anchor[1], g$method[1], g$k[1]); if (is.null(b)) return(NULL)
  bind_rows(lapply(seq_len(nrow(g)), \(i)
    bind_cols(g[i, c("anchor", "method", "ptype", "rep", "k", "part")], jac(b, fullv(g$idx[[i]], g$cls[[i]])))))
}, mc.cores = cores))

# ---- Class profiles (sizes, mean raw CBCL per wave) ------------------------------------------------
# Noise and null replicates cluster their own raw matrix; regenerate it (cheap, seeded).
Rcache <- list()
R_of <- function(ptype, rep) {
  key <- switch(ptype, noise90 = , noise80 = paste(ptype, rep), null = paste("null", rep),
                null_subsample = "null 1", "real")
  if (is.null(Rcache[[key]])) Rcache[[key]] <<- switch(key, real = base$R,
    if (startsWith(key, "null")) simulate_null(base, nm, as.integer(sub("null ", "", key)))
    else jitter_raw(base$R, base$sd_w, ptype, rep))
  Rcache[[key]]
}
profiles <- bind_rows(lapply(seq_len(nrow(sol)), function(i) {
  s <- sol[i, ]; R <- R_of(s$ptype, s$rep)[s$idx[[1]], , drop = FALSE]; cl <- s$cls[[1]]
  M <- rowsum(R, cl) / as.vector(table(cl))
  tibble(rep_id = s$rep_id, k = s$k, part = s$part, class = as.integer(rownames(M)),
         K = max(cl), n = as.vector(table(cl)), prop = n / length(cl),
         m_BL = M[, 1], m_Y1 = M[, 2], m_Y2 = M[, 3], m_Y3 = M[, 4], mean = rowMeans(M))
})) |> left_join(meta, by = "rep_id")
rm(Rcache)

# ---- k selection -----------------------------------------------------------------------------------
ksel <- fits |> filter(is.na(error), k_found == k) |>
  group_by(rep_id, anchor, method, ptype, variant, rep, part) |>
  summarise(k_asw = k[which.max(asw)],
            k_bic = if (all(is.na(bic))) NA_integer_ else k[which.min(bic)], .groups = "drop")

# ---- Poverty model per solution ------------------------------------------------------------------------
# Primary: race-adjusted (POV_F$adj, 00_setup.R); poverty_unadj: age/sex/site only.
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()
cov  <- wide[match(base$ids, wide$src_subject_id), ] |>
  transmute(income_halving = -log2(income), age_bl, sex = factor(sex), site = factor(site), race)
poverty_one <- function(idx, cls, f) {
  dat <- cov[idx, ]; dat$class <- factor(cls)
  dat <- pov_rows(dat, f); dat$site <- droplevels(dat$site); dat$race <- droplevels(dat$race)
  K <- nlevels(dat$class); if (K < 2) return(NULL)
  if (sum(dat$class == levels(dat$class)[K]) < MIN_TOP_N) return(NULL)   # not estimable (00_setup.R)
  fit <- tryCatch(multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  V <- tryCatch(vcov(fit), error = function(e) NULL); if (is.null(V)) return(NULL)
  term <- if (K == 2) "income_halving" else paste0(levels(dat$class)[K], ":income_halving")
  cf <- coef(fit); est <- if (K == 2) cf[["income_halving"]] else cf[levels(dat$class)[K], "income_halving"]
  se <- sqrt(max(V[term, term], 0))
  if (!is.finite(est) || !is.finite(se) || se <= 0 || se > MAX_SE) return(NULL)   # separation: not estimable (00_setup.R)
  tibble(K = K, est = est, se = se, mn_converged = fit$convergence == 0)
}
poverty_all <- function(f) {
  t0 <- Sys.time()
  pov <- mclapply(seq_len(nrow(sol)), \(i) poverty_one(sol$idx[[i]], sol$cls[[i]], f), mc.cores = cores)
  has <- !vapply(pov, is.null, logical(1))
  out <- bind_cols(sol[has, c("rep_id", "anchor", "method", "ptype", "variant", "rep", "k", "part")],
                   bind_rows(pov[has])) |>
    mutate(rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  say(sprintf("Poverty models (%s): %d of %d solutions (%.1f min).", deparse(f[[3]]), nrow(out), nrow(sol),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out
}
poverty <- poverty_all(POV_F$adj)
poverty_unadj <- poverty_all(POV_F$unadj)

# ---- Start recording ----------------------------------------------------------------------------------
starts <- bind_rows(lapply(startrec, function(r) bind_rows(lapply(r$fits, function(o) {
  s <- o$starts; if (is.null(s)) return(NULL)
  # Columns of full_cls are the converged runs, in converged_run order.
  run <- s$starts |> filter(!is.na(converged_run)) |> arrange(converged_run)
  okr <- run$conv_full %in% 1
  ib <- if (any(okr)) which(okr)[which.max(run$ll_full[okr])] else which.max(run$ll_full)
  best <- run$ll_full[ib]
  fc <- s$full_cls
  ari_full  <- apply(fc, 2, \(v) ari(v, fc[, ib])[1])
  ari_short <- apply(s$start_cls, 2, \(v) ari(v, fc[, ib])[1])
  b <- base_vec(r$rep$anchor, r$rep$method, o$fit$k[1])
  s$starts |> mutate(anchor = r$rep$anchor, method = r$rep$method, k = o$fit$k[1],
                     best_ll = best, at_best = abs(ll_full - best) < 0.01,
                     ari_full_vs_best = ari_full[converged_run],
                     ari_short_vs_best = ari_short,
                     best_vs_base = if (is.null(b)) NA else ari(fullv(seq_len(N), fc[, ib]), b)[1])
}))))

# ---- Settings variants ---------------------------------------------------------------------------------
var_full <- sol |> filter(ptype == "setting", rep == 0) |>
  mutate(ari_vs_base = map2_dbl(seq_along(rep_id), k, \(i, kk) {
    b <- base_vec(anchor[i], method[i], kk); if (is.null(b)) NA else ari(fullv(idx[[i]], cls[[i]]), b)[1]
  })) |> select(anchor, method, variant, k, ari_vs_base)
sub_stab <- pairs |> filter(j > 0, i <= 20, j <= 20,
                            (ptype == "setting") | (ptype == "subsample" & variant == "base")) |>
  group_by(anchor, method, variant, k) |>
  summarise(sub_ari_mean = mean(ari, na.rm = TRUE), n_pairs = n(), .groups = "drop")
settings <- full_join(var_full, sub_stab, by = c("anchor", "method", "variant", "k"))

# ---- Base (primary-run, full_refit) anchor fits, for reference --------------------------------------------------------
res <- readRDS(file.path(RES, "results_full_refit.rds"))   # primary run since 2026-11-07
basefits <- res$fits |> filter(universe_id %in% ba$universes$universe_id) |>
  left_join(select(ba$universes, universe_id, anchor), by = "universe_id") |>
  select(universe_id, anchor, method, k, bic, entropy, asw, k_found, min_prop, est, se, rrr, lo, hi, top_n, top_mean,
         rrr_unadj, lo_unadj, hi_unadj)
rm(res); invisible(gc())

saveRDS(list(fits = fits, pairs = pairs, site = site, child = child, child_hist = child_hist,
             jaccard = jaccard, profiles = profiles, ksel = ksel, poverty = poverty, poverty_unadj = poverty_unadj,
             starts = starts, settings = settings, base = basefits, meta = meta),
        file.path(SOUT, "summary.rds"))

# ---- Text summary -----------------------------------------------------------------------------------------
pt <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")
say("\nMean ARI between replicates (pairs among replicates; vs base = each replicate vs the primary-run solution):")
say(pt(pairs |> filter(variant == "base") |>
  group_by(anchor, ptype, method, k) |>
  summarise(pairs = mean(ari[j > 0], na.rm = TRUE), vs_base = mean(ari[j == 0], na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = k, values_from = c(pairs, vs_base), names_glue = "k{k}_{.value}") |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))))
say("\nSite split-half prediction strength (mean over splits and directions):")
say(pt(site |> group_by(anchor, method, k) |>
  summarise(ari = round(mean(ari, na.rm = TRUE), 2), ps = round(mean(ps, na.rm = TRUE), 2), .groups = "drop") |>
  pivot_wider(names_from = k, values_from = c(ari, ps), names_glue = "k{k}_{.value}")))
say("\nChild-level consistency, all real-data perturbations pooled (share of children):")
say(pt(child |> filter(ptype == "pooled") |>
  transmute(anchor, method, k, reliably_high = round(reliably_high, 3), boundary = round(boundary, 3),
            never_high = round(never_high, 3), base_top_reliable = round(base_top_reliable, 2))))
pov_tab <- function(p) pt(p |> filter(variant == "base", ptype != "starts") |>
  group_by(anchor, ptype, method) |>
  summarise(n = n(), median = round(median(rrr), 2), q25 = round(quantile(rrr, .25), 2),
            q75 = round(quantile(rrr, .75), 2), ci_above_1 = round(mean(lo > 1), 2),
            ci_excl_1 = round(mean(lo > 1 | hi < 1), 2), .groups = "drop"))
say("\nPoverty RRR across replicates (base variant), race-adjusted (primary): median, IQR, share with CI above 1:")
say(pov_tab(poverty))
say("\nSame, without race adjustment (sensitivity; age, sex, site only):")
say(pov_tab(poverty_unadj))
say("\nNull benchmark vs real data (mean silhouette, relative entropy):")
say(pt(bind_rows(
  basefits |> transmute(data = "real (base)", anchor, method, k, asw, entropy),
  fits |> filter(ptype == "null", is.na(error)) |> transmute(data = "null", anchor, method, k, asw, entropy)) |>
  group_by(data, anchor, method) |>
  summarise(asw = round(mean(asw, na.rm = TRUE), 3), entropy = round(mean(entropy, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = data, values_from = c(asw, entropy))))
if (nrow(starts)) {
  say("\nGBTM/GMM start recording: converged runs reaching the best log-likelihood (of those run):")
  say(pt(starts |> filter(!is.na(converged_run)) |> group_by(anchor, method, k) |>
    summarise(runs = n(), at_best = sum(at_best, na.rm = TRUE),
              short_starts_ari_mean = round(mean(ari_short_vs_best, na.rm = TRUE), 2),
              best_vs_base = round(first(best_vs_base), 2), .groups = "drop")))
}
say("\nSettings variants: ARI with base (full data); subsample stability, variant vs base method (reps 1-20):")
say(pt(settings |> mutate(across(where(is.double), \(x) round(x, 2)))))
say(sprintf("\nDone in %.1f min.", as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
