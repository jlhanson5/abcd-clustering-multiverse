# =============================================================================
# R/gckm_pilot.R: pilot of growth-curve k-means (GCKM; latrend::lcMethodGCKM),
# added 2026-11-14. Exploratory, not part of the pipeline.
#   Rscript R/gckm_pilot.R run [cores]   # restartable, one file per run
#   Rscript R/gckm_pilot.R summary
# Method (latrend's GCKM defaults, written out so it can be checked and predicted):
#   1. lmer(y ~ time + (time | id)), time = age - AGE_CENTER, REML;
#   2. each child's coefficients = fixed + random (intercept, slope);
#   3. standardise the two columns (scale());
#   4. k-means with N_STARTS starts (latrend's default is 1), k = 2-6.
# Variant "quad": random intercept, slope and quadratic (3 coefficients), to see
# whether the curvature adds anything with 4 waves.
# Anchors: raw and log1p(raw), one child per family, >= 3 waves, imputed (N = 8,331).
# Runs per anchor and variant: the full data (reference), 20 80% subsamples (draws
# 1-20) and 5 null datasets (draws 1-5), the same draws as the stability analysis.
# Stability rule for clusterBMA (proposed 2026-11-14): median pairwise subsample ARI
# >= .75 at every k on both anchors. Redundancy check: ARI of the full-data fit with
# the stability base solutions of k-means and kml (output/stability/base_assign.rds).
# Output: output/gckm/*.rds; output/summaries/summary_gckm_pilot.txt (aggregates only).
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages({ library(lme4); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
GDIR <- file.path(OUT, "gckm"); dir.create(GDIR, showWarnings = FALSE, recursive = TRUE)
SOUT <- file.path(OUT, "stability")
base <- readRDS(file.path(SOUT, "anchor.rds"))
nm   <- readRDS(file.path(SOUT, "null_model.rds"))
N    <- length(base$ids)

RUNS <- tidyr::expand_grid(anchor = c("raw", "log"), variant = c("lin", "quad"),
                           bind_rows(tibble(ptype = "full", rep = 0L),
                                     tibble(ptype = "subsample", rep = 1:20),
                                     tibble(ptype = "null", rep = 1:5))) |>
  mutate(file = file.path(GDIR, sprintf("%s_%s_%s_r%02d.rds", anchor, variant, ptype, rep)))

run_one <- function(i) {
  r <- RUNS[i, ]
  if (file.exists(r$file)) return(invisible(NULL))
  idx <- if (r$ptype == "subsample") draw_subsample(N, r$rep) else seq_len(N)
  R <- if (r$ptype == "null") simulate_null(base, nm, r$rep) else base$R
  X <- R[idx, , drop = FALSE]; if (r$anchor == "log") X <- log1p(X)
  A <- base$ages[idx, , drop = FALSE]; n <- length(idx)
  dl <- data.frame(id = factor(rep(seq_len(n), times = 4)), time = as.vector(A) - AGE_CENTER, y = as.vector(X))
  fml <- if (r$variant == "lin") y ~ time + (time | id) else y ~ time + I(time^2) + (time + I(time^2) | id)
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
recs <- lapply(RUNS$file[file.exists(RUNS$file)], readRDS)
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
