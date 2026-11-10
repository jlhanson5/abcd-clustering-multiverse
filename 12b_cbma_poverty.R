# =============================================================================
# 12b_cbma_poverty.R: poverty models on the stored clusterBMA consensus classes
# (output/cbma/*.rds from 12_cbma.R), race-adjusted (primary since 2026-11-09)
# and without race (age, sex, site; the model 12_cbma.R stored). No reclustering.
#   Rscript 12b_cbma_poverty.R          # after 12_cbma.R; before 13 and 14
# Every class vs the lowest (class 1), model-based SEs (one child per family).
# Checks that the unadjusted model reproduces 12_cbma.R's stored RRR.
# Output: output/results/cbma_poverty.rds, one row per fit x model x non-reference class
#   (top = highest class).
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(nnet))

wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()
files <- sort(list.files(file.path(OUT, "cbma"), "\\.rds$", full.names = TRUE))

pov_fit <- function(r, model) {
  f <- POV_F[[model]]
  w <- wide[match(r$ids, wide$src_subject_id), ]
  dat <- tibble(class = factor(r$labels), income_halving = -log2(w$income), age_bl = w$age_bl,
                sex = factor(w$sex), site = droplevels(factor(w$site)), race = w$race,
                cbcl_mean = rowMeans(as.matrix(w[paste0("imp_r_", WAVES)])))
  dat <- pov_rows(dat, f); dat$race <- droplevels(dat$race)
  m <- multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE)
  cf <- coef(m); V <- vcov(m)
  if (is.null(dim(cf))) cf <- matrix(cf, 1, dimnames = list(levels(dat$class)[2], names(cf)))
  cls <- rownames(cf)
  terms <- if (length(cls) == 1) "income_halving" else paste0(cls, ":income_halving")
  tibble(scale = r$scale, weighting = r$weighting, model = model, n_model = nrow(dat),
         class = as.integer(cls), K = nlevels(dat$class), est = cf[, "income_halving"],
         se = sqrt(pmax(diag(V)[terms], 0)),
         n = as.vector(table(dat$class)[cls]), class_mean = as.vector(tapply(dat$cbcl_mean, dat$class, mean)[cls]),
         converged = m$convergence == 0) |>
    mutate(top = class == K, rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

out <- bind_rows(lapply(files, function(fn) {
  r <- readRDS(fn)
  p <- bind_rows(pov_fit(r, "adj"), pov_fit(r, "unadj"))
  stored <- r$poverty$rrr
  got <- p$rrr[p$model == "unadj" & p$top]
  stopifnot(abs(got - stored) < 1e-6)                 # same model as 12_cbma.R
  p
}))
saveRDS(out, file.path(RES, "cbma_poverty.rds"))

top <- out |> filter(top) |> select(scale, weighting, model, rrr, lo, hi) |>
  mutate(across(rrr:hi, \(x) round(x, 2))) |>
  pivot_wider(names_from = model, values_from = c(rrr, lo, hi))
cat("Consensus fits:", length(files), " all unadjusted RRRs match 12_cbma.R.\n")
print(top, n = 50, width = 200)
