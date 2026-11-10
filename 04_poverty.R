# =============================================================================
# 04_poverty.R: collect fitted universes and run the poverty model in each.
#   Rscript 04_poverty.R pilot|full
#
# Classes are relabelled 1..K by mean raw CBCL across the 4 waves (imputed raw,
# whatever scale the universe clustered on), so class 1 = lowest and class K =
# highest. Empty classes are dropped first, so K = number of non-empty classes.
#
# Poverty model, every universe: multinomial logit (nnet::multinom)
#   primary (since 2026-11-09):  class ~ income_halving + age_bl + sex + site + race
#   sensitivity ("_unadj"):      class ~ income_halving + age_bl + sex + site
# reference = lowest class; formulas in 00_setup.R (POV_F). The one child with no
# race code is left out of the adjusted model only.
# income_halving = -log2(income), so exp(b) is the relative risk ratio (RRR)
# for one halving of income. The main spec curve uses highest vs lowest.
# SEs: model-based for one-per-family universes; for full-sample universes a
# family-clustered sandwich (bread = inverse Hessian from multinom).
#
# Output: output/results_<mode>.rds with
#   fits      : one row per universe (spec + fit indices + main contrast, adjusted;
#               est/se/rrr/lo/hi/p_unadj from the unadjusted model)
#   contrasts, contrasts_unadj : one row per non-reference class per universe
#   assign    : list of class assignments by universe_id (relabelled)
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(parallel) })

mode  <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "pilot"
cores <- 8

wide <- readRDS(file.path(DER, "wide.rds")) |> add_race() |>
  mutate(income_halving = -log2(income),
         cbcl_mean = rowMeans(across(all_of(paste0("imp_r_", WAVES)))),
         sex = factor(sex), site = factor(site))
univ <- readRDS(file.path(RES, "spec_universes.rds"))

files <- list.files(file.path(UNIV, mode), full.names = TRUE)
recs  <- unlist(lapply(files, \(f) readRDS(f)$fits), recursive = FALSE)
fits  <- bind_rows(lapply(recs, `[[`, "fit")) |>
  mutate(universe_id = paste0(task_id, "_k", k))
assign <- setNames(lapply(recs, `[[`, "assign"), fits$universe_id)
assign <- assign[!vapply(assign, is.null, logical(1))]

# ---- Relabel by class mean CBCL ------------------------------------------------
relabel <- function(a) {
  a <- left_join(a, select(wide, src_subject_id, cbcl_mean), by = "src_subject_id")
  ord <- a |> group_by(class) |> summarise(m = mean(cbcl_mean)) |> arrange(m)
  a |> mutate(class = match(class, ord$class)) |> select(src_subject_id, class)
}
assign <- lapply(assign, relabel)

# ---- Family-clustered sandwich for multinom ----------------------------------------
vcov_cluster <- function(fit, X, y, cluster) {
  B <- vcov(fit)
  lev <- fit$lev
  P <- as.matrix(fitted(fit)); if (ncol(P) == 1) P <- cbind(1 - P, P)
  colnames(P) <- lev
  S <- do.call(cbind, lapply(lev[-1], \(j) ((y == j) - P[, j]) * X))
  G <- length(unique(cluster))
  M <- crossprod(rowsum(S, cluster)) * G / (G - 1)
  B %*% M %*% B
}

# ---- Poverty model, one universe -------------------------------------------------
poverty_one <- function(uid, f) {
  a <- assign[[uid]]
  dat <- pov_rows(inner_join(a, wide, by = "src_subject_id"), f) |>
    mutate(class = factor(class), site = droplevels(site), race = droplevels(race))
  K <- nlevels(dat$class)
  cls_mean <- dat |> group_by(class) |> summarise(class_mean = mean(cbcl_mean), n = n(), .groups = "drop")
  if (K < 2) return(NULL)
  fit <- tryCatch(multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  full <- grepl("_full_", uid)
  V <- if (full) {
    X <- model.matrix(f, dat)
    tryCatch(vcov_cluster(fit, X, as.character(dat$class), dat$family), error = function(e) NULL)
  } else tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(V)) return(NULL)
  cf <- coef(fit); if (is.null(dim(cf))) cf <- matrix(cf, nrow = 1, dimnames = list(levels(dat$class)[2], names(cf)))
  terms <- paste0(rownames(cf), ":income_halving")
  if (K == 2) terms <- "income_halving"
  est <- cf[, "income_halving"]
  se  <- sqrt(pmax(diag(V)[terms], 0))
  tibble(universe_id = uid, class = as.integer(rownames(cf)), K = K,
         est = unname(est), se = unname(se)) |>
    left_join(mutate(cls_mean, class = as.integer(as.character(class))), by = "class") |>
    mutate(converged = fit$convergence == 0)
}

run_model <- function(f) {
  t0 <- Sys.time()
  out <- bind_rows(mclapply(names(assign), poverty_one, f = f, mc.cores = cores)) |>
    mutate(rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
           p = 2 * pnorm(-abs(est / se)))
  cat(sprintf("Poverty models (%s): %d universes, %.1f min\n", deparse(f[[3]]),
              n_distinct(out$universe_id), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out
}
contrasts <- run_model(POV_F$adj)
contrasts_unadj <- run_model(POV_F$unadj)

main <- contrasts |> filter(class == K) |>
  select(universe_id, K, est, se, rrr, lo, hi, p, top_mean = class_mean, top_n = n,
         mn_converged = converged)
main_unadj <- contrasts_unadj |> filter(class == K) |>
  select(universe_id, est, se, rrr, lo, hi, p, mn_converged = converged) |>
  rename_with(\(x) paste0(x, "_unadj"), -universe_id)
fits <- univ |>
  inner_join(fits, by = c("universe_id", "task_id", "k")) |>
  left_join(main, by = "universe_id") |>
  left_join(main_unadj, by = "universe_id") |>
  mutate(degenerate = !is.na(k_found) & k_found < k,
         top_too_small = !is.na(top_n) & top_n < MIN_TOP_N)
# Highest class under MIN_TOP_N children: main contrast not estimable (00_setup.R).
# The per-class contrasts tables keep every fitted value.
pov_cols <- c("est", "se", "rrr", "lo", "hi", "p")
fits[fits$top_too_small, c(pov_cols, paste0(pov_cols, "_unadj"))] <- NA

saveRDS(list(fits = fits, contrasts = contrasts, contrasts_unadj = contrasts_unadj, assign = assign,
             covariates = list(primary = deparse(POV_F$adj), sensitivity = deparse(POV_F$unadj),
                               min_top_n = MIN_TOP_N)),
        file.path(RES, paste0("results_", mode, ".rds")))
cat("Not estimable (highest class <", MIN_TOP_N, "children):", sum(fits$top_too_small), "\n")
cat("Universes:", nrow(fits), " errors:", sum(!is.na(fits$error)),
    " degenerate (empty class):", sum(fits$degenerate, na.rm = TRUE),
    " small class <5%:", sum(fits$small_class, na.rm = TRUE),
    " with poverty estimate:", sum(!is.na(fits$est)), "\n")
