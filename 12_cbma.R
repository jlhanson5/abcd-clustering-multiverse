# =============================================================================
# 12_cbma.R: clusterBMA consensus across methods and k, one scale and one weighting.
#   Rscript 12_cbma.R <scale> <weighting> [n_sub] [n_iter]
#     scale      raw | t | log | all (all three scales pooled; weighting equal or kbal)
#     weighting  ch (Calinski-Harabasz, the paper's default) | sdbw (S_Dbw, the
#                paper's alternative) | equal (no index) | avg (mean of the ch, sdbw
#                and equal weights; needs those three runs first) | kbal (each k
#                gets 1/5 of the weight, split equally among that k's inputs)
#     n_sub, n_iter  smoke tests only: a random subsample of children, fewer steps
# Decided 2026-11-07. Uses the native R port, ../clusterBMA_port/R/cbma.R
# (Forbes et al., 2023); no ABCD data goes into that repository.
#
# Inputs: the primary run's (full_refit) solutions for the anchor configuration
# (one child per family, >= 3 waves, imputed; n = 8,331): all 7 methods at k = 2-6,
# GBTM/GMM on the age metric. Left out: fits with an error or an empty class, GMM
# k = 5-6 and raw GMM k = 4 (not reliably estimated; README "GMM refit"). Inputs are
# hard (one-hot): only modal classes were saved. K_BMA = 6, the largest input k.
# The index is computed on the scale's wave matrix (raw, T-score or log1p(raw)),
# because index values do not compare across scales.
# avg (added 2026-11-07): the consensus matrix is linear in the weights, so the mean
# of the three runs' consensus matrices equals the consensus from their mean weights;
# this run averages over the three weighting schemes with equal prior weight.
# all: equal weights need no index, so the three scales' inputs can be pooled into
# one consensus across the whole anchor multiverse.
# kbal (added 2026-11-08): the leave-outs make k uneven (7 inputs at k = 2-3, 5-6 at
# k = 5-6; the missing ones are the model-based methods), so equal weighting gives
# low k more votes. kbal gives every k the same total weight, as a sensitivity check.
#
# Output: output/cbma/<scale>_<weighting>.rds with inputs, index values, weights,
# allocation matrix, labels (relabelled 1..K by mean raw CBCL), uncertainty, empty
# clusters, objective trace, ARI of each input with the consensus, and the poverty
# model on the consensus classes (04_poverty.R's model; one per family, so model SEs).
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(mclust) })
PORT <- normalizePath(file.path("..", "clusterBMA_port"))
source(file.path(PORT, "R", "cbma.R"))

args <- commandArgs(trailingOnly = TRUE)
scale <- args[1]; weighting <- args[2]
n_sub <- if (length(args) >= 3) as.integer(args[3]) else NA
n_iter <- if (length(args) >= 4) as.integer(args[4]) else 5000
stopifnot(scale %in% c("raw", "t", "log", "all"), weighting %in% c("ch", "sdbw", "equal", "avg", "kbal"),
          scale != "all" || weighting %in% c("equal", "kbal"))
scales <- if (scale == "all") c("raw", "t", "log") else scale
smoke <- !is.na(n_sub)
dest <- file.path(OUT, if (smoke) "cbma_smoke" else "cbma")
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
K_BMA <- 6

res <- readRDS(file.path(RES, "results_full_refit.rds"))
wide <- readRDS(file.path(DER, "wide.rds"))
anchor <- res$fits |>
  filter(family == "one", config == "ge3_imp", scale %in% scales, k %in% 2:6,
         time == ifelse(method %in% c("gbtm", "gmm"), "age", "wave")) |>
  mutate(keep = is.na(error) & !degenerate &
           !(method == "gmm" & k >= 5) & !(method == "gmm" & scale == "raw" & k == 4))
inputs <- filter(anchor, keep) |>
  arrange(match(scale, c("raw", "t", "log")), match(method, c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm")), k)

ids <- sort(res$assign[[inputs$universe_id[1]]]$src_subject_id)
if (smoke) { set.seed(SEED); ids <- sort(sample(ids, n_sub)) }
w <- wide[match(ids, wide$src_subject_id), ]
raw <- as.matrix(w[paste0("imp_r_", WAVES)])
X <- switch(scale, raw = raw, t = as.matrix(w[paste0("imp_t_", WAVES)]), log = log1p(raw), all = NULL)
A_list <- lapply(inputs$universe_id, function(u) {
  a <- res$assign[[u]]
  cls <- a$class[match(ids, a$src_subject_id)]
  stopifnot(!anyNA(cls))
  cbma_hard_to_prob(cls)
})
names(A_list) <- inputs$universe_id
cat(sprintf("[%s / %s] %d inputs (%d left out), N = %d, K_BMA = %d, %d steps\n", scale, weighting,
            length(A_list), sum(!anchor$keep), length(ids), K_BMA, n_iter))

t0 <- Sys.time()
index <- switch(weighting, ch = "Calinski_Harabasz", sdbw = "S_Dbw", "Calinski_Harabasz")
if (scale == "all") {
  values <- rep(NA_real_, length(A_list))           # no index across scales
} else {
  civi <- lapply(A_list, \(A) cbma_civi(X, A, index))
  values <- vapply(civi, `[[`, 0, "value")
}
weights <- switch(weighting,
  equal = rep(1 / length(A_list), length(A_list)),
  kbal = { nk <- table(inputs$k); as.vector(1 / (length(nk) * nk[as.character(inputs$k)])) },
  avg = {
    w3 <- sapply(c("ch", "sdbw", "equal"), function(wt) {
      r <- readRDS(file.path(dest, sprintf("%s_%s.rds", scale, wt)))
      r$inputs$weight[match(names(A_list), r$inputs$universe_id)]
    })
    stopifnot(!anyNA(w3))
    rowMeans(w3)
  },
  cbma_weights(values, civi[[1]]$direction))
C <- cbma_consensus(A_list, weights)
t_index <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
f <- cbma_ssmf(C, K_BMA, n_iter = n_iter, seed = 1L)
rm(C); invisible(gc())
t_ssmf <- as.numeric(difftime(Sys.time(), t0, units = "mins")) - t_index

# Relabel 1..K by mean raw CBCL (04_poverty.R's rule); empty clusters drop out.
lab0 <- max.col(f$alloc, ties.method = "first")
ord <- order(tapply(rowMeans(raw), factor(lab0, seq_len(K_BMA)), mean), na.last = NA)
labels <- match(lab0, ord)
alloc <- f$alloc[, ord, drop = FALSE]
uncert <- cbma_uncertainty(f$alloc)

# Poverty model on the consensus classes, as in 04_poverty.R.
dat <- w |> transmute(class = factor(labels), income_halving = -log2(income), age_bl,
                      sex = factor(sex), site = droplevels(factor(site)))
K <- nlevels(dat$class)
pov <- tryCatch({
  m <- multinom(class ~ income_halving + age_bl + sex + site, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000)
  V <- vcov(m); cf <- coef(m)
  term <- if (K == 2) "income_halving" else paste0(K, ":income_halving")
  est <- if (K == 2) cf[["income_halving"]] else cf[as.character(K), "income_halving"]
  se <- sqrt(V[term, term])
  tibble(K = K, est = est, se = se, rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
         top_n = sum(labels == K), top_mean = mean(rowMeans(raw)[labels == K]))
}, error = function(e) tibble(K = K, error = conditionMessage(e)))

agree <- inputs |> select(universe_id, input_scale = scale, method, k) |>
  mutate(index_value = values, weight = weights,
         ari = vapply(A_list, \(A) adjustedRandIndex(max.col(A, ties.method = "first"), labels), 0))
port_git <- tryCatch(system2("git", c("-C", PORT, "describe", "--always", "--dirty"), stdout = TRUE),
                     error = function(e) NA_character_)

out <- list(scale = scale, weighting = weighting, index = if (weighting %in% c("equal", "avg", "kbal")) NA else index,
            n = length(ids), ids = ids, inputs = agree, left_out = filter(anchor, !keep) |> select(universe_id, method, k),
            alloc = alloc, labels = labels, uncertainty = uncert,
            sizes = tabulate(labels, K), empty_clusters = length(f$empty_clusters),
            objective_trace = f$objective_trace, poverty = pov,
            settings = list(K_BMA = K_BMA, n_iter = n_iter, seed = 1L, lambda = 1, lr = 0.01),
            port = port_git, blas = basename(extSoftVersion()[["BLAS"]]),
            minutes = c(index_consensus = t_index, ssmf = t_ssmf), date = Sys.time())
fn <- file.path(dest, sprintf("%s_%s.rds", scale, weighting))
saveRDS(out, paste0(fn, ".tmp")); invisible(file.rename(paste0(fn, ".tmp"), fn))
cat(sprintf("Done in %.1f min (index + consensus %.1f, factorisation %.1f). %d of %d clusters used; sizes %s.\n",
            t_index + t_ssmf, t_index, t_ssmf, K, K_BMA, paste(out$sizes, collapse = "/")))
cat(sprintf("Weights: max %.3f, min %.3f. Poverty RRR (highest vs lowest): %s\n", max(weights), min(weights),
            if (!is.null(pov$rrr)) sprintf("%.2f (%.2f-%.2f)", pov$rrr, pov$lo, pov$hi) else pov$error))
