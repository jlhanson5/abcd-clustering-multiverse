# =============================================================================
# R/cbma_alt_consensus.R (2026-11-21): does the consensus depend on clusterBMA's
# factorisation step? For each finished clusterBMA fit, the same inputs and the same
# weights are combined two other ways, at the same number of classes K:
#   ea      evidence accumulation (Fred & Jain, 2005): average-linkage hierarchical
#           clustering of 1 - C, where C is the weighted co-assignment (consensus) matrix
#           that clusterBMA factorises; cut at K
#   median  weighted median partition (clue::cl_consensus, method "HE": the hard
#           partition with K classes minimising the weighted sum of squared Euclidean
#           distances to the inputs; Hornik, 2005)
# Each is relabelled 1..K by mean raw CBCL, compared with the clusterBMA labels (ARI),
# and given the race-adjusted poverty model, every class vs the lowest (12b's model).
#   Rscript R/cbma_alt_consensus.R [set] [n_sub]    set: rule (default) | four
# Output: output/results/cbma_alt_consensus[_four].rds and
# output/summaries/summary_cbma_alt_consensus[_four].txt. Aggregates only.
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(mclust); library(clue) })
source(file.path("..", "clusterBMA_port", "R", "cbma.R"))

args <- commandArgs(trailingOnly = TRUE)
SET <- if (length(args) >= 1) args[1] else "rule"
n_sub <- if (length(args) >= 2) as.integer(args[2]) else NA
stopifnot(SET %in% c("rule", "four"))
src <- file.path(OUT, if (SET == "rule") "cbma" else "cbma_four")
sfx <- if (SET == "rule") "" else "_four"

res <- readRDS(file.path(RES, "results_full_refit.rds"))
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()

pov_fit <- function(ids, labels, f = POV_F$adj) {     # as in R/cbma_poverty_ref2.R, reference = class 1
  w <- wide[match(ids, wide$src_subject_id), ]
  dat <- tibble(class = factor(labels), income_halving = -log2(w$income), age_bl = w$age_bl,
                sex = factor(w$sex), site = droplevels(factor(w$site)), race = w$race,
                cbcl_mean = rowMeans(as.matrix(w[paste0("imp_r_", WAVES)])))
  dat <- pov_rows(dat, f); dat$race <- droplevels(dat$race)
  m <- multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE)
  cf <- coef(m); V <- vcov(m)
  if (is.null(dim(cf))) cf <- matrix(cf, 1, dimnames = list(levels(dat$class)[2], names(cf)))
  cls <- rownames(cf)
  terms <- if (length(cls) == 1) "income_halving" else paste0(cls, ":income_halving")
  se <- sqrt(pmax(diag(V)[terms], 0))
  tibble(class = as.integer(cls), pct = as.vector(100 * table(dat$class)[cls] / nrow(dat)),
         class_mean = as.vector(tapply(dat$cbcl_mean, dat$class, mean)[cls]),
         rrr = exp(cf[, "income_halving"]), lo = exp(cf[, "income_halving"] - 1.96 * se),
         hi = exp(cf[, "income_halving"] + 1.96 * se))
}
relabel <- function(cl, raw) match(cl, order(tapply(rowMeans(raw), cl, mean)))

fits <- sub("\\.rds$", "", list.files(src, "\\.rds$"))
out <- list()
for (nm in fits) {
  t0 <- Sys.time()
  r <- readRDS(file.path(src, paste0(nm, ".rds")))
  ids <- r$ids; lab_bma <- r$labels
  if (!is.na(n_sub)) { set.seed(SEED); keep <- sort(sample(length(ids), n_sub)); ids <- ids[keep]; lab_bma <- lab_bma[keep] }
  K <- max(lab_bma)
  raw <- as.matrix(wide[match(ids, wide$src_subject_id), paste0("imp_r_", WAVES)])
  inp <- r$inputs
  cls <- lapply(inp$universe_id, function(u) { a <- res$assign[[u]]; v <- a$class[match(ids, a$src_subject_id)]; stopifnot(!anyNA(v)); v })
  w <- inp$weight

  # Evidence accumulation on the weighted consensus matrix
  C <- cbma_consensus(lapply(cls, cbma_hard_to_prob), w)
  hc <- hclust(as.dist(1 - C), method = "average")
  rm(C); invisible(gc())
  ea <- relabel(cutree(hc, K), raw)
  # How many average-linkage classes hold at least 1% of children at this cut, and the
  # smallest cut that gives K classes of at least 1% (average linkage peels off outliers)
  big <- function(cl) sum(tabulate(cl) >= 0.01 * length(cl))
  k_big <- NA_integer_
  for (kk in K:40) if (big(cutree(hc, kk)) >= K) { k_big <- kk; break }
  ea_big <- if (is.na(k_big)) NULL else {
    cl <- cutree(hc, k_big); keepc <- which(tabulate(cl) >= 0.01 * length(cl))
    v <- ifelse(cl %in% keepc, cl, NA_integer_); v
  }
  rm(hc); invisible(gc())

  # Weighted median partition
  ens <- cl_ensemble(list = lapply(cls, as.cl_hard_partition))
  set.seed(1)
  mp <- cl_consensus(ens, method = "HE", weights = w, control = list(k = K, nruns = 5))
  med <- relabel(as.vector(cl_class_ids(mp)), raw)

  ari <- function(a, b) { ok <- !is.na(a) & !is.na(b); adjustedRandIndex(a[ok], b[ok]) }
  out[[nm]] <- list(
    fit = nm, K = K, n = length(ids), minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
    sizes = list(bma = tabulate(lab_bma, K), ea = tabulate(ea, K), median = tabulate(med, K)),
    ari = c(bma_ea = ari(lab_bma, ea), bma_median = ari(lab_bma, med), ea_median = ari(ea, med),
            bma_ea_big = if (is.null(ea_big)) NA else ari(lab_bma, ea_big)),
    ea_big = list(cut = k_big, n_assigned = if (is.null(ea_big)) NA else sum(!is.na(ea_big))),
    poverty = bind_rows(bma = pov_fit(ids, lab_bma), ea = pov_fit(ids, ea), median = pov_fit(ids, med), .id = "consensus"))
  with(out[[nm]], cat(sprintf("%s (K = %d, %.1f min): ARI bma-ea %.2f, bma-median %.2f, ea-median %.2f; sizes bma %s, ea %s, median %s; EA with >= 1%% classes: cut %s, ARI %.2f\n",
    fit, K, minutes, ari[["bma_ea"]], ari[["bma_median"]], ari[["ea_median"]],
    paste(sizes$bma, collapse = "/"), paste(sizes$ea, collapse = "/"), paste(sizes$median, collapse = "/"),
    ea_big$cut, ari[["bma_ea_big"]])))
}
if (!is.na(n_sub)) quit(save = "no")

saveRDS(out, file.path(RES, paste0("cbma_alt_consensus", sfx, ".rds")))
over <- bind_rows(lapply(out, \(o) tibble(fit = o$fit, K = o$K,
  sizes_bma = paste(round(100 * o$sizes$bma / o$n), collapse = "/"),
  sizes_ea = paste(round(100 * o$sizes$ea / o$n), collapse = "/"),
  sizes_median = paste(round(100 * o$sizes$median / o$n), collapse = "/"),
  ari_bma_ea = round(o$ari[["bma_ea"]], 2), ari_bma_median = round(o$ari[["bma_median"]], 2),
  ari_ea_median = round(o$ari[["ea_median"]], 2), ea_cut_1pct = o$ea_big$cut,
  ari_bma_ea_1pct = round(o$ari[["bma_ea_big"]], 2))))
pv <- bind_rows(lapply(out, \(o) mutate(o$poverty, fit = o$fit))) |>
  mutate(class_mean = round(class_mean, 1), pct = round(pct, 1), ci = sprintf("%.2f-%.2f", lo, hi), rrr = round(rrr, 2)) |>
  select(fit, consensus, class, pct, class_mean, rrr, ci)
txt <- c(sprintf("Alternative consensus methods on the clusterBMA %s fits (%s). Same inputs and weights; K = the clusterBMA fit's number of classes.",
                 if (SET == "rule") "primary (k-means + kml, k 2-8)" else "four-method sensitivity", format(Sys.time(), "%Y-%m-%d %H:%M")),
         "ea: evidence accumulation (average linkage on 1 - weighted co-assignment), cut at K. median: weighted median partition (clue, HE).",
         "sizes: % of children per class, lowest to highest mean raw CBCL. ea_cut_1pct: the smallest cut with K classes of >= 1% (smaller classes set aside); ari_bma_ea_1pct on the children in those classes.", "",
         capture.output(print(over, n = 50, width = 220)), "",
         "Poverty RRR per halving of income, every class vs class 1 (the lowest), adjusted for age, sex, site and race/ethnicity:",
         capture.output(print(pv, n = 200, width = 200)))
writeLines(txt, file.path(SUMM, paste0("summary_cbma_alt_consensus", sfx, ".txt")))
cat("written\n")
