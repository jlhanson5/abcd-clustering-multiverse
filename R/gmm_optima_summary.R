# =============================================================================
# R/gmm_optima_summary.R: summarise R/gmm_optima.R.
#   Rscript R/gmm_optima_summary.R
# Per anchor x k:
#   optima    distinct converged optima (log-likelihoods within TOL are one
#             optimum): share of starts, classes, sizes, BIC vs the main run
#   reach     P(reaching the best optimum) with R starts, by resampling starts:
#             gridsearch = converge only the best 10-iteration start (what the
#             main run did, R = 20); converge_all = converge every start
#   poverty   RRR (highest vs lowest class, per halving of income) at the best
#             optimum vs the main run
# Output: output/gmm_optima_summary.rds, output/summary_gmm_optima.txt
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(nnet))

TOL <- 0.5            # log-likelihood units; converged runs this close are one optimum
R_GRID <- c(1, 5, 10, 20, 50, 100)
B <- 5000             # resamples per R
sumf <- file.path(OUT, "summary_gmm_optima.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt  <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

recs <- lapply(list.files(file.path(OUT, "gmm_optima"), "\\.rds$", full.names = TRUE), readRDS)
st <- tibble(anchor = map_chr(recs, "anchor"), k = map_int(recs, "k"), start = map_int(recs, "start"),
             ll_short = map_dbl(recs, "ll_short"), ll_full = map_dbl(recs, "ll_full"),
             conv = map_dbl(recs, "conv"), bic = map_dbl(recs, "bic"), mins = map_dbl(recs, "mins"),
             k_found = map_int(recs, \(r) sum(r$props > 0)),
             min_prop = map_dbl(recs, \(r) min(r$props[r$props > 0])),
             props = map_chr(recs, \(r) paste(sprintf("%.3f", sort(r$props[r$props > 0])), collapse = " ")))

# Main-run anchor fits (BIC, poverty) for comparison.
res <- readRDS(file.path(OUT, "results_full.rds"))
main <- res$fits |> filter(method == "gmm", family == "one", config == "ge3_imp", time == "age",
                           scale %in% c("raw", "log")) |>
  transmute(anchor = scale, k, bic_main = bic, k_found_main = k_found, rrr_main = rrr, lo_main = lo, hi_main = hi)
rm(res); invisible(gc())

# ---- Distinct optima --------------------------------------------------------------
# Only converged runs (conv == 1) count as optima. Within a cell, sort by
# log-likelihood and start a new optimum wherever the gap exceeds TOL.
opt <- st |> filter(conv == 1) |> group_by(anchor, k) |> arrange(desc(ll_full), .by_group = TRUE) |>
  mutate(optimum = cumsum(c(TRUE, -diff(ll_full) > TOL))) |> ungroup()
best <- opt |> filter(optimum == 1) |> group_by(anchor, k) |> summarise(best_ll = max(ll_full), .groups = "drop")
st <- st |> left_join(best, by = c("anchor", "k")) |>
  mutate(at_best = conv == 1 & ll_full > best_ll - TOL)

say(sprintf("GMM optima experiment (%s). %d starts in %d cells.", format(Sys.time(), "%Y-%m-%d %H:%M"),
            nrow(st), n_distinct(paste(st$anchor, st$k))))
say("\nStarts per cell, convergence, time per start:")
say(pt(st |> group_by(anchor, k) |>
  summarise(starts = n(), converged = sum(conv == 1), p_best = round(mean(at_best), 3),
            n_best = sum(at_best), mins_median = round(median(mins), 1), .groups = "drop")))

optima <- opt |> group_by(anchor, k, optimum) |>
  summarise(ll = round(max(ll_full), 1), bic = round(min(bic), 1), n_starts = n(),
            k_found = first(k_found), props = first(props), .groups = "drop") |>
  left_join(count(st, anchor, k, name = "total"), by = c("anchor", "k")) |>
  mutate(share = round(n_starts / total, 3)) |>
  left_join(select(main, anchor, k, bic_main, k_found_main), by = c("anchor", "k")) |>
  mutate(bic_vs_main = round(bic - bic_main, 1))
say("\nDistinct optima (converged runs; bic_vs_main < 0 = better than the main run):")
say(pt(optima |> select(anchor, k, optimum, ll, bic_vs_main, share, k_found, props)))

# ---- Starts needed ---------------------------------------------------------------------
set.seed(SEED)
reach <- st |> group_by(anchor, k) |> group_modify(function(g, key) {
  n <- nrow(g)
  map_dfr(R_GRID[R_GRID <= n], function(R) {
    hit_g <- hit_a <- logical(B)
    for (b in seq_len(B)) {
      s <- sample(n, R)
      hit_g[b] <- g$at_best[s[which.max(g$ll_short[s])]]
      hit_a[b] <- any(g$at_best[s])
    }
    tibble(R = R, gridsearch = mean(hit_g), converge_all = mean(hit_a))
  })
}) |> ungroup()
say("\nP(reaching the best optimum found) with R random starts:")
say(pt(reach |> pivot_longer(c(gridsearch, converge_all), names_to = "strategy") |>
  mutate(value = round(value, 2)) |>
  pivot_wider(names_from = R, values_from = value, names_prefix = "R")))
say("(gridsearch = converge only the best 10-iteration start, as in the main run with R = 20.)")

# ---- Poverty at the best optimum ----------------------------------------------------------
base <- readRDS(file.path(OUT, "stability", "anchor.rds"))
wide <- readRDS(file.path(DER, "wide.rds"))
cov  <- wide[match(base$ids, wide$src_subject_id), ] |>
  transmute(income_halving = -log2(income), age_bl, sex = factor(sex), site = factor(site))
rawmean <- rowMeans(base$R)
pov <- function(cls) {
  m <- tapply(rawmean, cls, mean); cl <- match(cls, as.integer(names(m))[order(m)])
  dat <- cov; dat$class <- factor(cl); K <- nlevels(dat$class)
  fit <- multinom(class ~ income_halving + age_bl + sex + site, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000)
  V <- vcov(fit)
  term <- if (K == 2) "income_halving" else paste0(K, ":income_halving")
  est <- if (K == 2) coef(fit)[["income_halving"]] else coef(fit)[as.character(K), "income_halving"]
  se <- sqrt(V[term, term])
  tibble(K = K, top_prop = mean(cl == K), top_mean = mean(rawmean[cl == K]),
         rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}
key <- st |> filter(at_best) |> group_by(anchor, k) |> slice_max(ll_full, n = 1, with_ties = FALSE) |> ungroup()
poverty <- bind_rows(lapply(seq_len(nrow(key)), function(i) {
  r <- recs[[which(map_chr(recs, "anchor") == key$anchor[i] & map_int(recs, "k") == key$k[i] &
                   map_int(recs, "start") == key$start[i])]]
  bind_cols(key[i, c("anchor", "k")], pov(r$cls))
})) |> left_join(select(main, anchor, k, k_found_main, rrr_main, lo_main, hi_main), by = c("anchor", "k"))
say("\nPoverty RRR at the best optimum vs the main run:")
say(pt(poverty |> mutate(across(where(is.double), \(x) round(x, 3)))))

saveRDS(list(starts = st, optima = optima, reach = reach, poverty = poverty, main = main),
        file.path(OUT, "gmm_optima_summary.rds"))
