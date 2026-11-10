# =============================================================================
# 18_sens_weights.R: sensitivity analysis, ABCD population weights (added 2026-11-09).
#   ABCD_REL=<folder with abcd-general/abcd_p_demo.csv> Rscript 18_sens_weights.R [cores] [n_test]
# Weight: acs_raked_propensity_score (abcd_p_demo, baseline; called acspsw03 up to
# release 4.0): ABCD's propensity weight to the American Community Survey, raked to ACS
# totals by age, sex and race/ethnicity. Rescaled to mean 1 within each model's sample.
# Model: the primary poverty model (race-adjusted multinomial, highest vs lowest, classes
# as in 04_poverty.R) fitted with these weights. SEs: sandwich (bread = inverse weighted
# Hessian, meat = weighted scores), clustered by family in full-sample universes and by
# child otherwise. Children without a weight are left out (count in the summary).
# The unweighted model is refitted on the same children with the same sandwich SEs, so the
# effect of weighting is separate from the change of SE type. Check: in full-sample
# universes this unweighted fit must reproduce the primary estimate and SE.
# Primary run (full_refit), main universes (k 2-6) with an estimable RRR.
# Output: data_derived/weights.rds; output/results/sens_weights.rds;
#         output/summaries/summary_sens_weights.txt; figures/7_sensitivity/s3_population_weights.png
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(parallel) })

args   <- commandArgs(trailingOnly = TRUE)
cores  <- as.integer(args[1]); if (is.na(cores)) cores <- 1
n_test <- as.integer(args[2])
figd <- fig_dir("7_sensitivity")
sumf <- file.path(SUMM, "summary_sens_weights.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

# ---- Weights ------------------------------------------------------------------------------
wfile <- file.path(DER, "weights.rds")
if (!file.exists(wfile)) {
  read_abcd("abcd-general/abcd_p_demo.csv", "acs_raked_propensity_score") |>
    filter(eventname == EV_BL) |>
    transmute(src_subject_id, weight = num(acs_raked_propensity_score)) |>
    filter(!is.na(weight)) |> saveRDS(wfile)
}
wts <- readRDS(wfile)

res  <- readRDS(file.path(RES, "results_full_refit.rds"))
f    <- POV_F$adj
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race() |>
  mutate(income_halving = -log2(income), sex = factor(sex), site = factor(site)) |>
  left_join(wts, by = "src_subject_id")
univ <- res$fits |> filter(!supplement, !is.na(est)) |>
  select(universe_id, method, scale, k, family_fork = family, K, est, se, rrr, lo, hi)
if (!is.na(n_test)) { set.seed(1); univ <- univ[sort(sample(nrow(univ), n_test)), ] }

# Sandwich for a (weighted) multinom: bread = vcov(fit) (inverse weighted Hessian),
# meat from weighted scores summed within clusters.
vcov_sandwich <- function(fit, X, y, w, cluster) {
  B <- vcov(fit); lev <- fit$lev
  P <- as.matrix(fitted(fit)); if (ncol(P) == 1) P <- cbind(1 - P, P)
  colnames(P) <- lev
  S <- do.call(cbind, lapply(lev[-1], \(j) w * ((y == j) - P[, j]) * X))
  G <- length(unique(cluster))
  B %*% (crossprod(rowsum(S, cluster)) * G / (G - 1)) %*% B
}
top_fit <- function(dat, K, w, cluster) {
  dat$w <- w
  fit <- tryCatch(multinom(f, data = dat, weights = w, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  V <- tryCatch(vcov_sandwich(fit, model.matrix(f, dat), as.character(dat$class), w, cluster), error = function(e) NULL)
  if (is.null(V)) return(c(NA, NA))
  cf <- coef(fit)
  term <- if (K == 2) "income_halving" else paste0(K, ":income_halving")
  est <- if (K == 2) cf[["income_halving"]] else cf[as.character(K), "income_halving"]
  c(est, sqrt(max(V[term, term], 0)))
}
fit_one <- function(i) {
  u <- univ[i, ]
  dat <- pov_rows(inner_join(res$assign[[u$universe_id]], wide, by = "src_subject_id"), f)
  n_all <- nrow(dat)
  dat <- filter(dat, !is.na(weight)) |>
    mutate(class = factor(class), site = droplevels(site), race = droplevels(race))
  top_n <- sum(dat$class == u$K)
  if (nlevels(dat$class) != u$K || top_n < MIN_TOP_N)
    return(tibble(universe_id = u$universe_id, n_model = nrow(dat), n_no_weight = n_all - nrow(dat), top_n = top_n,
                  est_unw = NA_real_, se_unw = NA_real_, est_w = NA_real_, se_w = NA_real_))
  cl <- if (u$family_fork == "full") dat$family else dat$src_subject_id
  unw <- top_fit(dat, u$K, rep(1, nrow(dat)), cl)
  wtd <- top_fit(dat, u$K, dat$weight / mean(dat$weight), cl)
  tibble(universe_id = u$universe_id, n_model = nrow(dat), n_no_weight = n_all - nrow(dat), top_n = top_n,
         est_unw = unw[1], se_unw = unw[2], est_w = wtd[1], se_w = wtd[2])
}

t0  <- Sys.time()
out <- bind_rows(mclapply(seq_len(nrow(univ)), fit_one, mc.cores = cores))
sep <- \(se) !is.finite(se) | se <= 0 | se > MAX_SE   # separation: not estimable (00_setup.R)
s <- univ |> inner_join(out, by = "universe_id") |>
  mutate(across(c(est_w, se_w), \(x) ifelse(sep(se_w), NA, x)),
         across(c(est_unw, se_unw), \(x) ifelse(sep(se_unw), NA, x)),
         rrr_unw = exp(est_unw), lo_unw = exp(est_unw - 1.96 * se_unw),
         rrr_w = exp(est_w), lo_w = exp(est_w - 1.96 * se_w), hi_w = exp(est_w + 1.96 * se_w),
         scale_lab = factor(c(raw = "Raw", t = "T-score", log = "log(1+raw)")[scale], c("Raw", "T-score", "log(1+raw)")))
saveRDS(s, file.path(RES, "sens_weights.rds"))

# Check: unweighted, family-clustered, all children with a weight. Equals the primary
# fit only when no child in the universe lacks a weight.
chk <- filter(s, family_fork == "full", n_no_weight == 0, !is.na(est_unw))
if (nrow(chk)) stopifnot(max(abs(chk$est_unw - chk$est)) < 1e-6, max(abs(chk$se_unw - chk$se)) < 1e-6)

say(sprintf("Population weights (%s, %.1f min). %d universes; weight %s; %d children in the main sample, %d without a weight.",
            format(Sys.time(), "%Y-%m-%d %H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(s),
            "acs_raked_propensity_score (baseline)", sum(!is.na(wide$income_halving)), sum(is.na(wide$weight))))
say(sprintf("Not estimable after dropping unweighted children: unweighted %d, weighted %d. Check vs primary: %d full-sample universes reproduced.",
            sum(is.na(s$est_unw)), sum(is.na(s$est_w)), nrow(chk)))
say("RRR columns: primary (stored); unweighted with sandwich SEs; weighted with sandwich SEs. Sandwich clustered by family",
    "(full sample) or child (one per family).")
summ <- function(d) d |> summarise(n = sum(!is.na(rrr_w)), median_primary = round(median(rrr), 3),
                                   median_unw = round(median(rrr_unw, na.rm = TRUE), 3), median_w = round(median(rrr_w, na.rm = TRUE), 3),
                                   above1_primary = round(mean(lo > 1), 3), above1_unw = round(mean(lo_unw > 1, na.rm = TRUE), 3),
                                   above1_w = round(mean(lo_w > 1, na.rm = TRUE), 3),
                                   ratio_w_unw = round(median(rrr_w / rrr_unw, na.rm = TRUE), 3),
                                   se_ratio_w_unw = round(median(se_w / se_unw, na.rm = TRUE), 2), .groups = "drop")
say(pt(bind_rows(summ(s) |> mutate(scale_lab = "All", .before = 1),
                 s |> group_by(scale_lab) |> summ() |> mutate(scale_lab = as.character(scale_lab)))))

S_COL <- c(Raw = "#2a78d6", `T-score` = "#eb6834", `log(1+raw)` = "#1baf7a")
INK2 <- "#52514e"; GRID <- "#ecebe8"
p <- ggplot(filter(s, !is.na(rrr_w), !is.na(rrr_unw)), aes(rrr_unw, rrr_w, color = scale_lab)) +
  geom_abline(color = INK2, linewidth = 0.4, linetype = "22") +
  geom_point(size = 1.3, alpha = 0.6) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = S_COL, guide = "none") +
  labs(x = "RRR, unweighted (same children)", y = "RRR, weighted to the ACS population",
       title = "ABCD population weights",
       subtitle = sprintf("Every universe (k 2-6, n = %d), highest vs lowest class, adjusted for age, sex, site and race/ethnicity. Dashed line: no change.",
                          sum(!is.na(s$rrr_w) & !is.na(s$rrr_unw)))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(figd, "s3_population_weights.png"), p, width = 10, height = 4.2, dpi = 300)
say("\nFigure written to ", figd)
