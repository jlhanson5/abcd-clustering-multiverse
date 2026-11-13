# =============================================================================
# 17_permutation.R: specification-curve permutation test (Simonsohn, Simmons &
# Nelson 2020, Nature Human Behaviour), added 2026-11-09.
#   Rscript 17_permutation.R <n_perm> [cores]      # restartable; perm 0 = observed
#   Rscript 17_permutation.R summary               # tables + figure from saved perms
# Null: income has no association with class, other links kept. Freedman-Lane:
# income_halving = fitted(age, sex, site, race) + residual; the residuals are permuted
# once per permutation over all 10,088 children, and the same permuted income is used
# in every universe, so the dependence between universes is kept. Classes are not
# refitted (the clustering does not use income).
# Every permutation reruns the primary poverty model (race-adjusted multinomial,
# highest vs lowest, sandwich SEs for full-sample universes as in 04_poverty.R) in the
# main universes with an estimable RRR (652 before the skew-aware GBTM was added on 2026-11-13). Statistics per curve: median RRR, share of
# universes with the CI above 1, Stouffer's Z (sum z / sqrt(n)). p = share of
# permutations at least as extreme as observed, counting the observed: (b + 1) / (B + 1).
# Output: output/permutation/perm_<id>.rds; output/results/permutation.rds;
#         output/summaries/summary_permutation.txt; figures/7_sensitivity/s2_permutation.png
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
PDIR <- file.path(OUT, "permutation"); dir.create(PDIR, showWarnings = FALSE)

res  <- readRDS(file.path(RES, "results_full_refit.rds"))
univ <- res$fits |> filter(!supplement, !is.na(est)) |> select(universe_id, K, est, se, rrr, lo, hi)
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race() |>
  mutate(income_halving = -log2(income), sex = factor(sex), site = factor(site))
f <- POV_F$adj

# ---- Freedman-Lane permuted income --------------------------------------------------------
ok <- !is.na(wide$race)
inc_fit <- lm(income_halving ~ age_bl + sex + site + race, data = wide[ok, ])
inc_hat <- fitted(inc_fit); inc_res <- resid(inc_fit)
perm_income <- function(b) {
  x <- wide$income_halving
  if (b > 0) { set.seed(SEED + b); x[ok] <- inc_hat + sample(inc_res) }
  x
}

# ---- Primary poverty model (as in 04_poverty.R), highest vs lowest only ---------------------
vcov_cluster <- function(fit, X, y, cluster) {
  B <- vcov(fit); lev <- fit$lev
  P <- as.matrix(fitted(fit)); if (ncol(P) == 1) P <- cbind(1 - P, P)
  colnames(P) <- lev
  S <- do.call(cbind, lapply(lev[-1], \(j) ((y == j) - P[, j]) * X))
  G <- length(unique(cluster))
  M <- crossprod(rowsum(S, cluster)) * G / (G - 1)
  B %*% M %*% B
}
top_one <- function(uid, K, w) {
  dat <- pov_rows(inner_join(res$assign[[uid]], w, by = "src_subject_id"), f) |>
    mutate(class = factor(class), site = droplevels(site), race = droplevels(race))
  fit <- tryCatch(multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA))
  V <- if (grepl("_full_", uid)) tryCatch(vcov_cluster(fit, model.matrix(f, dat), as.character(dat$class), dat$family),
                                         error = function(e) NULL) else vcov(fit)
  if (is.null(V)) return(c(NA, NA))
  cf <- coef(fit)
  if (K == 2) return(c(cf[["income_halving"]], sqrt(V["income_halving", "income_halving"])))
  term <- paste0(K, ":income_halving")
  c(cf[as.character(K), "income_halving"], sqrt(V[term, term]))
}
# A saved permutation is topped up with any universes it lacks (e.g. the skew-aware GBTM
# universes added 2026-11-13): the permuted income depends only on b, and each universe
# is fitted on its own, so this equals refitting the permutation from scratch.
run_perm <- function(b, cores) {
  fn <- file.path(PDIR, sprintf("perm_%04d.rds", b))
  old <- if (file.exists(fn)) readRDS(fn) else NULL
  todo <- if (is.null(old)) univ else filter(univ, !universe_id %in% old$universe_id)
  if (!nrow(todo)) return(invisible(NULL))
  w <- wide; w$income_halving <- perm_income(b)
  w <- select(w, src_subject_id, income_halving, age_bl, sex, site, race, family)
  r <- mclapply(seq_len(nrow(todo)), \(i) top_one(todo$universe_id[i], todo$K[i], w), mc.cores = cores)
  m <- do.call(rbind, r)
  out <- bind_rows(old, tibble(universe_id = todo$universe_id, est = m[, 1], se = m[, 2]))
  saveRDS(out, paste0(fn, ".tmp")); file.rename(paste0(fn, ".tmp"), fn)
}
curve_stats <- function(d) {
  d <- filter(d, is.finite(est), is.finite(se), se > 0)
  z <- d$est / d$se
  tibble(n = nrow(d), median_rrr = exp(median(d$est)), share_above = mean(z > 1.96),
         share_below = mean(z < -1.96), stouffer_z = sum(z) / sqrt(nrow(d)))
}

if (length(args) >= 1 && args[1] == "summary") {
  files <- list.files(PDIR, "^perm_\\d{4}\\.rds$", full.names = TRUE)
  ids <- as.integer(sub("perm_(\\d{4})\\.rds", "\\1", basename(files)))
  st <- bind_rows(lapply(seq_along(files), \(i) mutate(curve_stats(filter(readRDS(files[i]), universe_id %in% univ$universe_id)),
                                                       perm = ids[i], .before = 1)))
  obs <- filter(st, perm == 0); null <- filter(st, perm > 0); B <- nrow(null)
  p <- c(median_rrr = (sum(null$median_rrr >= obs$median_rrr) + 1) / (B + 1),
         share_above = (sum(null$share_above >= obs$share_above) + 1) / (B + 1),
         stouffer_z = (sum(null$stouffer_z >= obs$stouffer_z) + 1) / (B + 1))
  saveRDS(list(stats = st, p = p, B = B), file.path(RES, "permutation.rds"))
  sumf <- file.path(SUMM, "summary_permutation.txt"); unlink(sumf)
  say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
  say(sprintf("Specification-curve permutation test (%s): %d permutations, %d universes.",
              format(Sys.time(), "%Y-%m-%d %H:%M"), B, obs$n))
  say(sprintf("Observed: median RRR %.3f, CI above 1 in %.1f%%, below 1 in %.1f%%, Stouffer Z %.1f.",
              obs$median_rrr, 100 * obs$share_above, 100 * obs$share_below, obs$stouffer_z))
  say(sprintf("Null (permutations): median RRR %.3f (95%% range %.3f-%.3f); CI above 1 in %.1f%% (95%% range %.1f-%.1f%%); Stouffer Z %.2f (95%% range %.2f-%.2f).",
              median(null$median_rrr), quantile(null$median_rrr, .025), quantile(null$median_rrr, .975),
              100 * median(null$share_above), 100 * quantile(null$share_above, .025), 100 * quantile(null$share_above, .975),
              median(null$stouffer_z), quantile(null$stouffer_z, .025), quantile(null$stouffer_z, .975)))
  say(sprintf("p (one-sided, (b + 1) / (B + 1)): median RRR %.4f, share above %.4f, Stouffer Z %.4f.", p[1], p[2], p[3]))
  INK2 <- "#52514e"; GRID <- "#ecebe8"
  pl <- function(col, lab, fmt) {
    ggplot(null, aes(.data[[col]])) +
      geom_histogram(bins = 40, fill = "#b8b6b0", color = "white", linewidth = 0.2) +
      geom_vline(xintercept = obs[[col]], color = "#2a78d6", linewidth = 0.8) +
      annotate("text", x = obs[[col]], y = Inf, label = paste("Observed", sprintf(fmt, obs[[col]])),
               hjust = 1.05, vjust = 1.5, size = 3, color = "#2a78d6") +
      labs(x = lab, y = if (col == "median_rrr") "Permutations" else NULL) +
      theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID))
  }
  suppressPackageStartupMessages(library(patchwork))
  fig <- (pl("median_rrr", "Median RRR", "%.2f") | pl("share_above", "Share of universes with CI above 1", "%.2f") |
          pl("stouffer_z", "Stouffer's Z", "%.0f")) +
    plot_annotation(title = "The observed specification curve against curves with no income effect",
                    subtitle = sprintf("Grey: %d permutations (Freedman-Lane, income residuals shuffled once per permutation across all universes). Blue: observed.", B),
                    theme = theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = INK2),
                                  plot.background = element_rect(fill = "white", color = NA)))
  ggsave(file.path(fig_dir("7_sensitivity"), "s2_permutation.png"), fig, width = 11, height = 3.8, dpi = 300)
  quit(save = "no")
}

n_perm <- as.integer(args[1]); stopifnot(!is.na(n_perm))
cores  <- if (length(args) >= 2) as.integer(args[2]) else 8
logf <- file.path(LOGS, "permutation.log")
for (b in 0:n_perm) {
  t0 <- Sys.time(); run_perm(b, cores)
  if (b == 0) {   # perm 0 must reproduce the primary results
    o <- readRDS(file.path(PDIR, "perm_0000.rds")) |> inner_join(univ, by = "universe_id", suffix = c("", "_primary"))
    stopifnot(max(abs(o$est - o$est_primary)) < 1e-8, max(abs(o$se - o$se_primary)) < 1e-8)
  }
  cat(sprintf("%s perm %d done (%.1f min)\n", format(Sys.time(), "%m-%d %H:%M"), b,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = logf, append = TRUE)
}
