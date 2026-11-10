# =============================================================================
# 16_sens_site_re.R: sensitivity analysis, site as a random effect (added 2026-11-09).
#   Rscript 16_sens_site_re.R [cores] [n_test]          # default 8 cores, all universes
# The primary model has site as fixed effects (21 dummies) in a multinomial logit.
# Here the highest-vs-lowest contrast is fitted as a binary logit on the children in
# those two classes (Begg & Gray 1984: consistent for the multinomial contrast), with a
# random intercept for site (lme4::glmer, Laplace). Same covariates as the primary model
# (age, sex, race/ethnicity). The same binary model with site as fixed effects separates
# the binary-vs-multinomial step from the random-vs-fixed step.
# Site only, in every universe. A first run (2026-11-09) also gave full-sample universes
# a random intercept for family: most families have one child, the term was barely
# identified (24 fits failed, raw-scale SEs 16 times the fixed-effects SEs). That run is
# kept as output/results/sens_site_re_famRE_attempt.rds (+ _famRE_attempt summary, figure).
# Full-sample universes therefore ignore sibling clustering here (the primary model uses
# family-clustered sandwich SEs), so their CIs are slightly too narrow.
# Primary run (full_refit), main universes (k 2-6) with an estimable RRR.
# Output: output/results/sens_site_re.rds, output/summaries/summary_sens_site_re.txt,
#         output/figures/7_sensitivity/s1_site_random_effects.png
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(lme4); library(parallel) })

args  <- commandArgs(trailingOnly = TRUE)
cores <- as.integer(args[1]); if (is.na(cores)) cores <- 8
n_test <- as.integer(args[2])   # optional: run on a few universes only (smoke test)
figd <- fig_dir("7_sensitivity")
sumf <- file.path(SUMM, "summary_sens_site_re.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

res  <- readRDS(file.path(RES, "results_full_refit.rds"))
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race() |>
  transmute(src_subject_id, income_halving = -log2(income), age_bl, sex = factor(sex),
            site = factor(site), family = factor(family), race)
univ <- res$fits |> filter(!supplement, !is.na(est)) |> select(universe_id, method, scale, k, family_fork = family, K, rrr, lo, hi)
if (!is.na(n_test)) { set.seed(1); univ <- univ[sort(sample(nrow(univ), n_test)), ] }

fit_one <- function(i) {
  u <- univ[i, ]
  a <- res$assign[[u$universe_id]]
  dat <- inner_join(a, wide, by = "src_subject_id") |> filter(class %in% c(1, u$K), !is.na(race)) |>
    mutate(y = as.integer(class == u$K), site = droplevels(site), race = droplevels(race))
  re <- tryCatch(glmer(y ~ income_halving + age_bl + sex + race + (1 | site), data = dat, family = binomial, nAGQ = 1,
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                 error = function(e) NULL)
  fe <- glm(y ~ income_halving + age_bl + sex + race + site, data = dat, family = binomial)
  out <- tibble(universe_id = u$universe_id, n_binary = nrow(dat),
                est_fe_bin = coef(fe)[["income_halving"]],
                se_fe_bin = sqrt(vcov(fe)["income_halving", "income_halving"]))
  if (is.null(re)) return(mutate(out, est_re = NA_real_, se_re = NA_real_, site_sd = NA_real_, singular = NA, re_msg = "error"))
  vc <- as.data.frame(VarCorr(re))
  mutate(out, est_re = fixef(re)[["income_halving"]], se_re = sqrt(vcov(re)["income_halving", "income_halving"]),
         site_sd = vc$sdcor[vc$grp == "site"], singular = isSingular(re),
         re_msg = paste(re@optinfo$conv$lme4$messages, collapse = "; "))
}

t0 <- Sys.time()
out <- bind_rows(mclapply(seq_len(nrow(univ)), fit_one, mc.cores = cores))
s <- univ |> inner_join(out, by = "universe_id") |>
  mutate(rrr_re = exp(est_re), lo_re = exp(est_re - 1.96 * se_re), hi_re = exp(est_re + 1.96 * se_re),
         rrr_fe_bin = exp(est_fe_bin), lo_fe_bin = exp(est_fe_bin - 1.96 * se_fe_bin),
         scale_lab = factor(c(raw = "Raw", t = "T-score", log = "log(1+raw)")[scale], c("Raw", "T-score", "log(1+raw)")))
saveRDS(s, file.path(RES, "sens_site_re.rds"))

say(sprintf("Site as a random effect (%s, %.1f min). %d universes; glmer failed in %d, singular fit in %d, convergence messages in %d.",
            format(Sys.time(), "%Y-%m-%d %H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(s),
            sum(is.na(s$est_re)), sum(s$singular, na.rm = TRUE), sum(nchar(s$re_msg) > 0, na.rm = TRUE)))
say("Binary logit, highest vs lowest class, adjusted for age, sex and race/ethnicity. RRR columns: primary multinomial",
    "(site fixed); binary with site fixed; binary with site random. Full-sample universes: no family term here, see header.")
summ <- function(d) d |> summarise(n = sum(!is.na(rrr_re)), median_primary = round(median(rrr), 3),
                                   median_bin_fe = round(median(rrr_fe_bin), 3), median_re = round(median(rrr_re, na.rm = TRUE), 3),
                                   above1_primary = round(mean(lo > 1), 3), above1_bin_fe = round(mean(lo_fe_bin > 1), 3),
                                   above1_re = round(mean(lo_re > 1, na.rm = TRUE), 3),
                                   ratio_re_primary = round(median(rrr_re / rrr, na.rm = TRUE), 3),
                                   site_sd = round(median(site_sd, na.rm = TRUE), 3), .groups = "drop")
say(pt(bind_rows(summ(s) |> mutate(scale_lab = "All", .before = 1), s |> group_by(scale_lab) |> summ() |> mutate(scale_lab = as.character(scale_lab)))))
say("\nBy family fork (full = all children, siblings kept; no family term in this model):")
say(pt(s |> group_by(family_fork, scale_lab) |> summ()))

S_COL <- c(Raw = "#2a78d6", `T-score` = "#eb6834", `log(1+raw)` = "#1baf7a")
INK2 <- "#52514e"; GRID <- "#ecebe8"
p <- ggplot(filter(s, !is.na(rrr_re)), aes(rrr, rrr_re, color = scale_lab)) +
  geom_abline(color = INK2, linewidth = 0.4, linetype = "22") +
  geom_point(size = 1.3, alpha = 0.6) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = S_COL, guide = "none") +
  labs(x = "RRR, primary model (multinomial, site fixed)", y = "RRR, binary logit with site as a random effect",
       title = "Site as a random effect instead of fixed effects",
       subtitle = sprintf("Every universe (k 2-6, n = %d), highest vs lowest class, adjusted for age, sex and race/ethnicity. Dashed line: no change.",
                          sum(!is.na(s$rrr_re)))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", color = NA))
ggsave(file.path(figd, "s1_site_random_effects.png"), p, width = 10, height = 4.2, dpi = 300)
say("\nFigure written to ", figd)
