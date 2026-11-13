# =============================================================================
# 11_refit_figures.R: additional figures for the primary run (full_refit).
#   Rscript 11_refit_figures.R [mode]        # default full_refit
# Written to fig_dir(mode) (output/figures/1_spec_curve for full_refit) next to 06_figures.R's figures (main, k 2-6):
#   fig8_robustness        poverty RRR median and IQR, share with CI above 1, for
#                          subsets of the multiverse (by fork, without flagged GMM,
#                          and the main run for comparison)
#   fig9_rrr_by_method     every universe's RRR, by method and scale
#   fig10_top_class        size of the highest class vs RRR: why scale matters
#   fig11_trajectories     class mean raw CBCL by wave, anchor configuration, k = 3
#   fig12_gmm_refit        GMM k 2-4: main-run vs refit RRR, and how often starts
#                          converge (R/gmm_refit.R)
#   fig13_race_sensitivity every universe's RRR with race/ethnicity adjustment
#                          (primary) vs without (age, sex, site only)
# Tables: output/summaries/summary_<mode>_extra.txt
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(patchwork))

mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "full_refit"
res  <- readRDS(file.path(RES, paste0("results_", mode, ".rds")))
main_run <- readRDS(file.path(RES, "results_full.rds"))$fits
figd <- fig_dir(mode)
sumf <- file.path(SUMM, paste0("summary_", mode, "_extra.txt")); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm", "gbtm_link")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM", gbtm_link = "GBTM (link)")
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)")
# Reference palette (dataviz skill): categorical slots 1-3 for the three scales;
# CI status as in 06_figures.R; ordinal blue ramp steps 250-650 for k / class order.
S_COL <- c(Raw = "#2a78d6", `T-score` = "#eb6834", `log(1+raw)` = "#1baf7a")
SIG_COL <- c("CI above 1" = "#2a78d6", "CI includes 1" = "#9a9893", "CI below 1" = "#e34948")
ORD <- c("#86b6ef", "#3987e5", "#256abf", "#184f95", "#104281")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"

theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"),
        plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        legend.text = element_text(color = INK2), legend.title = element_text(color = INK2),
        plot.background = element_rect(fill = "white", color = NA)))
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
sig_of <- function(lo, hi) factor(case_when(lo > 1 ~ "CI above 1", hi < 1 ~ "CI below 1", TRUE ~ "CI includes 1"),
                                  names(SIG_COL))
prep <- function(f) f |> filter(!supplement, !is.na(est)) |>
  mutate(method = factor(method, METHODS), scale_lab = factor(S_LAB[scale], S_LAB), sig = sig_of(lo, hi))
ok <- prep(res$fits)
mr <- prep(main_run)

# ---- Fig 8: robustness across subsets ---------------------------------------------------
subsets <- list(
  "All universes" = ok,
  "Main run (before refits)" = mr,
  "Without GMM k 5-6" = filter(ok, !(method == "gmm" & k >= 5)),
  "GMM: replicated, no empty class only" = filter(ok, method != "gmm" | (gmm_flag %in% "replicated" & !degenerate)),
  "Without mixture models (GMM, Gaussian mixture)" = filter(ok, !method %in% c("gmm", "llpa")),
  "No class under 5%" = filter(ok, !small_class))
for (s in names(S_LAB)) subsets[[paste("Scale:", S_LAB[[s]])]] <- filter(ok, scale == s)
for (kk in 2:6) subsets[[paste0("k = ", kk)]] <- filter(ok, k == kk)
for (m in METHODS) subsets[[paste("Method:", M_LAB[[m]])]] <- filter(ok, method == m)
subsets[["Family: one per family"]] <- filter(ok, family == "one")
subsets[["Family: full sample"]] <- filter(ok, family == "full")
grp <- c(rep("Overall", 6), rep("Scale", 3), rep("k", 5), rep("Method", length(METHODS)), rep("Family", 2))
rob <- bind_rows(lapply(names(subsets), \(n) { d <- subsets[[n]]
  tibble(subset = n, n = nrow(d), med = median(d$rrr), q1 = quantile(d$rrr, .25), q3 = quantile(d$rrr, .75),
         above = mean(d$sig == "CI above 1"), below = mean(d$sig == "CI below 1")) })) |>
  mutate(group = factor(grp, unique(grp)), subset = factor(subset, rev(subset)))
say("Robustness of the poverty RRR across subsets (main universes, k 2-6):")
say(paste(capture.output(print(mutate(rob, across(c(med, q1, q3, above, below), \(x) round(x, 2))) |>
                                 select(group, subset, n, med, q1, q3, above, below), n = 40, width = 200)), collapse = "\n"))
p8 <- ggplot(rob, aes(y = subset)) +
  geom_vline(xintercept = 1, color = INK2, linewidth = 0.4) +
  geom_linerange(aes(xmin = q1, xmax = q3), color = "#2a78d6", linewidth = 1.6, alpha = 0.45) +
  geom_point(aes(x = med), color = "#2a78d6", size = 2.2) +
  geom_text(aes(x = 1.5, label = sprintf("%.0f%% above, %.0f%% below  (n = %d)", 100 * above, 100 * below, n)),
            hjust = 0, size = 2.9, color = INK2) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.9, 1, 1.1, 1.2, 1.3, 1.4)) +
  coord_cartesian(xlim = c(0.85, 2.3)) +
  labs(x = "RRR per halving of income (highest vs lowest class, log scale)", y = NULL,
       title = "The poverty gradient holds across subsets of the multiverse",
       subtitle = "Dot: median RRR. Bar: interquartile range. Text: share of universes whose 95% CI is above / below 1.") +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), panel.grid.major.y = element_blank(),
        plot.margin = margin(14, 10, 5, 5))
save(p8, "fig8_robustness", 9, 8.5)

# ---- Fig 9: every universe by method and scale ----------------------------------------------
# Axes are zoomed to the bulk; universes outside are named in the subtitle.
XL <- c(0.85, 1.8)
off <- filter(ok, rrr < XL[1] | rrr > XL[2])
off_txt <- if (nrow(off)) paste0("Off the scale: ", paste(sprintf("%s %s k=%d (RRR %.2f, top class %d children)",
  M_LAB[as.character(off$method)], S_LAB[off$scale], off$k, off$rrr, as.integer(off$top_n)), collapse = "; "), ".") else ""
set.seed(1)
p9 <- ggplot(ok, aes(rrr, method)) +
  geom_vline(xintercept = 1, color = INK2, linewidth = 0.4) +
  geom_jitter(aes(color = sig, shape = sig), height = 0.22, width = 0, size = 1.3, alpha = 0.75) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 6, color = INK) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_color_manual(values = SIG_COL, name = NULL) +
  scale_shape_manual(values = c("CI above 1" = 16, "CI includes 1" = 1, "CI below 1" = 17), name = NULL) +
  scale_x_log10(breaks = c(0.9, 1, 1.2, 1.4, 1.6)) + scale_y_discrete(labels = M_LAB, limits = rev) +
  coord_cartesian(xlim = XL) +
  labs(x = "RRR per halving of income (log scale)", y = NULL,
       title = "Poverty RRR in every universe, by method and scale",
       subtitle = paste0("Each point is one universe (k = 2-6, all data configurations). Black tick: median.\n", off_txt)) +
  theme(legend.position = "top", legend.justification = "left")
save(p9, "fig9_rrr_by_method", 11, 4.8)

# ---- Fig 10: highest-class size vs RRR -------------------------------------------------------
tc <- ok |> mutate(top_share = top_n / n)
say("\nHighest class, by scale: median share of children, median mean CBCL, median RRR:")
say(paste(capture.output(print(tc |> group_by(scale_lab) |>
  summarise(n = n(), top_share = round(median(top_share), 3), top_mean = round(median(top_mean), 1),
            rrr = round(median(rrr), 2), .groups = "drop"))), collapse = "\n"))
p10 <- ggplot(tc, aes(top_share, rrr, color = scale_lab, shape = scale_lab)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.4) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_x_log10(labels = scales::percent_format(accuracy = 1), breaks = c(.01, .03, .1, .3, .6, .9)) +
  scale_y_log10(breaks = c(0.9, 1, 1.2, 1.4, 1.6)) +
  coord_cartesian(xlim = c(0.01, 0.95), ylim = XL) +
  scale_color_manual(values = S_COL, name = "Scale") + scale_shape_manual(values = c(16, 15, 17), name = "Scale") +
  labs(x = "Size of the highest-symptom class (share of children, log scale)",
       y = "RRR per halving of income (log scale)",
       title = "Small severe classes carry the poverty gradient",
       subtitle = paste0("Each point is one universe. On log(1+raw) the \"highest\" class is often a majority group, and the RRR falls to about 1.\n", off_txt)) +
  theme(legend.position = "top", legend.justification = "left")
save(p10, "fig10_top_class", 8, 5.5)

# ---- Fig 11: class trajectories, anchor configuration, k = 3 --------------------------------
wide <- readRDS(file.path(DER, "wide.rds"))
R <- wide |> select(src_subject_id, all_of(paste0("imp_r_", WAVES)))
anc <- res$fits |> filter(family == "one", config == "ge3_imp", k == 3,
                          time == ifelse(method %in% c("gbtm", "gmm", "gbtm_link"), "age", "wave"))
tr <- bind_rows(lapply(seq_len(nrow(anc)), function(i) {
  a <- res$assign[[anc$universe_id[i]]]; if (is.null(a)) return(NULL)
  inner_join(a, R, by = "src_subject_id") |>
    group_by(class) |> summarise(n = n(), across(starts_with("imp_r_"), mean), .groups = "drop") |>
    mutate(share = n / sum(n), method = anc$method[i], scale = anc$scale[i])
})) |>
  pivot_longer(starts_with("imp_r_"), names_to = "wave", values_to = "cbcl") |>
  mutate(wave = factor(sub("imp_r_", "", wave), WAVES), method = factor(method, METHODS),
         scale_lab = factor(S_LAB[scale], S_LAB), class = factor(class))
lab11 <- tr |> filter(wave == WAVES[length(WAVES)])
p11 <- ggplot(tr, aes(wave, cbcl, group = class, color = class)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  geom_text(data = lab11, aes(label = sprintf("%.0f%%", 100 * share)), hjust = -0.25, size = 2.4,
            color = INK2, show.legend = FALSE, check_overlap = TRUE) +
  facet_grid(scale_lab ~ method, labeller = labeller(method = M_LAB)) +
  scale_color_manual(values = ORD[c(1, 3, 5)], name = "Class (1 = lowest)") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
  labs(x = "Wave", y = "Mean raw CBCL total",
       title = "What the three classes look like, by method and scale",
       subtitle = "One child per family, >= 3 waves (imputed), k = 3. Labels: class size. Raw CBCL shown for every scale.") +
  theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(size = 6.5))
save(p11, "fig11_trajectories", 13, 6.5)

# ---- Fig 12: GMM refit ----------------------------------------------------------------------
g <- inner_join(mr |> filter(method == "gmm", k <= 4) |> select(universe_id, k, scale_lab, rrr_main = rrr),
                ok |> filter(method == "gmm", k <= 4) |> select(universe_id, rrr_refit = rrr, gmm_flag),
                by = "universe_id") |>
  mutate(kf = factor(paste0("k = ", k)), flag = factor(ifelse(gmm_flag == "replicated", "Replicated", "Best found once")))
p12a <- ggplot(g, aes(rrr_main, rrr_refit)) +
  geom_abline(color = INK2, linewidth = 0.4, linetype = "22") +
  geom_hline(yintercept = 1, color = GRID) + geom_vline(xintercept = 1, color = GRID) +
  geom_point(aes(color = kf, shape = flag), size = 1.8, alpha = 0.85) +
  facet_wrap(~ scale_lab, nrow = 1, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = ORD[c(1, 3, 5)], name = NULL) +
  scale_shape_manual(values = c(Replicated = 16, `Best found once` = 1), name = NULL) +
  labs(x = "RRR, main run (gridsearch)", y = "RRR, refit (best of 50 converged starts)",
       title = "GMM: how the refit moved each universe",
       subtitle = "k = 2-4, 36 data configurations each. Dashed line: no change.") +
  theme(legend.position = "top", legend.justification = "left")
st <- bind_rows(lapply(list.files(file.path(OUT, "gmm_refit"), "\\.rds$", full.names = TRUE), \(f) {
  r <- readRDS(f); tibble(task_id = r$task_id, k = r$k, conv = r$conv) }))
cv <- st |> mutate(scale_lab = factor(S_LAB[sub("gmm_([a-z]+)_.*", "\\1", task_id)], S_LAB),
                   time = ifelse(grepl("_age_", task_id), "Age", "Wave")) |>
  group_by(scale_lab, time, k) |> summarise(n = n(), conv = sum(conv %in% 1), .groups = "drop") |>
  mutate(share = conv / n, row = paste(scale_lab, "by", tolower(time)))
say("\nGMM refit: converged starts by scale, time metric and k:")
say(paste(capture.output(print(select(cv, scale_lab, time, k, conv, n))), collapse = "\n"))
p12b <- ggplot(cv, aes(factor(k), row, fill = share)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%d of %d", conv, n), color = share > 0.55), size = 3) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = c("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281"), limits = c(0, 1),
                       labels = scales::percent, name = "Converged") +
  scale_y_discrete(limits = rev) +
  labs(x = "k", y = NULL, title = "How often a GMM start converges",
       subtitle = "50 starts per data configuration, 6 configurations per row") +
  theme(panel.grid = element_blank())
save(p12a / p12b + plot_layout(heights = c(1.15, 1)), "fig12_gmm_refit", 10, 8.5)

# ---- Fig 13: race/ethnicity adjustment (primary) vs age, sex, site only ----------------------
rs <- ok |> filter(!is.na(rrr_unadj)) |> mutate(ratio = rrr / rrr_unadj, sig_unadj = sig_of(lo_unadj, hi_unadj))
say("\nRace/ethnicity adjustment: RRR adjusted / RRR without race, by scale (median, IQR), and CI status:")
say(paste(capture.output(print(rs |> group_by(scale_lab) |>
  summarise(n = n(), median_adj = round(median(rrr), 3), median_unadj = round(median(rrr_unadj), 3),
            ratio_median = round(median(ratio), 3), ratio_q25 = round(quantile(ratio, .25), 3),
            ratio_q75 = round(quantile(ratio, .75), 3), above1_adj = round(mean(lo > 1), 3),
            above1_unadj = round(mean(lo_unadj > 1), 3), .groups = "drop"), width = 200)), collapse = "\n"))
say(paste(capture.output(print(table(adjusted = rs$sig, unadjusted = rs$sig_unadj))), collapse = "\n"))
p13 <- ggplot(rs, aes(rrr_unadj, rrr)) +
  geom_abline(color = INK2, linewidth = 0.4, linetype = "22") +
  geom_hline(yintercept = 1, color = GRID) + geom_vline(xintercept = 1, color = GRID) +
  geom_point(aes(color = scale_lab), size = 1.3, alpha = 0.6) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = S_COL, guide = "none") +
  labs(x = "RRR, adjusted for age, sex, site (sensitivity)", y = "RRR, also adjusted for race/ethnicity (primary)",
       title = "Adjusting for race/ethnicity moves the poverty estimate up, not down",
       subtitle = sprintf("Every universe (k 2-6, n = %d), highest vs lowest class. Dashed line: no change.", nrow(rs)))
save(p13, "fig13_race_sensitivity", 10, 4.2)

say("\nFigures written to ", figd)
