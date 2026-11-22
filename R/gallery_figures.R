# =============================================================================
# R/gallery_figures.R (2026-11-22): a CBCL results gallery to choose from, matching the
# cognition and positive affect galleries (../ABCD_cognition_affect_clustering_20261114/
# 12_gallery.R). Draws the figure types this project lacked and copies the strongest
# existing ones, numbered in reading order.
# New: stability by k (real vs null), consensus classes with 95% CI, predicted class
# membership across income, class make-up by income, uncertainty by class, inputs vs
# the consensus. Aggregates only.
#   Rscript R/gallery_figures.R      (after 10, 13, 14, 15, R/main_figures.R, R/cbma_check_figures.R)
# Output: output/figures/gallery_cbcl/
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(patchwork) })

figd <- file.path(FIG, "gallery_cbcl"); dir.create(figd, recursive = TRUE, showWarnings = FALSE)
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        legend.position = "top", legend.justification = "left",
        plot.background = element_rect(fill = "white", color = NA)))
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Long. LPA", kml = "kml", gbtm = "GBTM",
           gbtm_link = "GBTM (link)", gmm = "GMM", gmmv = "Gaussian mixture")
M_COL <- c(`k-means` = "#2563b8", PAM = "#c0611f", Ward = "#7a7975", `Long. LPA` = "#5fa8a0", kml = "#0d3266",
           GBTM = "#b8860b", `GBTM (link)` = "#d4a84a", GMM = "#a33b5c", `Gaussian mixture` = "#8a5bb5")
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)", all = "All scales pooled")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged")
SHOW <- c("raw_ch", "t_ch", "log_ch", "all_equal")          # one fit per scale for the per-fit figures
fit_lab <- function(nm) { p <- strsplit(nm, "_")[[1]]; paste(S_LAB[p[1]], W_LAB[p[2]], sep = ", ") }
REL <- function(aes = "color") {
  f <- if (aes == "color") scale_color_gradientn else scale_fill_gradientn
  f(colours = c("#0d3266", "#5b9be6", "#e0a96d", "#c0611f"), limits = c(0, 1), breaks = c(0, 1),
    labels = c("lowest", "highest"), name = "Class (CBCL)")
}
inc_lab <- function(x) ifelse(x >= 1000, paste0("$", round(x / 1000), "k"), paste0("$", x))

wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()
cbma <- lapply(setNames(file.path(OUT, "cbma", paste0(SHOW, ".rds")), SHOW), readRDS)
child <- function(nm) { r <- cbma[[nm]]; wide[match(r$ids, wide$src_subject_id), ] |> mutate(class = r$labels, uncertainty = r$uncertainty) }

# ---- n06 stability by k, real vs null (k = 2-6 from summary.rds, 7-8 from summary_k78.rds) --
S <- bind_rows(readRDS(file.path(OUT, "stability", "summary.rds"))$pairs,
               readRDS(file.path(OUT, "stability", "summary_k78.rds"))$pairs) |>
  filter(ptype %in% c("subsample", "null_subsample"), variant == "base", method %in% c("kmeans", "kml", "pam", "gbtm"))   # the clusterBMA candidates
st <- S |> group_by(anchor, method, ptype, k) |>
  summarise(med = median(ari, na.rm = TRUE), q1 = quantile(ari, .25, na.rm = TRUE), q3 = quantile(ari, .75, na.rm = TRUE), .groups = "drop") |>
  mutate(Data = factor(ifelse(ptype == "subsample", "Real data", "Null data"), c("Real data", "Null data")),
         Anchor = factor(S_LAB[anchor], S_LAB), Method = factor(M_LAB[method], M_LAB))
p <- ggplot(st, aes(k, med, color = Method, linetype = Data)) +
  geom_hline(yintercept = 0.75, color = "#c0611f", linetype = "22", linewidth = 0.4) +
  geom_ribbon(data = filter(st, Data == "Real data"), aes(ymin = q1, ymax = q3, fill = Method), alpha = 0.1, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.1) + facet_wrap(~ Anchor, nrow = 1) +
  scale_color_manual(values = M_COL, name = NULL) + scale_fill_manual(values = M_COL, guide = "none") +
  scale_linetype_manual(values = c(`Real data` = "solid", `Null data` = "22"), name = NULL) +
  scale_x_continuous(breaks = 2:8) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = "k", y = "Median pairwise subsample ARI", title = "CBCL: stability by k, real vs null data",
       subtitle = "The four clusterBMA candidates. Band: interquartile range over subsample pairs (real data). Dashed orange: the .75 rule. Null: no subgroups (none for GBTM).")
save(p, "06_stability_by_k", 12, 5.5)

# ---- n08 consensus classes with 95% CI ------------------------------------------------------
m <- bind_rows(lapply(SHOW, function(nm) {
  d <- child(nm)
  d |> select(class, all_of(paste0("imp_r_", WAVES))) |> pivot_longer(-class, names_to = "wave", values_to = "y") |>
    group_by(class, wave) |> summarise(mean = mean(y), se = sd(y) / sqrt(n()), n = n(), .groups = "drop") |>
    mutate(share = n / nrow(d), fit = fit_lab(nm), K = max(class))
})) |> mutate(wave = factor(sub("imp_r_", "", wave), WAVES), fit = factor(fit, vapply(SHOW, fit_lab, "")),
              rel = (class - 1) / pmax(1, K - 1))
p <- ggplot(m, aes(wave, mean, group = class, color = rel, fill = rel)) +
  geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  geom_text(data = filter(m, wave == WAVES[length(WAVES)]), aes(label = sprintf("%.0f%%", 100 * share)), hjust = -0.35, size = 2.8, color = INK2) +
  facet_wrap(~ fit, nrow = 1) + REL("color") + REL("fill") + guides(fill = "none") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.8))) +
  labs(x = "Wave", y = "Mean raw CBCL total", title = "CBCL: consensus classes with 95% CI",
       subtitle = "Primary input set (k-means + kml, k = 2-8); one fit per scale. Labels: share of children.")
save(p, "08_consensus_classes_ci", 13, 5)

# ---- n10 / n11 / n12: income and uncertainty ---------------------------------------------------
pp <- list(); comp <- list(); unc <- list()
for (nm in SHOW) {
  d <- child(nm) |> mutate(income_halving = -log2(income), class = factor(class), sex = factor(sex), site = factor(site)) |>
    filter(if_all(all_of(setdiff(all.vars(POV_F$adj), "class")), \(v) !is.na(v)))
  mod <- multinom(POV_F$adj, data = d, trace = FALSE, maxit = 1000, MaxNWts = 5000)
  inc <- d |> distinct(income_halving, income) |> arrange(income)
  pp[[nm]] <- bind_rows(lapply(seq_len(nrow(inc)), function(i) {
    nd <- d; nd$income_halving <- inc$income_halving[i]; pr <- colMeans(predict(mod, newdata = nd, type = "probs"))
    tibble(income = inc$income[i], class = as.integer(names(pr)), p = pr)
  })) |> mutate(fit = fit_lab(nm))
  comp[[nm]] <- d |> dplyr::count(income, class) |> group_by(income) |> mutate(p = n / sum(n), N = sum(n)) |> ungroup() |>
    filter(N >= 30) |> mutate(class = as.integer(as.character(class)), fit = fit_lab(nm))
  unc[[nm]] <- d |> group_by(class) |> summarise(q10 = quantile(uncertainty, .1), q25 = quantile(uncertainty, .25), med = median(uncertainty),
                                             q75 = quantile(uncertainty, .75), q90 = quantile(uncertainty, .9), mean = mean(uncertainty), .groups = "drop") |>
    mutate(class = as.integer(as.character(class)), fit = fit_lab(nm))
}
relf <- function(x) x |> group_by(fit) |> mutate(rel = (class - 1) / pmax(1, max(class) - 1)) |> ungroup() |>
  mutate(fit = factor(fit, vapply(SHOW, fit_lab, "")))
pp <- relf(bind_rows(pp))
p <- ggplot(pp, aes(income, p, color = rel, group = class)) + geom_line(linewidth = 1) + geom_point(size = 1.3) +
  geom_text(data = filter(pp, income == max(income)), aes(label = paste("class", class)), hjust = -0.15, size = 2.7, show.legend = FALSE) +
  facet_wrap(~ fit, nrow = 1) + REL("color") +
  scale_x_log10(labels = inc_lab, expand = expansion(mult = c(0.03, 0.2))) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Family income (log scale)", y = "Predicted share in class",
       title = "CBCL: predicted class membership across family income",
       subtitle = "Race-adjusted multinomial, averaged over the sample's age, sex, site and race/ethnicity. Class 1 = lowest CBCL.")
save(p, "10_predicted_by_income", 14, 5)
comp <- relf(bind_rows(comp))
p <- ggplot(comp, aes(factor(inc_lab(income), inc_lab(sort(unique(income)))), p, fill = rel, group = -class)) +
  geom_col(width = 0.85, color = "white", linewidth = 0.2) + facet_wrap(~ fit, ncol = 2) +
  scale_y_continuous(labels = scales::percent) + REL("fill") +
  labs(x = "Family income (ABCD bracket)", y = "Share of children", title = "CBCL: class make-up by family income",
       subtitle = "Observed shares; primary input set, one fit per scale.") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save(p, "11_classes_by_income", 12, 8)
unc <- relf(bind_rows(unc))
p <- ggplot(unc, aes(factor(class), color = rel)) +
  geom_linerange(aes(ymin = q10, ymax = q90), linewidth = 0.6) + geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 2.6) +
  geom_point(aes(y = med), color = "white", size = 1.6) + geom_point(aes(y = mean), shape = 4, color = INK, size = 1.8) +
  facet_wrap(~ fit, nrow = 1, scales = "free_x") + REL("color") + guides(color = "none") +
  labs(x = "Consensus class (1 = lowest CBCL)", y = "Uncertainty (1 - highest class probability)", title = "CBCL: how sure is each child's class?",
       subtitle = "Thick bar: interquartile range; thin: 10th-90th percentile; white dot: median; cross: mean.")
save(p, "12_uncertainty", 12, 4.5)

# ---- n13 inputs vs the consensus ---------------------------------------------------------
inp <- bind_rows(lapply(SHOW[SHOW != "all_equal"], function(nm) cbma[[nm]]$inputs |> mutate(fit = fit_lab(nm)))) |>
  mutate(Method = factor(M_LAB[method], M_LAB), family = recode(family, one = "One per family", full = "All children"),
         fit = factor(fit, vapply(SHOW, fit_lab, "")))
p <- ggplot(inp, aes(k, ari, color = Method, shape = family, size = weight)) +
  geom_line(aes(group = interaction(Method, family)), linewidth = 0.4, alpha = 0.6) + geom_point(alpha = 0.85) +
  facet_wrap(~ fit, nrow = 1) + scale_color_manual(values = M_COL, name = NULL) + scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_size_area(max_size = 3.5, name = "CH weight") + scale_x_continuous(breaks = 2:8) +
  labs(x = "k of the input solution", y = "ARI with the consensus", title = "CBCL: which inputs does the consensus follow?",
       subtitle = "Each point is one input solution; size is its Calinski-Harabasz weight.")
save(p, "13_inputs_vs_consensus", 12, 5)

# ---- copies of the strongest existing figures --------------------------------------------
copies <- c(`01_spec_curve` = "1_spec_curve/fig1_spec_curve", `02_decision_importance` = "1_spec_curve/fig7_decision_importance",
            `03_rrr_by_method` = "1_spec_curve/fig9_rrr_by_method", `04_top_class` = "1_spec_curve/fig10_top_class",
            `05_stability_main` = "0_main/main4_stability", `07_clusterbma_main` = "0_main/main5_clusterbma",
            `09_gradient_lines` = "6_class_gradient/g1_gradient_lines", `14_consensus_methods_poverty` = "8_consensus_checks/c3_methods_poverty",
            `15_consensus_agreement` = "8_consensus_checks/c4_agreement", `16_permutation` = "7_sensitivity/s2_permutation")
for (nm in names(copies)) file.copy(file.path(FIG, paste0(copies[[nm]], ".png")), file.path(figd, paste0(nm, ".png")), overwrite = TRUE)
cat("Gallery written to", figd, "\n")
