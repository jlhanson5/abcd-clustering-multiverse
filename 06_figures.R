# =============================================================================
# 06_figures.R: figures and summary tables.
#   Rscript 06_figures.R pilot|full|full_refit
# Main (k 2-6):
#   fig1_spec_curve          poverty RRR (highest vs lowest class, per halving of
#                            income), every universe, sorted; decision grid below
#   fig2_agreement_methods   mean ARI between methods, matched on data config and k
#   fig3_agreement_full      ARI between every pair of universes
#   fig4_alluvial            same children, 8 methods, one data config at k = 4
#   fig5_k_selection         k chosen by BIC / silhouette, with k-max 5 vs 6
#   fig6_diagnostics         failures, empty classes, smallest class (5% line)
#   fig7_decision_importance share of variance in log RRR explained by each fork
# Supplement:
#   figS1_class_gradient     every class's log RRR vs its mean CBCL, by method
#   figS2_k_extension        poverty RRR at k 2-8 for the five cheap methods
# Tables: output/summaries/summary_<mode>.txt
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(patchwork); library(ggalluvial) })

mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "pilot"
res  <- readRDS(file.path(RES, paste0("results_", mode, ".rds")))
agr  <- readRDS(file.path(RES, paste0("agreement_", mode, ".rds")))
figd <- fig_dir(mode)
sumf <- file.path(SUMM, paste0("summary_", mode, ".txt")); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm", "gbtm_link")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM", gbtm_link = "GBTM (link)")
BLUE <- "#2a78d6"; GRAY <- "#9a9893"; INK <- "#0b0b0b"; INK2 <- "#52514e"
SEQ  <- c("#e8f1fb", "#b5d2f2", "#7fb0e8", "#4a8ddc", "#2a78d6", "#1a5aa8", "#0f3d75")

theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ecebe8", linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"),
        plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        plot.background = element_rect(fill = "white", color = NA)))
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)

fits <- res$fits |> mutate(method = factor(method, METHODS))
main <- fits |> filter(!supplement)
ok   <- main |> filter(!is.na(est))

say(sprintf("ABCD clustering multiverse, %s run (%s)", mode, format(Sys.time(), "%Y-%m-%d %H:%M")))
say(sprintf("Main universes (k 2-6): %d. Fit errors: %d. Empty-class (degenerate): %d. Smallest class < 5%%: %d. With poverty estimate: %d.",
            nrow(main), sum(!is.na(main$error)), sum(main$degenerate, na.rm = TRUE),
            sum(main$small_class, na.rm = TRUE), nrow(ok)))

# ---- Fig 1: specification curve ------------------------------------------------
sc <- ok |>
  arrange(est) |>
  mutate(rank = row_number(),
         sig = case_when(lo > 1 ~ "CI above 1", hi < 1 ~ "CI below 1", TRUE ~ "CI includes 1"))
say(sprintf("Poverty RRR, highest vs lowest class, per halving of income: median %.2f (IQR %.2f-%.2f), range %.2f-%.2f.",
            median(sc$rrr), quantile(sc$rrr, .25), quantile(sc$rrr, .75), min(sc$rrr), max(sc$rrr)))
say(sprintf("CI above 1: %d of %d (%.0f%%). CI below 1: %d. CI includes 1: %d.",
            sum(sc$sig == "CI above 1"), nrow(sc), 100 * mean(sc$sig == "CI above 1"),
            sum(sc$sig == "CI below 1"), sum(sc$sig == "CI includes 1")))
for (km in 5:6) {
  s <- filter(sc, k <= km)
  say(sprintf("  k-max %d: %d universes, median RRR %.2f, CI above 1 in %.0f%%.",
              km, nrow(s), median(s$rrr), 100 * mean(s$sig == "CI above 1")))
}

top <- ggplot(sc, aes(rank, rrr)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi, color = sig), linewidth = 0.25, alpha = 0.6) +
  geom_point(aes(color = sig), size = 0.7) +
  scale_color_manual(values = c("CI above 1" = BLUE, "CI includes 1" = GRAY, "CI below 1" = "#e34948"), name = NULL) +
  scale_y_log10() +
  labs(x = NULL, y = "RRR per halving of income\n(highest vs lowest class, log scale)",
       title = "Poverty and membership in the highest-symptom class, across the multiverse",
       subtitle = sprintf("%d universes (k = 2-6), sorted by effect size. Multinomial logit, adjusted for age, sex, site and race/ethnicity.", nrow(sc))) +
  theme(axis.text.x = element_blank(), legend.position = "top", legend.justification = "left")

grid_df <- sc |>
  transmute(rank, Method = M_LAB[as.character(method)], k = paste0("k=", k),
            Scale = c(raw = "raw", t = "T-score", log = "log(1+raw)")[scale],
            Time = time, Family = c(one = "one/family", full = "full sample")[family],
            Sample = c(ge3_imp = ">=3 waves, imputed", ge3_native = ">=3 waves, native", all4 = "all 4 waves")[config],
            `Small class` = ifelse(small_class, "<5% class", "none <5%")) |>
  pivot_longer(-rank, names_to = "fork", values_to = "level") |>
  mutate(fork = factor(fork, c("Method", "k", "Scale", "Time", "Family", "Sample", "Small class")))
lvl_order <- grid_df |> distinct(fork, level) |> arrange(fork, level) |> pull(level) |> unique()
grid_df$level <- factor(grid_df$level, rev(lvl_order))
bot <- ggplot(grid_df, aes(rank, level)) +
  geom_point(shape = 15, size = 0.6, color = INK2) +
  facet_grid(fork ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Universe (sorted by RRR)", y = NULL) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
        axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 6.5))
save(top / bot + plot_layout(heights = c(1, 2.2)), "fig1_spec_curve", 11, 10)

# ---- Fig 2/3: agreement ------------------------------------------------------------
U <- agr$universes |> mutate(method = factor(method, METHODS))
mainU <- which(!U$supplement)
ARI <- agr$ARI[mainU, mainU]; Um <- U[mainU, ]
key <- paste(Um$scale, Um$family, Um$config, Um$k)
pairs <- expand_grid(i = seq_len(nrow(Um)), j = seq_len(nrow(Um))) |>
  filter(key[i] == key[j], i != j) |>
  mutate(m1 = Um$method[i], m2 = Um$method[j], ari = ARI[cbind(i, j)])
mm <- pairs |> group_by(m1, m2) |> summarise(ari = mean(ari), .groups = "drop")
nodiag <- setdiff(METHODS, as.character(mm$m1[mm$m1 == mm$m2]))
mm <- bind_rows(mm, tibble(m1 = factor(nodiag, METHODS), m2 = factor(nodiag, METHODS), ari = NA))
p2 <- ggplot(mm, aes(m1, m2, fill = ari)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(is.na(ari), "", sprintf("%.2f", ari)),
                color = ari > 0.55), size = 3.2) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = SEQ, limits = c(0, 1), na.value = "white", name = "Mean ARI") +
  scale_x_discrete(labels = M_LAB) + scale_y_discrete(labels = M_LAB) +
  labs(x = NULL, y = NULL, title = "How often do two methods put the same children together?",
       subtitle = "Adjusted Rand index between methods, averaged over matched data configurations and k (2-6).\nWithin a method, GBTM/GMM's two time metrics are compared with each other.") +
  coord_equal() + theme(panel.grid = element_blank())
same <- pairs |> filter(m1 == m2) |> group_by(m1) |> summarise(ari = mean(ari))
save(p2, "fig2_agreement_methods", 7.5, 6.5)
say("Mean ARI between methods (matched config and k):")
say(paste(capture.output(print(mm |> filter(!is.na(ari), as.integer(m1) < as.integer(m2)) |>
                                 arrange(ari) |> mutate(ari = round(ari, 2)), n = 30)), collapse = "\n"))

ord <- order(Um$method, Um$k, Um$scale, Um$family, Um$config, Um$time)
hm <- as_tibble(as.table(ARI[ord, ord]), .name_repair = ~ c("a", "b", "ari")) |>
  mutate(a = factor(a, rownames(ARI)[ord]), b = factor(b, rownames(ARI)[ord]))
brk <- cumsum(table(Um$method[ord])); brk <- brk[-length(brk)] + 0.5
mids <- (c(0, head(cumsum(table(Um$method[ord])), -1)) + cumsum(table(Um$method[ord]))) / 2
p3 <- ggplot(hm, aes(a, b, fill = ari)) + geom_raster() +
  geom_vline(xintercept = brk, color = "white", linewidth = 0.6) +
  geom_hline(yintercept = brk, color = "white", linewidth = 0.6) +
  scale_fill_gradientn(colors = SEQ, limits = c(-0.05, 1), name = "ARI") +
  scale_x_discrete(breaks = rownames(ARI)[ord][round(mids)], labels = M_LAB[METHODS[table(Um$method) > 0]]) +
  scale_y_discrete(breaks = rownames(ARI)[ord][round(mids)], labels = M_LAB[METHODS[table(Um$method) > 0]]) +
  labs(x = NULL, y = NULL, title = "Agreement across the full multiverse",
       subtitle = "ARI for every pair of universes (k = 2-6), grouped by method, then k") +
  coord_equal() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
save(p3, "fig3_agreement_full", 8.5, 8)

# ---- Fig 4: alluvial ------------------------------------------------------------------
al_ids <- Um |> filter(k == 4, scale == "raw", family == "one", config == "ge3_imp", time == "wave") |>
  arrange(method)
al <- bind_rows(lapply(seq_len(nrow(al_ids)), \(i)
        res$assign[[al_ids$universe_id[i]]] |> mutate(method = al_ids$method[i]))) |>
  mutate(method = factor(M_LAB[as.character(method)], M_LAB))
first <- al |> filter(method == M_LAB[as.character(al_ids$method[1])]) |> select(src_subject_id, start = class)
al <- left_join(al, first, by = "src_subject_id") |> mutate(class = factor(class), start = factor(start))
p4 <- ggplot(al, aes(x = method, stratum = class, alluvium = src_subject_id, fill = start)) +
  geom_flow(alpha = 0.55, width = 0.25) +
  geom_stratum(width = 0.25, fill = "white", color = INK2, linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = INK) +
  scale_fill_manual(values = SEQ[c(2, 4, 5, 7)], na.translate = FALSE,
                    name = paste0(M_LAB[as.character(al_ids$method[1])], " class")) +
  labs(x = NULL, y = "Children",
       title = "The same children, sorted into 4 classes by each method",
       subtitle = paste0("Raw CBCL total, one child per family, >=3 waves (imputed), k = 4. Classes numbered 1 (lowest) to 4 (highest) mean CBCL.\nFlows are colored by the child's k-means class.",
                         if (length(missing_m <- setdiff(METHODS, as.character(al_ids$method))))
                           paste0(" Not shown (empty class or failed fit at k = 4): ", paste(M_LAB[missing_m], collapse = ", "), ".") else "")) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.x = element_blank())
save(p4, "fig4_alluvial", 11, 6)

# ---- Fig 5: k selection ------------------------------------------------------------
sel <- main |> filter(is.na(error)) |>
  group_by(method, task_id) |>
  group_modify(\(d, g) {
    out <- list()
    for (km in 5:6) {
      dd <- filter(d, k <= km, !degenerate)
      if (nrow(dd) == 0) next
      out[[length(out) + 1]] <- tibble(criterion = "Silhouette (max)", kmax = km, k = dd$k[which.max(dd$asw)])
      if (any(!is.na(dd$bic))) out[[length(out) + 1]] <- tibble(criterion = "BIC (min)", kmax = km, k = dd$k[which.min(dd$bic)])
    }
    bind_rows(out)
  }) |> ungroup()
selc <- sel |> count(method, criterion, kmax, k) |>
  group_by(method, criterion, kmax) |> mutate(share = n / sum(n)) |> ungroup() |>
  mutate(panel = paste0(criterion, ", k-max ", kmax))
p5 <- ggplot(selc, aes(factor(k), method, fill = share)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n, color = share > 0.55), size = 3) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = SEQ, limits = c(0, 1), name = "Share of\nconfigs") +
  scale_y_discrete(labels = M_LAB, limits = rev) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Number of classes selected", y = NULL,
       title = "Which k does each criterion pick?",
       subtitle = "Count of data configurations selecting each k. BIC exists only for the likelihood-based methods.\nEmpty-class solutions are excluded from selection.") +
  theme(panel.grid = element_blank())
save(p5, "fig5_k_selection", 12, 4.5)
say("k selected (count of configs), by method x criterion x k-max:")
say(paste(capture.output(print(selc |> select(method, criterion, kmax, k, n) |>
                                 pivot_wider(names_from = k, values_from = n, values_fill = 0, names_prefix = "k"), n = 40)), collapse = "\n"))

# ---- Fig 6: diagnostics ---------------------------------------------------------------
dg <- main |> group_by(method, k) |>
  summarise(n = n(), failed = mean(!is.na(error) | (!is.na(converged) & !converged)),
            empty = mean(degenerate, na.rm = TRUE),
            small = mean(small_class, na.rm = TRUE), .groups = "drop") |>
  pivot_longer(c(failed, empty, small), names_to = "what", values_to = "share") |>
  mutate(what = factor(c(failed = "Failed or not converged", empty = "Empty class (k not reached)",
                         small = "Smallest class < 5%")[what],
                       c("Failed or not converged", "Empty class (k not reached)", "Smallest class < 5%")))
p6a <- ggplot(dg, aes(factor(k), method, fill = share)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * share), color = share > 0.55), size = 2.8) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = SEQ, limits = c(0, 1), name = "Share of\nuniverses") +
  scale_y_discrete(labels = M_LAB, limits = rev) + facet_wrap(~ what, nrow = 1) +
  labs(x = "k", y = NULL, title = "Where the methods struggle") + theme(panel.grid = element_blank())
p6b <- ggplot(filter(main, is.na(error)), aes(factor(k), min_prop)) +
  geom_hline(yintercept = CLASS_FLOOR, color = "#e34948", linewidth = 0.4, linetype = "22") +
  geom_jitter(width = 0.2, height = 0, size = 0.6, alpha = 0.5, color = BLUE) +
  facet_wrap(~ method, nrow = 1, labeller = as_labeller(M_LAB)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "k", y = "Smallest class", subtitle = "Dashed line: 5% floor (solutions below it are flagged, not dropped)")
save(p6a / p6b + plot_layout(heights = c(1, 1)), "fig6_diagnostics", 12, 7.5)

# ---- Fig 7: decision importance -----------------------------------------------------------
di <- ok |> mutate(across(c(method, scale, time, family, config), factor), kf = factor(k),
                   method_time = interaction(method, time, drop = TRUE))
fit_all <- lm(est ~ method + kf + scale + time + family + config, data = di)
ss_tot <- sum((di$est - mean(di$est))^2)
imp <- tibble(fork = c("method", "kf", "scale", "time", "family", "config")) |>
  mutate(share = map_dbl(fork, \(f) {
    red <- update(fit_all, as.formula(paste(". ~ . -", f)))
    (deviance(red) - deviance(fit_all)) / ss_tot
  }),
  label = c(method = "Method", kf = "Number of classes (k)", scale = "Scale",
            time = "Time metric", family = "Family rule", config = "Sample / missing data")[fork])
say(sprintf("Decision importance (unique share of variance in log RRR; all forks R2 = %.2f):", summary(fit_all)$r.squared))
say(paste(sprintf("  %-24s %.3f", imp$label, imp$share), collapse = "\n"))
p7 <- ggplot(imp, aes(share, reorder(label, share))) +
  geom_col(fill = BLUE, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * share)), hjust = -0.15, size = 3, color = INK) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Unique share of variance in log RRR", y = NULL,
       title = "Which analytic choices move the poverty estimate?",
       subtitle = sprintf("Drop-one comparison from a linear model of log RRR on all forks (R² = %.2f)", summary(fit_all)$r.squared)) +
  theme(panel.grid.major.y = element_blank())
save(p7, "fig7_decision_importance", 7.5, 4)

# ---- GMM refit: sensitivity to flagged solutions (refit modes only) -----------------------
if ("gmm_flag" %in% names(ok)) {
  ci <- \(d) sprintf("%d universes, median RRR %.2f, CI above 1 in %.0f%%, below 1 in %.0f%%",
                     nrow(d), median(d$rrr), 100 * mean(d$lo > 1), 100 * mean(d$hi < 1))
  g <- filter(ok, method == "gmm")
  say("GMM refit, poverty RRR by flag:")
  for (fl in sort(unique(g$gmm_flag))) say(sprintf("  %-26s %s", fl, ci(filter(g, gmm_flag == fl))))
  say(sprintf("  %-26s %s", "all GMM", ci(g)))
  rep_only <- filter(g, gmm_flag == "replicated", !degenerate)
  say(sprintf("  %-26s %s", "replicated, no empty class", ci(rep_only)))
  say(sprintf("Whole multiverse: all %s.", ci(ok)))
  say(sprintf("Whole multiverse without GMM k 5-6: %s.", ci(filter(ok, !(method == "gmm" & k >= 5)))))
  say(sprintf("Whole multiverse, GMM replicated only: %s.", ci(filter(ok, method != "gmm" | (gmm_flag == "replicated" & !degenerate)))))
}

# ---- Supplement S1: class gradient ------------------------------------------------------------
cg <- res$contrasts |>
  inner_join(select(fits, universe_id, method, supplement), by = "universe_id") |>
  filter(!supplement, is.finite(est), se < 5)
pS1 <- ggplot(cg, aes(class_mean, est)) +
  geom_hline(yintercept = 0, color = INK2, linewidth = 0.4) +
  geom_point(size = 0.6, alpha = 0.35, color = BLUE) +
  geom_smooth(method = "loess", se = FALSE, color = INK, linewidth = 0.5, formula = y ~ x) +
  facet_wrap(~ method, nrow = 2, labeller = as_labeller(M_LAB)) +
  labs(x = "Class mean CBCL total (raw, averaged over BL-Y3)",
       y = "log RRR per halving of income\n(vs lowest class)",
       title = "Does the poverty effect follow symptom level, whatever the method?",
       subtitle = "Each point is one non-reference class in one universe (k = 2-6). Points with SE > 5 (separation) omitted.")
save(pS1, "figS1_class_gradient", 11, 6)

# ---- Supplement S2: k extension ------------------------------------------------------------------
ext <- fits |> filter(method %in% c("kmeans", "pam", "ward", "llpa", "kml"), !is.na(est))
pS2 <- ggplot(ext, aes(factor(k), rrr)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.4) +
  geom_jitter(aes(color = supplement), width = 0.15, height = 0, size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c(`FALSE` = GRAY, `TRUE` = BLUE), labels = c("main (k 2-6)", "supplement (k 7-8)"), name = NULL) +
  scale_y_log10() +
  facet_wrap(~ method, nrow = 1, labeller = as_labeller(M_LAB)) +
  labs(x = "k", y = "RRR per halving of income\n(highest vs lowest, log scale)",
       title = "Extending k to 8 for the five cheap methods") +
  theme(legend.position = "top", legend.justification = "left")
save(pS2, "figS2_k_extension", 11, 4.5)

say("Figures written to ", figd)
