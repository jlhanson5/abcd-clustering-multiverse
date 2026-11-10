# =============================================================================
# 15_class_gradient.R: the full poverty gradient across classes (added 2026-11-09).
#   Rscript 15_class_gradient.R
# The primary contrast stays highest vs lowest-mean class. This script reports every
# class against the lowest, so a dip below the lowest class (a near-zero reference
# group whose families are poorer than the mild classes) is visible, not hidden.
# Race-adjusted model (primary). Primary run (full_refit), main universes (k 2-6);
# classes with fewer than MIN_TOP_N children are left out (separation).
# Figures in output/figures/6_class_gradient/:
#   g1_gradient_lines   each universe's RRR vs the lowest class across its classes,
#                       against class mean raw CBCL; scale x k
#   g2_gradient_shape   share of universes (k 3-6) whose gradient rises monotonically,
#                       rises with a dip, or dips below the lowest class; scale x method
#   g3_consensus_gradient  every class of the 16 clusterBMA consensus fits, 95% CI
# Table: output/summaries/summary_class_gradient.txt
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(patchwork))

figd <- fig_dir("6_class_gradient")
sumf <- file.path(SUMM, "summary_class_gradient.txt"); unlink(sumf)
say  <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt   <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM")
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)", all = "All scales")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged", kbal = "k-balanced")
SHAPE <- c("Rises steadily", "Rises, not steadily", "Middle class below lowest", "Middle class below lowest (CI)")
SHAPE_COL <- setNames(c("#2a78d6", "#9cc3f0", "#f0a07a", "#e34948"), SHAPE)
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        legend.text = element_text(color = INK2), legend.title = element_text(color = INK2),
        plot.background = element_rect(fill = "white", color = NA)))
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)

# ---- Data ---------------------------------------------------------------------------------------
res  <- readRDS(file.path(RES, "results_full_refit.rds"))
wide <- readRDS(file.path(DER, "wide.rds")) |>
  mutate(cbcl_mean = rowMeans(across(all_of(paste0("imp_r_", WAVES)))))
main <- res$fits |> filter(!supplement, !top_too_small, universe_id %in% res$contrasts$universe_id) |>
  select(universe_id, method, scale, k)
# Lowest class (the reference): size and mean raw CBCL, from the relabelled assignments
ref <- bind_rows(lapply(main$universe_id, \(u) {
  a <- res$assign[[u]]; m <- wide$cbcl_mean[match(a$src_subject_id[a$class == 1], wide$src_subject_id)]
  tibble(universe_id = u, class = 1L, class_mean = mean(m), n = length(m), rrr = 1, lo = 1, hi = 1)
}))
cx_all <- res$contrasts |> filter(universe_id %in% main$universe_id)
cx <- cx_all |> filter(n >= MIN_TOP_N) |> select(universe_id, class, K, class_mean, n, rrr, lo, hi)
g <- bind_rows(cx, ref |> left_join(distinct(cx, universe_id, K), by = "universe_id")) |>
  inner_join(main, by = "universe_id") |>
  mutate(scale_lab = factor(S_LAB[scale], S_LAB[1:3]), method = factor(method, METHODS)) |>
  arrange(universe_id, class)
say(sprintf("Class gradient, primary run (full_refit), race-adjusted (%s). %d main universes, %d classes vs the lowest.",
            format(Sys.time(), "%Y-%m-%d %H:%M"), n_distinct(g$universe_id), nrow(cx)))
say(sprintf("Left out (class under %d children): %d classes in %d universes.", MIN_TOP_N,
            sum(cx_all$n < MIN_TOP_N), n_distinct(cx_all$universe_id[cx_all$n < MIN_TOP_N])))

# ---- Shape of each universe's gradient (k 3-6: at least one middle class) ------------------------
shape <- g |> filter(K >= 3) |> group_by(universe_id, scale_lab, method, k) |>
  summarise(monotone = all(diff(rrr) >= 0),
            dip = any(rrr[class > 1 & class < max(class)] < 1),
            dip_ci = any(hi[class > 1 & class < max(class)] < 1), .groups = "drop") |>
  mutate(shape = factor(case_when(dip_ci ~ SHAPE[4], dip ~ SHAPE[3], monotone ~ SHAPE[1], TRUE ~ SHAPE[2]), SHAPE))
say("\nShape of the gradient (k 3-6). Steady: RRR never falls from one class to the next. Below lowest:",
    "a middle class has RRR < 1 against the lowest class (CI: its 95% CI is below 1).")
say(pt(shape |> count(scale_lab, shape) |> group_by(scale_lab) |> mutate(pct = round(100 * n / sum(n))) |>
         ungroup() |> select(-n) |> pivot_wider(names_from = shape, values_from = pct, values_fill = 0)))
say("\nBy scale and method (% of universes with a middle class below the lowest, point estimate):")
say(pt(shape |> group_by(scale_lab, method) |> summarise(n = n(), dip_pct = round(100 * mean(dip)),
                                                          dip_ci_pct = round(100 * mean(dip_ci)), .groups = "drop") |>
         pivot_wider(names_from = scale_lab, values_from = c(n, dip_pct, dip_ci_pct))))
say("\nDips below the lowest class by k (all scales):")
say(pt(shape |> group_by(k) |> summarise(n = n(), dip_pct = round(100 * mean(dip)), dip_ci_pct = round(100 * mean(dip_ci)))))

# ---- Gradient by class position --------------------------------------------------------------------
say("\nMedian RRR vs the lowest class, by class position (1 = lowest) and k, per scale:")
say(pt(g |> filter(class > 1) |> group_by(scale_lab, k, class) |>
         summarise(median_rrr = round(median(rrr), 2), .groups = "drop") |>
         pivot_wider(names_from = class, names_prefix = "class_", values_from = median_rrr)))

# ---- g1: gradient lines -------------------------------------------------------------------------------
g1d <- g |> left_join(select(shape, universe_id, shape), by = "universe_id") |>
  mutate(shape = factor(ifelse(K == 2, "k = 2 (one contrast)", as.character(shape)), c(SHAPE, "k = 2 (one contrast)")),
         kf = factor(paste0("k = ", k)))
p1 <- ggplot(g1d, aes(class_mean, rrr, group = universe_id, color = shape)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.3) +
  geom_line(linewidth = 0.35, alpha = 0.55) +
  geom_point(size = 0.7, alpha = 0.7) +
  facet_grid(scale_lab ~ kf) +
  scale_x_log10(breaks = c(3, 10, 30)) + scale_y_log10() +
  scale_color_manual(values = c(SHAPE_COL, `k = 2 (one contrast)` = "#9a9893"), name = NULL, drop = FALSE) +
  labs(x = "Class mean raw CBCL total (log scale)", y = "Poverty RRR vs the lowest class (log scale)",
       title = "The poverty gradient across every class, not just highest vs lowest",
       subtitle = paste("One line per universe (primary run, k 2-6), from the lowest class (RRR = 1) to the highest.",
                        "Adjusted for age, sex, site and race/ethnicity.")) +
  theme(legend.position = "top", legend.justification = "left") +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))
p1 <- p1 + theme(panel.spacing.x = unit(1, "lines"))
save(p1, "g1_gradient_lines", 12, 8)

# ---- g2: shape shares -------------------------------------------------------------------------------
g2d <- shape |> count(scale_lab, method, shape) |> group_by(scale_lab, method) |> mutate(share = n / sum(n)) |> ungroup()
p2 <- ggplot(g2d, aes(share, method, fill = shape)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.4) +
  facet_wrap(~ scale_lab, nrow = 1) +
  scale_x_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_y_discrete(labels = M_LAB, limits = rev) +
  scale_fill_manual(values = SHAPE_COL, name = NULL, drop = FALSE) +
  labs(x = "Share of universes (k 3-6)", y = NULL,
       title = "Where the lowest class is not the least poor",
       subtitle = "Rises steadily: RRR never falls from one class to the next. Below lowest: a middle class has RRR < 1 against the lowest.") +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.y = element_blank(), panel.spacing.x = unit(1.5, "lines"))
save(p2, "g2_gradient_shape", 11, 4.5)

# ---- g3: consensus fits ---------------------------------------------------------------------------------
cp <- readRDS(file.path(RES, "cbma_poverty.rds")) |> filter(model == "adj")
runs <- lapply(sort(list.files(file.path(OUT, "cbma"), "\\.rds$", full.names = TRUE)), readRDS)
cref <- bind_rows(lapply(runs, \(r) {
  m <- wide$cbcl_mean[match(r$ids[r$labels == 1], wide$src_subject_id)]
  tibble(scale = r$scale, weighting = r$weighting, class = 1L, class_mean = mean(m), n = length(m), rrr = 1, lo = 1, hi = 1)
}))
g3d <- bind_rows(select(cp, scale, weighting, class, class_mean, n, rrr, lo, hi), cref) |>
  mutate(scale_lab = factor(S_LAB[scale], S_LAB), weighting = factor(W_LAB[weighting], W_LAB),
         class = factor(class), pct = 100 * n / 8331) |>
  arrange(scale_lab, weighting, class_mean)
say("\nclusterBMA consensus fits, every class vs the lowest (race-adjusted):")
say(pt(g3d |> transmute(scale_lab, weighting, class, pct = round(pct), class_mean = round(class_mean, 1),
                        rrr = round(rrr, 2), ci = ifelse(class == "1", "ref", sprintf("%.2f-%.2f", lo, hi)))))
ORD <- c("#a6c8f0", "#5b9be6", "#2563b8", "#0d3266", "#08203f", "#04101f")
p3 <- ggplot(g3d, aes(class_mean, rrr)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.3) +
  geom_line(aes(group = weighting), color = GRID, linewidth = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, color = "#5b9be6", linewidth = 0.5) +
  geom_point(aes(fill = class, size = pct), shape = 21, color = "white", stroke = 0.4) +
  facet_grid(scale_lab ~ weighting) +
  scale_x_log10() + scale_y_log10() +
  scale_fill_manual(values = ORD, name = "Class (1 = lowest, reference)") +
  guides(fill = guide_legend(override.aes = list(size = 3.5))) +
  scale_size_area(max_size = 4.5, name = "% of children", breaks = c(5, 20, 40)) +
  labs(x = "Class mean raw CBCL total (log scale)", y = "Poverty RRR vs the lowest class (95% CI, log scale)",
       title = "The full gradient in every clusterBMA consensus fit",
       subtitle = "Adjusted for age, sex, site and race/ethnicity. On log(1+raw), equal and k-balanced weights make the lowest class a near-zero group.") +
  theme(legend.position = "top", legend.justification = "left")
save(p3, "g3_consensus_gradient", 12, 8.5)
say("\nFigures written to ", figd)
