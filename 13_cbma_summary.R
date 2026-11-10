# =============================================================================
# 13_cbma_summary.R: compare the clusterBMA consensus fits from 12_cbma.R.
#   Rscript 13_cbma_summary.R [dir]          # default output/cbma
# Tables: output/summaries/summary_cbma.txt (summary_cbma_smoke.txt for the smoke dir)
# Figures in output/figures/4_clusterbma/:
#   cbma1_weights        weight of every input, by method and k, per scale x weighting
#   cbma2_trajectories   consensus classes: mean raw CBCL by wave, with class sizes
#   cbma3_uncertainty    per-child allocation uncertainty, per scale x weighting
#   cbma4_agreement      ARI of each method's inputs with the consensus
#   cbma5_poverty        poverty RRR (95% CI) for every consensus fit
#   cbma6_weighting_agreement  ARI between weightings on the same scale
#   cbma7_uncertainty_by_class  mean allocation uncertainty per consensus class
# Runs: scale raw / t / log / all x weighting ch / sdbw / equal / avg / kbal (see 12_cbma.R).
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(patchwork); library(mclust) })

dir <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(dir)) dir <- file.path(OUT, "cbma")
smoke <- grepl("smoke", dir)
figd <- fig_dir(if (smoke) "cbma_smoke" else "cbma")
sumf <- file.path(SUMM, if (smoke) "summary_cbma_smoke.txt" else "summary_cbma.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM")
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)", all = "All scales")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged", kbal = "k-balanced")
ORD <- c("#a6c8f0", "#5b9be6", "#2563b8", "#0d3266", "#08203f", "#04101f")  # light to dark; wide steps so 4 classes read apart
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        plot.background = element_rect(fill = "white", color = NA)))
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)

files <- list.files(dir, "\\.rds$", full.names = TRUE)
runs <- setNames(lapply(files, readRDS), sub("\\.rds$", "", basename(files)))
lab <- function(r) tibble(scale = factor(S_LAB[r$scale], S_LAB), weighting = factor(W_LAB[r$weighting], W_LAB))
# Poverty: race-adjusted (primary since 2026-11-09) and unadjusted, from 12b_cbma_poverty.R.
# The smoke dir keeps the unadjusted model stored by 12_cbma.R.
if (!smoke) {
  cp <- readRDS(file.path(RES, "cbma_poverty.rds")) |> filter(top)
  top_of <- function(r, m) filter(cp, scale == r$scale, weighting == r$weighting, model == m) |>
    transmute(K, est, se, rrr, lo, hi, top_n = n, top_mean = class_mean)
  for (nm in names(runs)) {
    runs[[nm]]$poverty <- top_of(runs[[nm]], "adj")
    runs[[nm]]$poverty_unadj <- top_of(runs[[nm]], "unadj")
    stopifnot(nrow(runs[[nm]]$poverty) == 1, nrow(runs[[nm]]$poverty_unadj) == 1)
  }
}
say(sprintf("clusterBMA consensus fits (%s): %d runs. Port %s, %s.", format(Sys.time(), "%Y-%m-%d %H:%M"),
            length(runs), runs[[1]]$port, runs[[1]]$blas))

# ---- Overview ----------------------------------------------------------------------------
ov <- bind_rows(lapply(runs, function(r) bind_cols(lab(r), tibble(
  inputs = nrow(r$inputs), eff_inputs = round(1 / sum(r$inputs$weight^2), 1),
  max_weight = round(max(r$inputs$weight), 3), clusters = length(r$sizes),
  sizes = paste(round(100 * r$sizes / r$n), collapse = "/"),
  mean_uncert = round(mean(r$uncertainty), 3), uncert_over_20 = round(mean(r$uncertainty > .2), 3),
  rrr = if (!is.null(r$poverty$rrr)) round(r$poverty$rrr, 2) else NA,
  ci = if (!is.null(r$poverty$rrr)) sprintf("%.2f-%.2f", r$poverty$lo, r$poverty$hi) else NA,
  top_pct = if (!is.null(r$poverty$top_n)) round(100 * r$poverty$top_n / r$n, 1) else NA,
  top_mean = if (!is.null(r$poverty$top_mean)) round(r$poverty$top_mean, 1) else NA,
  rrr_unadj = if (!is.null(r$poverty_unadj$rrr)) round(r$poverty_unadj$rrr, 2) else NA,
  minutes = round(sum(r$minutes)))))) |> arrange(scale, weighting)
say("\nOverview. eff_inputs = 1 / sum(w^2), the number of equally weighted inputs with the same spread.",
    "sizes: % of children per consensus class, lowest to highest mean CBCL. rrr: poverty RRR, highest vs lowest class,",
    "adjusted for age, sex, site and race/ethnicity (primary); rrr_unadj: without race/ethnicity.")
say(pt(ov))

# ---- Weights by method -----------------------------------------------------------------------
wt <- bind_rows(lapply(runs, \(r) bind_cols(lab(r), r$inputs))) |> mutate(method = factor(method, METHODS))
say("\nTotal weight by method:")
say(pt(wt |> group_by(scale, weighting, method) |> summarise(w = round(sum(weight), 3), .groups = "drop") |>
         pivot_wider(names_from = method, values_from = w, values_fill = 0)))
w_dark <- 0.6 * max(wt$weight)   # white text on the darker tiles
p1 <- ggplot(wt, aes(factor(k), method, fill = weight)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", weight), color = weight > w_dark), size = 2.2) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = c("#e8f1fb", "#86b6ef", "#3987e5", "#184f95"), name = "Weight") +
  scale_y_discrete(labels = M_LAB, limits = rev) +
  facet_grid(scale ~ weighting) +
  labs(x = "k", y = NULL, title = "How much each input counts in the consensus",
       subtitle = "Normalised weights (sum to 1 within each panel). Blank cells: inputs left out.") +
  theme(panel.grid = element_blank())
save(p1, "cbma1_weights", 11, 9)

# ---- Consensus trajectories ---------------------------------------------------------------------
wide <- readRDS(file.path(DER, "wide.rds"))
tr <- bind_rows(lapply(runs, function(r) {
  R <- as.matrix(wide[match(r$ids, wide$src_subject_id), paste0("imp_r_", WAVES)])
  M <- rowsum(R, r$labels) / as.vector(table(r$labels))
  bind_cols(lab(r)[rep(1, nrow(M)), ], tibble(class = factor(as.integer(rownames(M))),
            share = as.vector(table(r$labels)) / r$n), as_tibble(M, .name_repair = ~ WAVES))
})) |> pivot_longer(all_of(WAVES), names_to = "wave", values_to = "cbcl") |> mutate(wave = factor(wave, WAVES))
p2 <- ggplot(tr, aes(wave, cbcl, group = class, color = class)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  geom_text(data = filter(tr, wave == WAVES[length(WAVES)]), aes(label = sprintf("%.0f%%", 100 * share)),
            hjust = -0.25, size = 2.4, color = INK2, show.legend = FALSE, check_overlap = TRUE) +
  facet_grid(scale ~ weighting) +
  scale_color_manual(values = ORD, name = "Consensus class (1 = lowest)") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
  labs(x = "Wave", y = "Mean raw CBCL total", title = "Consensus classes across methods and k",
       subtitle = "clusterBMA, K_BMA = 6; clusters nobody is allocated to drop out. Labels: class size.") +
  theme(legend.position = "top", legend.justification = "left")
save(p2, "cbma2_trajectories", 10, 8)

# ---- Uncertainty --------------------------------------------------------------------------------
un <- bind_rows(lapply(runs, \(r) bind_cols(lab(r)[rep(1, r$n), ], tibble(u = r$uncertainty))))
p3 <- ggplot(un, aes(u)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "#2a78d6", color = "white", linewidth = 0.2) +
  facet_grid(scale ~ weighting) +
  labs(x = "Allocation uncertainty (1 - probability of the assigned consensus class)", y = "Children",
       title = "How sure is the consensus about each child?",
       subtitle = "With hard inputs, uncertainty reflects disagreement between methods and k, not within-model uncertainty.")
save(p3, "cbma3_uncertainty", 10, 7)

# ---- Agreement of inputs with the consensus -------------------------------------------------------
say("\nARI of each input with its consensus, mean by method:")
say(pt(wt |> group_by(scale, weighting, method) |> summarise(ari = round(mean(ari), 2), .groups = "drop") |>
         pivot_wider(names_from = method, values_from = ari)))
p4 <- ggplot(wt, aes(ari, method, color = factor(k))) +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_grid(scale ~ weighting) +
  scale_color_manual(values = c("#b9d4f3", "#78acec", "#3987e5", "#1f5fae", "#0d3266"), name = "k") +
  scale_y_discrete(labels = M_LAB, limits = rev) +
  labs(x = "ARI with the consensus classes", y = NULL, title = "Which inputs does the consensus resemble?") +
  theme(legend.position = "top", legend.justification = "left")
save(p4, "cbma4_agreement", 10, 7)

# ---- Consensus vs consensus within a scale -----------------------------------------------------------
cc <- bind_rows(lapply(names(S_LAB), function(s) {
  rs <- runs[vapply(runs, \(r) r$scale == s, logical(1))]
  if (length(rs) < 2) return(NULL)
  pr <- combn(names(rs), 2)
  tibble(scale = S_LAB[[s]], a = W_LAB[vapply(pr[1, ], \(n) rs[[n]]$weighting, "")],
         b = W_LAB[vapply(pr[2, ], \(n) rs[[n]]$weighting, "")],
         ari = round(apply(pr, 2, \(p) adjustedRandIndex(rs[[p[1]]]$labels, rs[[p[2]]]$labels)), 2))
}))
if (nrow(cc)) { say("\nARI between consensus fits on the same scale (different weightings):"); say(pt(cc)) }

# ---- Poverty RRR per consensus fit --------------------------------------------------------------
pv <- ov |> filter(!is.na(rrr)) |>
  mutate(lo = as.numeric(sub("-.*", "", ci)), hi = as.numeric(sub(".*-", "", ci)),
         lab = sprintf("%.2f   highest %s%%, lowest %s%%", rrr, round(top_pct), sub("/.*", "", sizes)))
p5 <- ggplot(pv, aes(rrr, weighting)) +
  geom_vline(xintercept = 1, color = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", color = "#5b9be6", linewidth = 0.7) +
  { if (!all(is.na(pv$rrr_unadj))) geom_point(aes(x = rrr_unadj), shape = 21, size = 2.2, fill = "white", color = INK2, stroke = 0.6) } +
  geom_point(size = 2.4, color = "#0d3266") +
  geom_text(aes(x = max(hi) * 1.02, label = lab), hjust = 0, size = 2.5, color = INK2) +
  facet_grid(scale ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.6))) +
  scale_y_discrete(limits = rev) +
  labs(x = "Poverty RRR per halving of income, highest vs lowest consensus class (95% CI, log scale)", y = NULL,
       title = "Poverty gradient in every consensus fit",
       subtitle = paste("Filled: adjusted for age, sex, site and race/ethnicity, with 95% CI. Open: without race/ethnicity.",
                        "Text: adjusted RRR, then the share of children in the highest and lowest classes.")) +
  theme(strip.text.y = element_text(angle = 0))
save(p5, "cbma5_poverty", 9, 6.5)

# ---- Agreement between weightings ----------------------------------------------------------------
if (nrow(cc)) {
  cm <- bind_rows(cc, rename(cc, a = b, b = a),
                  distinct(bind_rows(select(cc, scale, a), select(cc, scale, a = b))) |> mutate(b = a, ari = 1)) |>
    filter(scale != S_LAB[["all"]]) |>
    mutate(scale = factor(scale, S_LAB), a = factor(a, W_LAB), b = factor(b, W_LAB))
  p6 <- ggplot(cm, aes(a, b, fill = ari)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", ari), color = ari > .75), size = 2.6) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
    scale_fill_gradientn(colors = c("#e8f1fb", "#86b6ef", "#3987e5", "#184f95"), limits = c(min(cm$ari), 1), name = "ARI") +
    facet_wrap(~ scale, nrow = 1) +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL, title = "Do the weightings give the same consensus?",
         subtitle = "ARI between consensus classes on the same scale. The pooled all-scales runs are left out (two fits).") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
  save(p6, "cbma6_weighting_agreement", 11, 4)
}

# ---- Uncertainty by consensus class -----------------------------------------------------------------
uc <- bind_rows(lapply(runs, \(r) bind_cols(lab(r)[rep(1, r$n), ], tibble(class = factor(r$labels), u = r$uncertainty)))) |>
  group_by(scale, weighting, class) |>
  summarise(mean_u = mean(u), over20 = mean(u > .2), n = n(), .groups = "drop") |>
  group_by(scale, weighting) |> mutate(share = n / sum(n)) |> ungroup()
say("\nUncertainty by consensus class (1 = lowest mean CBCL):")
say(pt(uc |> mutate(across(c(mean_u, over20, share), \(x) round(x, 3)))))
p7 <- ggplot(uc, aes(class, mean_u, fill = class)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * share)), vjust = -0.4, size = 2.3, color = INK2) +
  facet_grid(scale ~ weighting) +
  scale_fill_manual(values = ORD, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Consensus class (1 = lowest mean CBCL)", y = "Mean allocation uncertainty",
       title = "Which consensus classes are the methods unsure about?",
       subtitle = "Mean of 1 - probability of the assigned class, per class. Labels: class size.")
save(p7, "cbma7_uncertainty_by_class", 10, 7)
say("\nFigures written to ", figd)
