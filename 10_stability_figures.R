# =============================================================================
# 10_stability_figures.R: stability figures (see docs/STABILITY_PLAN.md).
#   Rscript 10_stability_figures.R full|smoke
# Main:
#   stab1_heatmap      mean ARI by method x k; facets perturbation x anchor. Seeds
#                      for GBTM/GMM = converged recorded starts vs the best start;
#                      Ward is deterministic (no seeds)
#   stab2_children     child-level consistency: share of children reliably in the
#                      highest class (>= .90 of replicates) vs on the boundary
#   stab3_null         real vs null: silhouette and entropy by k; k selected
#   stab4_poverty      poverty RRR across subsamples / site halves / noise, the
#                      primary-run (full_refit) estimate with its CI, and the null distribution
# Supplement:
#   stabS1_starts      GBTM/GMM: log-likelihood of each converged start vs the best
#   stabS2_settings    settings variants: ARI with base; subsample stability
# =============================================================================

source("00_setup.R")
source("R/perturb.R")
suppressPackageStartupMessages(library(patchwork))

mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "full"
SOUT <- file.path(OUT, if (mode == "smoke") "stability_smoke" else "stability")
S    <- readRDS(file.path(SOUT, "summary.rds"))
figd <- fig_dir(if (mode == "smoke") "stability_smoke" else "stability")
dir.create(figd, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm", "gbtm_link")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM", gbtm_link = "GBTM (link)")
A_LAB <- c(raw = "Raw CBCL", log = "log(1+raw)")
P_LAB <- c(seed = "Seeds / starts", subsample = "80% subsamples", site = "Site split-half",
           noise90 = "Noise, r = .90", noise80 = "Noise, r = .80")
BLUE <- "#2a78d6"; GRAY <- "#9a9893"; INK <- "#0b0b0b"; INK2 <- "#52514e"; RED <- "#e34948"
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
fm <- function(x) factor(x, METHODS)
fa <- function(x) factor(A_LAB[x], A_LAB)

# ---- Stab-1: stability heatmap ------------------------------------------------------------
h_pairs <- S$pairs |> filter(variant == "base", j > 0, ptype %in% names(P_LAB)) |>
  group_by(anchor, method, ptype, k) |> summarise(ari = mean(ari, na.rm = TRUE), .groups = "drop")
h_site <- S$site |> group_by(anchor, method, k) |>
  summarise(ari = mean(ari, na.rm = TRUE), .groups = "drop") |> mutate(ptype = "site")
h_starts <- if (nrow(S$starts)) S$starts |> filter(!is.na(converged_run)) |>
  group_by(anchor, method, k) |> summarise(ari = mean(ari_full_vs_best, na.rm = TRUE), .groups = "drop") |>
  mutate(ptype = "seed")
h <- bind_rows(h_pairs, h_site, h_starts) |>
  complete(anchor, method = METHODS, ptype = names(P_LAB), k = K_STAB) |>
  mutate(method = fm(method), anchor = fa(anchor), ptype = factor(P_LAB[ptype], P_LAB),
         lab = case_when(!is.na(ari) ~ sprintf("%.2f", ari),
                         method == "ward" & ptype == P_LAB[["seed"]] ~ "det.", TRUE ~ ""))
p1 <- ggplot(h, aes(factor(k), method, fill = ari)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = lab, color = !is.na(ari) & ari > 0.55), size = 2.4) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = SEQ, limits = c(0, 1), na.value = "#f6f5f3", name = "Mean ARI") +
  scale_y_discrete(labels = M_LAB, limits = rev) +
  facet_grid(anchor ~ ptype) +
  labs(x = "k", y = NULL, title = "Do the same classes come back?",
       subtitle = paste("Mean adjusted Rand index between replicate solutions (subsamples: on shared children; site split: training half's",
                        "classes predicted\nfor the test half vs its own fit). GBTM/GMM seeds: converged recorded starts vs the best start.",
                        "det. = deterministic (Ward), no seed variation.")) +
  theme(panel.grid = element_blank())
save(p1, "stab1_heatmap", 13, 6)

# ---- Stab-2: child-level consistency -------------------------------------------------------
cc <- S$child |> filter(ptype == "pooled") |>
  select(anchor, method, k, `Reliably high (>= .90)` = reliably_high, `Boundary (.10-.90)` = boundary) |>
  pivot_longer(-c(anchor, method, k), names_to = "group", values_to = "share") |>
  mutate(method = fm(method), anchor = fa(anchor),
         group = factor(group, c("Boundary (.10-.90)", "Reliably high (>= .90)")))
p2 <- ggplot(cc, aes(factor(k), share, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`Reliably high (>= .90)` = SEQ[6], `Boundary (.10-.90)` = SEQ[3]), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
  labs(x = "k", y = "Share of children",
       title = "Who is reliably in the highest-symptom class?",
       subtitle = paste("Per child: share of replicates (seeds, subsamples, site halves, noise) that include the child and place them",
                        "in the highest class.\nChildren never high (<= .10) are the remainder.")) +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.x = element_blank())
save(p2, "stab2_children", 13, 5.5)

# ---- Stab-3: null benchmark ------------------------------------------------------------------
ix <- bind_rows(
  S$base |> transmute(data = "Real (primary run)", anchor, method, k, Silhouette = asw, `Relative entropy` = entropy),
  S$fits |> filter(ptype == "null", is.na(error)) |>
    transmute(data = "Null (no subgroups)", anchor, method, k, Silhouette = asw, `Relative entropy` = entropy)) |>
  pivot_longer(c(Silhouette, `Relative entropy`), names_to = "index", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(method = fm(method), anchor = fa(anchor))
p3a <- ggplot(ix, aes(k, value, color = data)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
  stat_summary(fun.data = \(y) data.frame(y = mean(y), ymin = min(y), ymax = max(y)), geom = "pointrange", size = 0.15) +
  scale_color_manual(values = c(`Real (primary run)` = BLUE, `Null (no subgroups)` = GRAY), name = NULL) +
  facet_grid(index + anchor ~ method, scales = "free_y", labeller = labeller(method = M_LAB)) +
  labs(x = "k", y = NULL, title = "Would the methods find classes in data with no subgroups?",
       subtitle = "Null: CBCL simulated from a one-class growth model with the real children's ages. Points: mean and range over null datasets.") +
  theme(legend.position = "top", legend.justification = "left")
ks <- bind_rows(
  S$base |> filter(k_found == k) |> group_by(anchor, method) |>
    summarise(k_sel = k[which.max(asw)], .groups = "drop") |> mutate(data = "Real (primary run)"),
  S$ksel |> filter(ptype == "null") |> transmute(anchor, method, k_sel = k_asw, data = "Null (no subgroups)")) |>
  count(data, anchor, method, k_sel) |> group_by(data, anchor, method) |> mutate(share = n / sum(n)) |> ungroup() |>
  mutate(method = fm(method), anchor = fa(anchor))
p3b <- ggplot(ks, aes(factor(k_sel), share, fill = data)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = c(`Real (primary run)` = BLUE, `Null (no subgroups)` = GRAY), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
  labs(x = "k selected by silhouette", y = "Share", subtitle = "k chosen by maximum silhouette") +
  theme(legend.position = "none", panel.grid.major.x = element_blank())
save(p3a / p3b + plot_layout(heights = c(2.2, 1)), "stab3_null", 13, 11)

# ---- Stab-4: poverty RRR -----------------------------------------------------------------------
pv <- S$poverty |> filter(variant == "base", ptype %in% c("subsample", "site", "noise90", "noise80", "null")) |>
  mutate(ptype = factor(c(P_LAB, null = "Null (no subgroups)")[ptype],
                        c(P_LAB[c("subsample", "site", "noise90", "noise80")], "Null (no subgroups)")),
         method = fm(method), anchor = fa(anchor))
bref <- S$base |> filter(!is.na(rrr)) |> mutate(method = fm(method), anchor = fa(anchor))
p4 <- ggplot(pv, aes(factor(k), rrr)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.4) +
  geom_boxplot(aes(fill = ptype), outlier.size = 0.3, linewidth = 0.25, width = 0.75,
               position = position_dodge(0.85)) +
  geom_pointrange(data = bref, aes(x = factor(k), y = rrr, ymin = lo, ymax = hi), inherit.aes = FALSE,
                  color = RED, size = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = c(SEQ[c(2, 4, 5, 7)], GRAY), name = NULL) +
  scale_y_log10() +
  facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
  labs(x = "k", y = "RRR per halving of income\n(highest vs lowest class, log scale)",
       title = "How much does the poverty estimate move when the data move?",
       subtitle = "Boxes: replicate RRRs. Red: primary-run estimate with 95% CI. Gray: null datasets, where income is unrelated to CBCL by construction.") +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.x = element_blank())
save(p4, "stab4_poverty", 14, 6.5)

# ---- Stab-S1: start recording ---------------------------------------------------------------------
if (nrow(S$starts)) {
  st <- S$starts |> filter(!is.na(converged_run)) |>
    mutate(gap = best_ll - ll_full, method = fm(method), anchor = fa(anchor))
  pS1 <- ggplot(st, aes(factor(k), pmax(gap, 1e-3))) +
    geom_hline(yintercept = 0.01, color = RED, linewidth = 0.4, linetype = "22") +
    geom_jitter(aes(color = conv_full %in% 1), width = 0.15, height = 0, size = 1.4) +
    scale_color_manual(values = c(`TRUE` = BLUE, `FALSE` = GRAY), labels = c(`TRUE` = "converged", `FALSE` = "not converged"), name = NULL) +
    scale_y_log10() +
    facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
    labs(x = "k", y = "Best log-likelihood minus this run's\n(log scale; floored at .001)",
         title = "Is the best log-likelihood replicated?",
         subtitle = sprintf("The top %d of %d random starts (by 10-iteration log-likelihood),\neach run to convergence. Dashed line: within 0.01 of the best.",
                            N_CONVERGE, N_STARTS)) +
    theme(legend.position = "top", legend.justification = "left")
  save(pS1, "stabS1_starts", 8, 6)
}

# ---- Stab-S2: settings variants ---------------------------------------------------------------------
sv <- S$settings |> mutate(method = fm(method), anchor = fa(anchor), lab = paste(M_LAB[as.character(method)], variant))
pS2a <- ggplot(filter(sv, variant != "base"), aes(factor(k), ari_vs_base, color = variant, group = variant)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.2) +
  scale_y_continuous(limits = c(NA, 1)) +
  facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
  labs(x = "k", y = "ARI with base settings", color = "Variant",
       title = "Do method settings change the solution?",
       subtitle = "Full anchor data: each variant's solution vs the base-settings solution")
pS2b <- ggplot(filter(sv, !is.na(sub_ari_mean)), aes(factor(k), sub_ari_mean, color = variant, group = variant)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.2) +
  facet_grid(anchor ~ method, labeller = labeller(method = M_LAB)) +
  labs(x = "k", y = "Mean ARI between\nsubsample fits", color = "Settings",
       subtitle = "Subsample stability (the same 20 subsamples) for base and variant settings, cheap methods")
save(pS2a / pS2b, "stabS2_settings", 11, 9)

cat("Figures written to", figd, "\n")
