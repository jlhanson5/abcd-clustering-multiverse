# =============================================================================
# 14_extra_figures.R: figures that put the single-method solutions (full_refit,
# anchor configuration: wave, one per family, >=3 waves imputed, k 2-6) next to
# the clusterBMA consensus fits (12_cbma.R).
#   Rscript 14_extra_figures.R
# Figures in output/figures/5_single_vs_consensus/:
#   x1_rrr_inputs_vs_consensus  poverty RRR of every input solution and every consensus fit
#   x2_class_map                size and mean CBCL of every class; highest classes marked
#   x3_crosswalk                consensus (equal weights) classes vs each method's k = 3 classes
#   x4_child_uncertainty        per-child consensus uncertainty against the child's mean CBCL
#   x5_class_gradient           RRR of every class vs the lowest, against class mean CBCL
#   x6_highest_share            share of input solutions that put each child in the highest class
#   x7_consensus_agreement      ARI between all 16 consensus fits, across scales
# Table: output/summaries/summary_extra_figures.txt
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(mclust) })

figd <- fig_dir("extra")
sumf <- file.path(SUMM, "summary_extra_figures.txt"); unlink(sumf)
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = sumf, append = TRUE) }
pt <- function(x) paste(capture.output(print(x, n = 200, width = 200)), collapse = "\n")

METHODS <- c("kmeans", "pam", "ward", "llpa", "kml", "gbtm", "gmm", "gbtm_link")
M_LAB <- c(kmeans = "k-means", pam = "PAM", ward = "Ward", llpa = "Gaussian mixture",
           kml = "kml", gbtm = "GBTM", gmm = "GMM", gbtm_link = "GBTM (link)")
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)", all = "All scales")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged", kbal = "k-balanced")
# Palette as in 11_refit_figures.R / 13_cbma_summary.R
SIG_COL <- c("CI above 1" = "#2a78d6", "CI includes 1" = "#9a9893", "CI below 1" = "#e34948")
ORD <- c("#a6c8f0", "#5b9be6", "#2563b8", "#0d3266", "#08203f", "#04101f")   # class / k order, light to dark
K_COL <- c("#b9d4f3", "#78acec", "#3987e5", "#1f5fae", "#0d3266")
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"; GRAY <- "#b8b6b0"
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        legend.text = element_text(color = INK2), legend.title = element_text(color = INK2),
        plot.background = element_rect(fill = "white", color = NA)))
save <- function(p, name, w, h) ggsave(file.path(figd, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
sig_of <- function(lo, hi) factor(case_when(lo > 1 ~ "CI above 1", hi < 1 ~ "CI below 1", TRUE ~ "CI includes 1"),
                                  names(SIG_COL))

# ---- Data -------------------------------------------------------------------------------------
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race() |>
  mutate(income_halving = -log2(income), cbcl_mean = rowMeans(across(all_of(paste0("imp_r_", WAVES)))),
         sex = factor(sex), site = factor(site))
res <- readRDS(file.path(RES, "results_full_refit.rds"))
files <- list.files(file.path(OUT, "cbma"), "[.]rds$", full.names = TRUE)
runs <- setNames(lapply(files, readRDS), sub("[.]rds$", "", basename(files)))
runs <- runs[order(match(vapply(runs, `[[`, "", "scale"), names(S_LAB)),
                   match(vapply(runs, `[[`, "", "weighting"), names(W_LAB)))]
lab <- function(r) tibble(scale = factor(S_LAB[r$scale], S_LAB), weighting = factor(W_LAB[r$weighting], W_LAB))
ids <- runs[[1]]$ids
stopifnot(all(vapply(runs, \(r) identical(r$ids, ids), logical(1))))

# Input solutions: the 32 per scale used by clusterBMA (the "all" run pools them)
inp <- bind_rows(lapply(runs[vapply(runs, \(r) r$scale != "all", logical(1))], \(r) r$inputs)) |>
  distinct(universe_id, input_scale, method, k) |>
  left_join(select(res$fits, universe_id, rrr, lo, hi, top_n), by = "universe_id") |>
  mutate(scale = factor(S_LAB[input_scale], S_LAB), method = factor(method, METHODS))
say(sprintf("Extra figures (%s): %d input solutions, %d consensus fits, n = %d.",
            format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(inp), length(runs), length(ids)))
lab_of <- function(uid) res$assign[[uid]]$class[match(ids, res$assign[[uid]]$src_subject_id)]
L <- sapply(inp$universe_id, lab_of)              # n x inputs, classes relabelled 1..K by mean CBCL
stopifnot(!anyNA(L))
cm <- wide$cbcl_mean[match(ids, wide$src_subject_id)]

# ---- Poverty model on every consensus class (as in 04_poverty.R / 12b_cbma_poverty.R) ----------------
# Race-adjusted (primary since 2026-11-09); input RRRs from results_full_refit.rds are adjusted too.
pov_all <- function(labels, f = POV_F$adj) {
  dat <- tibble(src_subject_id = ids, class = factor(labels)) |> inner_join(wide, by = "src_subject_id") |>
    pov_rows(f) |> mutate(site = droplevels(site), race = droplevels(race))
  m <- multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE)
  cf <- coef(m); V <- vcov(m)
  if (is.null(dim(cf))) cf <- matrix(cf, 1, dimnames = list(levels(dat$class)[2], names(cf)))
  cls <- rownames(cf)
  se <- sqrt(diag(V)[paste0(cls, ":income_halving")])
  means <- tapply(dat$cbcl_mean, dat$class, mean); ns <- table(dat$class)
  tibble(class = as.integer(cls), est = cf[, "income_halving"], se = se,
         class_mean = as.vector(means[cls]), n = as.vector(ns[cls])) |>
    mutate(rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}
cpov <- bind_rows(lapply(runs, \(r) bind_cols(lab(r)[1, ], pov_all(r$labels))))
cp12b <- readRDS(file.path(RES, "cbma_poverty.rds")) |> filter(model == "adj")
chk <- cpov |> transmute(sc = names(S_LAB)[as.integer(scale)], wt = names(W_LAB)[as.integer(weighting)], class, rrr) |>
  left_join(cp12b, by = c(sc = "scale", wt = "weighting", "class"))
stopifnot(nrow(chk) == nrow(cp12b), max(abs(chk$rrr.x - chk$rrr.y)) < 1e-6)   # same model as 12b_cbma_poverty.R
chk_u <- vapply(runs, \(r) { p <- pov_all(r$labels, POV_F$unadj); p$rrr[p$class == max(p$class)] - r$poverty$rrr }, 1)
stopifnot(max(abs(chk_u)) < 1e-6)                                            # unadjusted = 12_cbma.R's stored model

# ---- x1: poverty RRR, inputs vs consensus ------------------------------------------------------------
top <- cpov |> group_by(scale, weighting) |> slice_max(class, n = 1) |> ungroup()
d1 <- bind_rows(
  inp |> transmute(scale, row = M_LAB[as.character(method)], k = factor(k), rrr, lo, hi, what = "input"),
  top |> filter(scale != "All scales") |> transmute(scale, row = "Consensus", k = NA, rrr, lo, hi, what = as.character(weighting))) |>
  mutate(row = factor(row, rev(c(M_LAB[METHODS], "Consensus"))))
all_eq <- top$rrr[top$scale == "All scales"]
say("\nPoverty RRR, input solutions by scale (median, range) and consensus fits:")
say(pt(inp |> group_by(scale) |> summarise(n = sum(!is.na(rrr)), median = round(median(rrr, na.rm = TRUE), 2),
                                            min = round(min(rrr, na.rm = TRUE), 2), max = round(max(rrr, na.rm = TRUE), 2))))
say(pt(top |> transmute(scale, weighting, rrr = round(rrr, 2), ci = sprintf("%.2f-%.2f", lo, hi), top_pct = round(100 * n / length(ids), 1))))
p1 <- ggplot(d1, aes(rrr, row)) +
  geom_vline(xintercept = 1, color = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = all_eq, color = "#2563b8", linewidth = 0.4, linetype = "dotted") +
  geom_point(data = filter(d1, what == "input"), aes(color = k), size = 2, alpha = 0.9,
             position = position_jitter(height = 0.12, width = 0, seed = 1)) +
  geom_point(data = filter(d1, what != "input"), shape = 23, size = 2.6, fill = "#0d3266", color = "white", stroke = 0.4) +
  facet_wrap(~ scale, ncol = 1) +
  scale_color_manual(values = K_COL, name = "Input k", na.translate = FALSE) +
  scale_x_log10() +
  labs(x = "Poverty RRR per halving of income, highest vs lowest class (log scale)", y = NULL,
       title = "Single-method solutions and the consensus",
       subtitle = "Dots: the inputs clusterBMA combines (full_refit, anchor configuration). Diamonds: the five consensus weightings.\nDotted line: the pooled all-scales consensus (equal weights). Adjusted for age, sex, site and race/ethnicity.") +
  theme(legend.position = "top", legend.justification = "left")
save(p1, "x1_rrr_inputs_vs_consensus", 8, 9)

# ---- x2: class map ---------------------------------------------------------------------------------
d2_in <- bind_rows(lapply(seq_len(nrow(inp)), function(j) {
  tibble(scale = inp$scale[j], method = inp$method[j], k = inp$k[j], class = L[, j], cm = cm) |>
    group_by(scale, method, k, class) |> summarise(share = n() / length(ids), mean = mean(cm), .groups = "drop") |>
    mutate(highest = class == max(class))
}))
d2_cb <- bind_rows(lapply(runs, \(r) bind_cols(lab(r)[1, ], tibble(class = r$labels, cm = cm) |>
  group_by(class) |> summarise(share = n() / length(ids), mean = mean(cm)) |> mutate(highest = class == max(class)))))
d2_cb <- bind_rows(filter(d2_cb, scale != "All scales"),
                   bind_rows(lapply(S_LAB[1:3], \(s) filter(d2_cb, scale == "All scales") |> mutate(scale = s)))) |>
  mutate(scale = factor(scale, S_LAB[1:3]))
p2 <- ggplot(mapping = aes(share, mean)) +
  geom_point(data = d2_in, aes(color = highest), size = 1.4, alpha = 0.7) +
  geom_point(data = filter(d2_cb, highest), shape = 23, size = 2.6, fill = "#0d3266", color = "white", stroke = 0.4) +
  geom_point(data = filter(d2_cb, !highest), shape = 23, size = 2.2, fill = "#9a9893", color = "white", stroke = 0.4) +
  facet_wrap(~ scale, nrow = 1) +
  scale_color_manual(values = c(`TRUE` = "#3987e5", `FALSE` = GRAY), labels = c(`TRUE` = "highest class", `FALSE` = "other classes"), name = "Input solutions") +
  scale_x_log10(labels = scales::percent) +
  labs(x = "Class size (share of children, log scale)", y = "Class mean raw CBCL (BL-Y3)",
       title = "What the \"highest\" class is on each scale",
       subtitle = "Every class of every input solution (k 2-6). Diamonds: consensus classes (dark = highest), all weightings including the pooled run.") +
  theme(legend.position = "top", legend.justification = "left")
save(p2, "x2_class_map", 10, 4.5)

# ---- x3: crosswalk, consensus (equal) vs each method's k = 3 ------------------------------------------
d3 <- bind_rows(lapply(names(runs)[vapply(runs, \(r) r$weighting == "equal" && r$scale != "all", logical(1))], function(nm) {
  r <- runs[[nm]]
  js <- which(inp$input_scale == r$scale & inp$k == 3)
  bind_rows(lapply(js, \(j) as_tibble(as.data.frame(table(consensus = r$labels, method_class = L[, j]))) |>
    group_by(consensus) |> mutate(p = Freq / sum(Freq)) |> ungroup() |>
    mutate(scale = factor(S_LAB[r$scale], S_LAB), method = factor(M_LAB[as.character(inp$method[j])], M_LAB))))
}))
p3 <- ggplot(d3, aes(method_class, consensus, fill = p)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = round(100 * p), color = p > .6), size = 2.3) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = c("#f4f8fd", "#86b6ef", "#3987e5", "#184f95"), limits = c(0, 1), labels = scales::percent,
                       name = "Share of the\nconsensus class") +
  facet_grid(scale ~ method, scales = "free_y", space = "free_y") +
  labs(x = "Method's class at k = 3 (1 = lowest mean CBCL)", y = "Consensus class (equal weights)",
       title = "Where each consensus class comes from",
       subtitle = "Row percentages: of the children in a consensus class, the % in each class of the method's k = 3 solution.") +
  theme(panel.grid = element_blank())
save(p3, "x3_crosswalk", 11, 5.5)

# ---- x4: child uncertainty vs mean CBCL ---------------------------------------------------------------
d4 <- bind_rows(lapply(runs[vapply(runs, \(r) r$weighting == "equal", logical(1))], \(r)
  bind_cols(lab(r)[rep(1, length(ids)), ], tibble(cm = cm, u = r$uncertainty, class = factor(r$labels)))))
p4 <- ggplot(d4, aes(cm, u, color = class)) +
  geom_point(size = 0.35, alpha = 0.35) +
  facet_wrap(~ scale, nrow = 1) +
  scale_color_manual(values = ORD, name = "Consensus class") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 5, 10, 20, 50, 100)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "Child's mean raw CBCL, BL-Y3 (log(1+x) axis)", y = "Allocation uncertainty",
       title = "Which children the methods disagree about",
       subtitle = "Equal-weight consensus fits. Uncertainty peaks at the borders between classes.") +
  theme(legend.position = "top", legend.justification = "left")
save(p4, "x4_child_uncertainty", 11, 4)

# ---- x5: every class's RRR vs class mean ----------------------------------------------------------------
d5_in <- res$contrasts |> semi_join(inp, by = "universe_id") |>
  left_join(select(inp, universe_id, scale), by = "universe_id")
d5_cb <- cpov |> filter(scale != "All scales") |> mutate(scale = factor(scale, S_LAB[1:3]))
say("\nConsensus contrasts, every class vs the lowest:")
say(pt(cpov |> transmute(scale, weighting, class, n, class_mean = round(class_mean, 1), rrr = round(rrr, 2),
                         ci = sprintf("%.2f-%.2f", lo, hi))))
p5 <- ggplot(mapping = aes(class_mean, rrr)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_point(data = d5_in, color = GRAY, size = 1.2, alpha = 0.7) +
  geom_linerange(data = d5_cb, aes(ymin = lo, ymax = hi, color = weighting), linewidth = 0.5) +
  geom_point(data = d5_cb, aes(color = weighting), size = 1.8, alpha = 0.9) +
  facet_wrap(~ scale, nrow = 1, scales = "free_x") +
  scale_color_manual(values = setNames(c("#a6c8f0", "#5b9be6", "#2563b8", "#0d3266", "#e8a33a"), W_LAB), name = NULL) +
  scale_y_log10() +
  labs(x = "Class mean raw CBCL (BL-Y3)", y = "Poverty RRR vs the lowest class (log scale)",
       title = "The poverty gradient rises with class severity",
       subtitle = "Grey: every non-reference class of the input solutions. Coloured: consensus classes, 95% CI.\nWhen the lowest class is a small near-zero group (log scale, equal and k-balanced), every contrast shifts.") +
  theme(legend.position = "top", legend.justification = "left")
save(p5, "x5_class_gradient", 11, 4.5)

# ---- x6: how often is a child in the highest class ------------------------------------------------------
hi_mat <- sweep(L, 2, apply(L, 2, max), `==`)
d6 <- bind_rows(lapply(runs[vapply(runs, \(r) r$weighting == "equal", logical(1))], function(r) {
  w <- r$inputs$weight; H <- hi_mat[, r$inputs$universe_id, drop = FALSE]
  bind_cols(lab(r)[rep(1, length(ids)), ], tibble(share = as.vector(H %*% w) / sum(w), class = factor(r$labels)))
}))
say("\nChildren put in the highest class by at least half / all of the inputs (equal weights):")
say(pt(d6 |> group_by(scale) |> summarise(half = round(mean(share >= .5), 3), all = round(mean(share > .999), 3),
                                          never = round(mean(share == 0), 3))))
p6 <- ggplot(filter(d6, share > 0), aes(share, fill = class)) +
  geom_histogram(binwidth = 1 / 32, boundary = 0, color = "white", linewidth = 0.15) +
  facet_wrap(~ scale, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = ORD, name = "Consensus class") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Share of input solutions that put the child in their highest class", y = "Children",
       title = "Is a child \"severe\" in most solutions or only a few?",
       subtitle = "Equal-weight consensus fits; children never in a highest class are left out (counts in the summary table).") +
  theme(legend.position = "top", legend.justification = "left")
save(p6, "x6_highest_share", 11, 4)

# ---- x7: agreement between all consensus fits ------------------------------------------------------------
nm <- vapply(runs, \(r) sprintf("%s · %s", S_LAB[[r$scale]], W_LAB[[r$weighting]]), "")
A <- outer(seq_along(runs), seq_along(runs), Vectorize(\(i, j) adjustedRandIndex(runs[[i]]$labels, runs[[j]]$labels)))
d7 <- expand_grid(a = factor(nm, nm), b = factor(nm, nm)) |> mutate(ari = as.vector(A))
p7 <- ggplot(d7, aes(a, b, fill = ari)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", ari), color = ari > .7), size = 2) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colors = c("#f4f8fd", "#86b6ef", "#3987e5", "#184f95"), limits = c(0, 1), name = "ARI") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, title = "Do the consensus fits agree across scales?",
       subtitle = "ARI between the class assignments of all 16 clusterBMA fits.") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 50, hjust = 1))
save(p7, "x7_consensus_agreement", 9, 8)

say("\nFigures written to ", figd)
