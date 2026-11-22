# =============================================================================
# R/cbma_check_figures.R (2026-11-22): figures for the two consensus robustness checks.
#   c1  poverty gradient, every class vs the lowest: primary vs four-method sensitivity set
#   c2  consensus class trajectories: primary vs sensitivity set
#   c3  poverty gradient by consensus method: clusterBMA vs evidence accumulation vs
#       weighted median partition (R/cbma_alt_consensus.R; run that first)
#   c4  agreement (ARI) between the three consensus methods, and between the two sets
# Race-adjusted multinomial as in 12b (age, sex, site, race/ethnicity). Aggregates only.
#   Rscript R/cbma_check_figures.R
# Output: output/figures/8_consensus_checks/
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages({ library(nnet); library(mclust) })

figd <- fig_dir("8_consensus_checks"); dir.create(figd, recursive = TRUE, showWarnings = FALSE)
wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()
S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)", all = "All scales pooled")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged")
fit_lab <- function(nm) { p <- strsplit(nm, "_")[[1]]; paste(S_LAB[p[1]], W_LAB[p[2]], sep = ", ") }
FIT_ORD <- unname(vapply(c("raw_ch", "raw_equal", "t_ch", "t_sdbw", "t_equal", "t_avg", "log_ch", "log_equal", "all_equal"), fit_lab, ""))

pov_fit <- function(ids, labels, f = POV_F$adj) {    # 12b's model; class 1 (lowest) is the reference
  w <- wide[match(ids, wide$src_subject_id), ]
  dat <- tibble(class = factor(labels), income_halving = -log2(w$income), age_bl = w$age_bl,
                sex = factor(w$sex), site = droplevels(factor(w$site)), race = w$race,
                cbcl_mean = rowMeans(as.matrix(w[paste0("imp_r_", WAVES)])))
  dat <- pov_rows(dat, f); dat$race <- droplevels(dat$race)
  m <- multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE)
  cf <- coef(m); V <- vcov(m)
  if (is.null(dim(cf))) cf <- matrix(cf, 1, dimnames = list(levels(dat$class)[2], names(cf)))
  cls <- rownames(cf)
  terms <- if (length(cls) == 1) "income_halving" else paste0(cls, ":income_halving")
  se <- sqrt(pmax(diag(V)[terms], 0))
  bind_rows(
    tibble(class = 1L, rrr = 1, lo = NA_real_, hi = NA_real_,
           class_mean = mean(dat$cbcl_mean[dat$class == levels(dat$class)[1]]),
           pct = 100 * mean(dat$class == levels(dat$class)[1])),
    tibble(class = as.integer(cls), rrr = exp(cf[, "income_halving"]),
           lo = exp(cf[, "income_halving"] - 1.96 * se), hi = exp(cf[, "income_halving"] + 1.96 * se),
           class_mean = as.vector(tapply(dat$cbcl_mean, dat$class, mean)[cls]),
           pct = as.vector(100 * table(dat$class)[cls] / nrow(dat))))
}

INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#ecebe8"
ORD <- c("#a6c8f0", "#5b9be6", "#2563b8", "#0d3266", "#08203f")
theme_set(theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = GRID, linewidth = 0.3),
        axis.text = element_text(color = INK2), axis.title = element_text(color = INK2),
        plot.title = element_text(color = INK, face = "bold"), plot.subtitle = element_text(color = INK2),
        strip.text = element_text(color = INK, face = "bold", hjust = 0),
        legend.position = "top", legend.justification = "left",
        plot.background = element_rect(fill = "white", color = NA)))
YB <- c(0.8, 0.9, 1, 1.1, 1.25, 1.5, 1.75)
gradient_plot <- function(d, col, title, subtitle, values) {
  ggplot(d, aes(class_mean, rrr, color = .data[[col]])) +
    geom_hline(yintercept = 1, color = INK2, linewidth = 0.3) +
    geom_line(linewidth = 0.6, position = position_dodge(width = 1.5)) +
    geom_pointrange(aes(ymin = lo, ymax = hi, shape = class == 1), size = 0.3, linewidth = 0.55,
                    position = position_dodge(width = 1.5)) +
    facet_wrap(~ fit, ncol = 3) +
    scale_color_manual(values = values, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
    scale_y_log10(breaks = YB) +
    labs(x = "Class mean raw CBCL (all waves)", y = "RRR per halving of income (log scale)", title = title, subtitle = subtitle)
}

# ---- c1, c2: primary vs sensitivity set ---------------------------------------------------
SET_LAB <- c(cbma = "Primary: k-means + kml", cbma_four = "Sensitivity: k-means, kml, PAM, GBTM")
fits <- sub("\\.rds$", "", list.files(file.path(OUT, "cbma"), "\\.rds$"))
runs <- list(); for (s in names(SET_LAB)) for (nm in fits) {
  f <- file.path(OUT, s, paste0(nm, ".rds")); if (file.exists(f)) runs[[paste(s, nm)]] <- readRDS(f)
}
key_set <- function(k) factor(SET_LAB[sub(" .*", "", k)], SET_LAB)
key_fit <- function(k) factor(fit_lab(sub(".* ", "", k)), FIT_ORD)

pv <- bind_rows(lapply(names(runs), \(k) mutate(pov_fit(runs[[k]]$ids, runs[[k]]$labels), set = key_set(k), fit = key_fit(k))))
p1 <- gradient_plot(pv, "set", "Poverty gradient, every class vs the lowest: primary vs sensitivity input set",
  "clusterBMA, k = 2-8, K_BMA = 8. Race-adjusted multinomial, 95% CI. Open point: the lowest class (reference).",
  c("#2563b8", "#b0aea8"))
ggsave(file.path(figd, "c1_sets_poverty.png"), p1, width = 12, height = 11, dpi = 300)

tr <- bind_rows(lapply(names(runs), function(k) {
  r <- runs[[k]]
  R <- as.matrix(wide[match(r$ids, wide$src_subject_id), paste0("imp_r_", WAVES)])
  M <- rowsum(R, r$labels) / as.vector(table(r$labels))
  as_tibble(M, .name_repair = ~ WAVES) |>
    mutate(class = factor(seq_len(nrow(M))), share = as.vector(table(r$labels)) / r$n, set = key_set(k), fit = key_fit(k))
})) |> pivot_longer(all_of(WAVES), names_to = "wave", values_to = "cbcl") |> mutate(wave = factor(wave, WAVES))
p2 <- ggplot(tr, aes(wave, cbcl, group = class, color = class)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  geom_text(data = filter(tr, wave == WAVES[length(WAVES)]), aes(label = sprintf("%.0f%%", 100 * share)),
            hjust = -0.3, size = 2.6, color = INK2, show.legend = FALSE) +
  facet_grid(fit ~ set) + scale_color_manual(values = ORD, name = "Consensus class (1 = lowest)") +
  scale_x_discrete(expand = expansion(add = c(0.3, 0.9))) +
  labs(x = "Wave", y = "Mean raw CBCL total", title = "Consensus classes: primary vs sensitivity input set",
       subtitle = "clusterBMA, K_BMA = 8; unused clusters drop out. Labels: share of children.")
ggsave(file.path(figd, "c2_sets_trajectories.png"), p2, width = 9, height = 3 + 2.4 * length(fits), dpi = 300, limitsize = FALSE)

# ---- c3, c4: consensus method -------------------------------------------------------------
alt <- readRDS(file.path(RES, "cbma_alt_consensus.rds"))
M_LAB <- c(bma = "clusterBMA", ea = "Evidence accumulation", median = "Median partition")
# The stored models give classes 2..K against class 1; the reference itself is the line at 1.
pm <- bind_rows(lapply(alt, \(o) mutate(o$poverty, fit = factor(fit_lab(o$fit), FIT_ORD)))) |>
  mutate(consensus = factor(M_LAB[consensus], M_LAB)) |>
  filter(pct >= 1)   # evidence accumulation can peel off a handful of children (3 in all_equal)
p3 <- ggplot(pm, aes(class_mean, rrr, color = consensus)) +
  geom_hline(yintercept = 1, color = INK2, linewidth = 0.3) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 1.5)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3, linewidth = 0.55, position = position_dodge(width = 1.5)) +
  facet_wrap(~ fit, ncol = 3) +
  scale_color_manual(values = c("#2563b8", "#c0611f", "#7a7975"), name = NULL) +
  scale_y_log10(breaks = YB) +
  labs(x = "Class mean raw CBCL (all waves)", y = "RRR per halving of income vs the lowest class (log scale)",
       title = "Poverty gradient by consensus method",
       subtitle = paste("Same inputs, weights and number of classes; only the step that turns the consensus into classes differs.",
                        "Race-adjusted, 95% CI; the lowest class is the reference (RRR = 1). Classes under 1% of children not shown.", sep = "\n"))
ggsave(file.path(figd, "c3_methods_poverty.png"), p3, width = 12, height = 10, dpi = 300)

ag <- bind_rows(
  bind_rows(lapply(alt, \(o) tibble(fit = fit_lab(o$fit), pair = c("clusterBMA vs\nevidence accumulation", "clusterBMA vs\nmedian partition", "Evidence accumulation\nvs median partition"),
                                    ari = unname(o$ari[c("bma_ea", "bma_median", "ea_median")])))),
  bind_rows(lapply(fits, function(nm) {
    a <- runs[[paste("cbma", nm)]]; b <- runs[[paste("cbma_four", nm)]]
    if (is.null(b)) return(NULL)
    tibble(fit = fit_lab(nm), pair = "Primary vs\nsensitivity set", ari = adjustedRandIndex(a$labels, b$labels[match(a$ids, b$ids)]))
  }))) |> mutate(fit = factor(fit, rev(FIT_ORD)), pair = factor(pair, unique(pair)))
p4 <- ggplot(ag, aes(pair, fit, fill = ari)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", ari), color = ari > 0.7), size = 3.3, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient(low = "#e8f0fb", high = "#0d3266", limits = c(0, 1), name = "ARI") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "How much do the consensus classes depend on the method and the input set?",
       subtitle = "Adjusted Rand index between consensus partitions of the same 8,331 children.") +
  theme(panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(figd, "c4_agreement.png"), p4, width = 9, height = 5.5, dpi = 300)
cat("Figures written to", figd, "\n")
