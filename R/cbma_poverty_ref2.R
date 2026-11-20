# =============================================================================
# R/cbma_poverty_ref2.R: supplementary table (decided 2026-11-20, docs/DECISIONS_memory.md).
# On the log and T-score scales every clusterBMA consensus fit splits off a near-zero
# lowest class, and the next class up is less poor than it (a U-shaped gradient). The
# main text keeps every class vs the lowest (12b_cbma_poverty.R). This table refits the
# same multinomial models on the stored consensus classes with class 2 as the reference,
# for the log and T fits only. Race-adjusted (primary) and without race. No reclustering.
#   Rscript R/cbma_poverty_ref2.R              # primary set (output/cbma)
#   CBMA_SET=four Rscript R/cbma_poverty_ref2.R  # sensitivity set (output/cbma_four)
# Output: output/results/cbma_poverty_ref2<sfx>.rds (one row per fit x model x class other
#   than 2), output/summaries/summary_cbma_poverty_ref2<sfx>.txt. Aggregates only.
# =============================================================================

source("00_setup.R")
source("R/cbma_set.R")
suppressPackageStartupMessages(library(nnet))

wide <- readRDS(file.path(DER, "wide.rds")) |> add_race()
files <- cbma_files()
files <- files[sub("_.*", "", basename(files)) %in% c("log", "t")]

pov_ref2 <- function(r, model) {
  f <- POV_F[[model]]
  w <- wide[match(r$ids, wide$src_subject_id), ]
  dat <- tibble(class = factor(r$labels), income_halving = -log2(w$income), age_bl = w$age_bl,
                sex = factor(w$sex), site = droplevels(factor(w$site)), race = w$race,
                cbcl_mean = rowMeans(as.matrix(w[paste0("imp_r_", WAVES)])))
  dat <- pov_rows(dat, f); dat$race <- droplevels(dat$race)
  stopifnot(nlevels(dat$class) >= 3)
  dat$class <- relevel(dat$class, ref = "2")
  m <- multinom(f, data = dat, trace = FALSE, maxit = 1000, MaxNWts = 5000, Hess = TRUE)
  cf <- coef(m); V <- vcov(m)
  cls <- rownames(cf)
  terms <- paste0(cls, ":income_halving")
  ref_mean <- mean(dat$cbcl_mean[dat$class == "2"])
  tibble(scale = r$scale, weighting = r$weighting, model = model, n_model = nrow(dat),
         K = nlevels(dat$class), class = as.integer(cls), est = cf[, "income_halving"],
         se = sqrt(pmax(diag(V)[terms], 0)),
         n = as.vector(table(dat$class)[cls]), class_mean = as.vector(tapply(dat$cbcl_mean, dat$class, mean)[cls]),
         ref_n = sum(dat$class == "2"), ref_mean = ref_mean, converged = m$convergence == 0) |>
    mutate(rrr = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

out <- bind_rows(lapply(files, function(fn) {
  r <- readRDS(fn)
  bind_rows(pov_ref2(r, "adj"), pov_ref2(r, "unadj"))
}))
stopifnot(all(out$converged))
saveRDS(out, file.path(RES, paste0("cbma_poverty_ref2", CBMA_SFX, ".rds")))

S_LAB <- c(raw = "Raw", t = "T-score", log = "log(1+raw)")
W_LAB <- c(ch = "Calinski-Harabasz", sdbw = "S_Dbw", equal = "Equal", avg = "Averaged")
tab <- out |>
  mutate(scale = factor(S_LAB[scale], S_LAB), weighting = factor(W_LAB[weighting], W_LAB)) |>
  arrange(scale, weighting, model, class) |>
  transmute(scale, weighting, model, class, pct = round(100 * n / n_model, 1), class_mean = round(class_mean, 1),
            rrr = round(rrr, 2), ci = sprintf("%.2f-%.2f", lo, hi),
            ref_pct = round(100 * ref_n / n_model, 1), ref_mean = round(ref_mean, 1))
sumf <- file.path(SUMM, paste0("summary_cbma_poverty_ref2", CBMA_SFX, ".txt"))
txt <- c(sprintf("clusterBMA consensus fits (%s), log and T-score scales: poverty RRR per halving of income, every class vs class 2 (%s).",
                 CBMA_SET_LAB, format(Sys.time(), "%Y-%m-%d %H:%M")),
         "Class 1 is the near-zero lowest class; class 2 is the reference here. adj: age, sex, site and race/ethnicity (primary); unadj: without race/ethnicity.",
         "pct: % of the children in the model; class_mean / ref_mean: mean raw CBCL over waves. Classes vs the lowest are in 12b_cbma_poverty.R.", "",
         capture.output(print(tab, n = 200, width = 200)))
writeLines(txt, sumf); cat(txt, sep = "\n")
