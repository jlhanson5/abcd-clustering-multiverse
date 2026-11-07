# =============================================================================
# 07_stability_spec.R: the stability replicate table (see STABILITY_PLAN.md).
#   Rscript 07_stability_spec.R full|smoke
# One row per replicate: anchor x method x perturbation x replicate (x settings
# variant). A replicate fits every k (2-6), like a task in 02_spec.R.
#   output/stability/spec_reps.rds    replicate table with seeds
#   output/stability/anchor.rds       anchor children (ids, raw matrix, ages, site)
#   output/stability/null_model.rds   1-class growth model parameters (aggregate)
#   output/stability/base_assign.rds  the primary run's (full_refit) anchor solutions (relabelled)
# smoke: a 600-child sample, 2 replicates per perturbation, 5 starts, written to
# output/stability_smoke/. For checking the pipeline end to end only.
# =============================================================================

source("00_setup.R")
source("R/perturb.R")

mode <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(mode)) mode <- "full"
stopifnot(mode %in% c("full", "smoke"))
smoke <- mode == "smoke"
SOUT <- file.path(OUT, if (smoke) "stability_smoke" else "stability")
dir.create(SOUT, recursive = TRUE, showWarnings = FALSE)

wide <- readRDS(file.path(DER, "wide.rds"))
base <- anchor_base(wide, smoke_n = if (smoke) 600 else NULL)
cat(sprintf("[%s] anchor: %d children, %d sites\n", mode, length(base$ids), n_distinct(base$site)))
cat("Raw SD by wave:", round(base$sd_w, 2), " SEM at r = .90 / .80:",
    round(mean(base$sd_w) * sqrt(1 - NOISE_R), 2), "\n")

nm <- fit_null_model(base)
cat(sprintf("Null model (log scale): beta = %s, SD(int, slope) = %s, residual SD = %.3f\n",
            paste(round(nm$beta, 3), collapse = ", "),
            paste(round(sqrt(diag(nm$Sigma)), 3), collapse = ", "), nm$sigma))

# ---- Replicate table -------------------------------------------------------------
reps_for <- function(ptype, method) {
  r <- REPS[REPS$ptype == ptype, ]
  n <- if (method %in% TRAJ) r$traj else r$cheap
  if (ptype == "seed" && !method %in% STOCH) n <- 0
  if (smoke) n <- min(n, 2)
  if (n == 0) integer(0) else seq_len(n)
}

grid <- expand_grid(anchor = STAB_ANCHORS, method = c(TRAJ, CHEAP),
                    ptype = setdiff(REPS$ptype, "setting")) |>
  mutate(rep = map2(ptype, method, reps_for)) |>
  unnest_longer(rep) |>
  mutate(variant = "base")

settings <- VARIANTS |>
  expand_grid(anchor = STAB_ANCHORS) |>
  mutate(ptype = "setting",
         rep = map(method, \(m) c(0L, if (m %in% CHEAP) reps_for("setting", m)))) |>
  unnest_longer(rep) |>
  select(anchor, method, ptype, rep, variant, args)

reps <- bind_rows(grid, settings) |>
  mutate(args = map(args, \(a) if (is.null(a)) list() else a),
         fit_seed = map2_int(ptype, rep, fit_seed),
         rep_id = sprintf("%s_%s_%s_%s_r%03d", anchor, method, ptype, variant, rep)) |>
  arrange(match(method, c("gmm", "gbtm", "kml", "pam", "llpa", "kmeans", "ward")), ptype, rep, anchor)
stopifnot(!anyDuplicated(reps$rep_id))

# ---- Base solutions from the main run ----------------------------------------------
# The anchor universes: one per family, ge3_imp, age metric for GBTM/GMM.
res <- readRDS(file.path(OUT, if (smoke) "results_full.rds" else "results_full_refit.rds"))   # primary run since 2026-11-07
uid <- expand_grid(anchor = STAB_ANCHORS, method = c(TRAJ, CHEAP), k = K_STAB) |>
  mutate(universe_id = sprintf("%s_%s_%s_one_ge3_imp_k%d", method, anchor,
                               ifelse(method %in% TRAJ, "age", "wave"), k))
base_assign <- lapply(uid$universe_id, function(u) {
  a <- res$assign[[u]]
  if (is.null(a)) return(NULL)
  a <- a[a$src_subject_id %in% base$ids, ]
  list(idx = match(a$src_subject_id, base$ids), cls = as.integer(a$class))
})
names(base_assign) <- uid$universe_id
rm(res); invisible(gc())
cat("Base solutions found:", sum(!vapply(base_assign, is.null, logical(1))), "of", nrow(uid), "\n")

saveRDS(reps, file.path(SOUT, "spec_reps.rds"))
saveRDS(base, file.path(SOUT, "anchor.rds"))
saveRDS(nm, file.path(SOUT, "null_model.rds"))
saveRDS(list(universes = uid, assign = base_assign), file.path(SOUT, "base_assign.rds"))

cat("\nReplicates:", nrow(reps), "\n")
print(count(reps, method, ptype) |> pivot_wider(names_from = method, values_from = n, values_fill = 0), n = 20)
