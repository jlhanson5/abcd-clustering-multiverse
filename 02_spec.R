# =============================================================================
# 02_spec.R: the specification table.
#   output/results/spec_tasks.rds     : one row per fitting task (method x data config);
#                               a task fits every k for that method at once
#   output/results/spec_universes.rds : one row per universe (task x k)
#
# Data configs (sample x missing-data handling):
#   ge3_imp    >= 3 of 4 waves, missing wave filled by the shared single imputation
#   ge3_native >= 3 of 4 waves, method's own handling: lcmm uses all observed
#              rows (FIML). For the vector methods this equals complete cases,
#              which is all4, so it is not a separate config for them.
#   all4       all 4 waves observed (no missing; imputed == native)
# Time metric (wave vs age) only changes GBTM/GMM. The vector methods and kml
# treat the 4 waves as equally spaced positions, so they get time = "wave" only.
# =============================================================================

source("00_setup.R")

VECTOR <- c("kmeans", "pam", "ward", "llpa", "kml")
TRAJ   <- c("gbtm", "gmm")        # gbtm_link: see lnk below
K_MAIN <- 2:6
K_SUPP <- 7:8            # supplement: the cheap methods only

vec <- expand_grid(method = VECTOR, scale = c("raw", "t", "log"), time = "wave",
                   family = c("one", "full"), config = c("ge3_imp", "all4"))
trj <- expand_grid(method = TRAJ, scale = c("raw", "t", "log"), time = c("wave", "age"),
                   family = c("one", "full"), config = c("ge3_imp", "ge3_native", "all4"))

# Skew-aware GBTM (lcmm with an estimated link; added 2026-11-13): raw scores only, since
# the link estimates the transformation. Appended after the original tasks so their
# task indices, and with them their seeds (R/fit_task.R), do not change.
lnk <- expand_grid(method = "gbtm_link", scale = "raw", time = c("wave", "age"),
                   family = c("one", "full"), config = c("ge3_imp", "ge3_native", "all4"))

tasks <- bind_rows(vec, trj, lnk) |>
  mutate(ks = map(method, \(m) if (m %in% VECTOR) c(K_MAIN, K_SUPP) else K_MAIN),
         task_id = paste(method, scale, time, family, config, sep = "_"))

universes <- tasks |>
  unnest_longer(ks, values_to = "k") |>
  mutate(supplement = k %in% K_SUPP,
         universe_id = paste0(task_id, "_k", k))

dir.create(RES, recursive = TRUE, showWarnings = FALSE)
saveRDS(tasks, file.path(RES, "spec_tasks.rds"))
saveRDS(universes, file.path(RES, "spec_universes.rds"))

cat("Tasks:", nrow(tasks), "\n")
cat("Universes, main (k 2-6):", sum(!universes$supplement),
    " supplement (k 7-8):", sum(universes$supplement), "\n")
print(count(universes, method, supplement))
