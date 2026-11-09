# Archive (moved 2026-11-09)

Code that is no longer part of the pipeline. Nothing in the 01-13 scripts or `R/` sources
these files. They expect to be run from the project root and write to their old `output/`
paths, so move outputs back (or edit the paths) before rerunning one.

| File | What it was |
|---|---|
| `code/check_predict_lcmm.R` | Simulated-data check that `predict_lcmm()` in R/predict_class.R works (STABILITY_PLAN "Open items") |
| `code/clear_failed_age.R` | One-off (2026-10-23): deleted full-run age-metric task files whose fits all failed on the child with no interview_age, so 03_run.R refit them |
| `code/gmm_optima.R`, `code/gmm_optima_summary.R` | GMM optima experiment (2026-11-03), which led to `R/gmm_refit.R`. Outputs: `output/archive/gmm_optima_experiment/` |

Superseded outputs (not in git) are in `output/archive/`: `pilot/`, `smoke/` (stability and
clusterBMA smoke tests, plus the 10,000-step clusterBMA convergence check),
`stability_archive_2026-10-25/` and `gmm_optima_experiment/`. Old logs are in
`output/logs/archive/`.
