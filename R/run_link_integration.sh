#!/bin/zsh
# R/run_link_integration.sh: add the skew-aware GBTM (gbtm_link, added 2026-11-13) to the
# primary run (full_refit), the stability analysis and clusterBMA, after
#   Rscript 03_run.R full 6 gbtm_link
# has written output/universes/full/gbtm_link_*.rds. Run from the project folder:
#   caffeinate -i zsh R/run_link_integration.sh > output/logs/link_integration.log 2>&1
# Before anything is overwritten, the current results, summaries, figures, stability
# tables and the clusterBMA fits that change (raw scale and pooled) are moved or copied to
# output/archive/pre_link_2026-11-13/.
# Phase 1 (sequential): 04, 05, 06 for full_refit; 11, 15, 16; 07 (stability table and base
#   solutions, now with gbtm_link on the raw anchor).
# Phase 2 (in parallel, at most 7 cores):
#   A  08 for gbtm_link (5 cores), then 09 and 10
#   B  clusterBMA raw ch/sdbw/equal/avg/kbal and all equal (two at a time), 12b, 13, 14;
#      then the permutation top-up (17, 3 cores) and its summary; then 18
setopt err_exit
ARC=output/archive/pre_link_2026-11-13
stamp() { echo "$(date '+%m-%d %H:%M')  $*"; }

[[ $(ls output/universes/full/gbtm_link_*.rds | wc -l) -eq 12 ]] || { echo "gbtm_link fits missing"; exit 1; }
mkdir -p $ARC/cbma $ARC/stability
cp -Rp output/results output/summaries output/figures $ARC/
cp -p output/stability/spec_reps.rds output/stability/base_assign.rds output/stability/null_model.rds $ARC/stability/
for f in raw_ch raw_sdbw raw_equal raw_avg raw_kbal all_equal; do
  [[ -f output/cbma/$f.rds ]] && mv output/cbma/$f.rds $ARC/cbma/
done
cp -p output/universes/full/gbtm_link_*.rds output/universes/full_refit/
stamp "archived; gbtm_link copied into full_refit"

# ---- Phase 1 -------------------------------------------------------------------------------
for s in "04_poverty.R full_refit" "05_agreement.R full_refit" "06_figures.R full_refit" \
         "11_refit_figures.R" "15_class_gradient.R" "16_sens_site_re.R 8" "07_stability_spec.R full"; do
  Rscript --vanilla ${=s} > output/logs/link_${${s%% *}%.R}.log 2>&1
  stamp "$s done"
done
# 07 refits the null model; it must equal the one the stability replicates used.
Rscript -e 'a <- readRDS("output/stability/null_model.rds"); b <- readRDS("output/archive/pre_link_2026-11-13/stability/null_model.rds"); stopifnot(isTRUE(all.equal(a, b)))'
stamp "null model unchanged"

# ---- Phase 2 -------------------------------------------------------------------------------
chain_a() {
  Rscript --vanilla 08_stability_run.R full 5 method=gbtm_link > output/logs/link_08.log 2>&1; stamp "08 gbtm_link done"
  Rscript --vanilla 09_stability_summary.R full 5 > output/logs/link_09.log 2>&1; stamp "09 done"
  Rscript --vanilla 10_stability_figures.R full > output/logs/link_10.log 2>&1; stamp "10 done"
}
cbma() {
  Rscript --vanilla 12_cbma.R $1 $2 > output/logs/cbma/cbma_$1_$2.log 2>&1
  stamp "cbma $1 $2: $(tail -2 output/logs/cbma/cbma_$1_$2.log | tr '\n' ' ')"
}
chain_b() {
  cbma raw ch & cbma raw sdbw & wait
  cbma raw equal & cbma all equal & wait
  cbma raw avg & cbma raw kbal & wait
  Rscript --vanilla 12b_cbma_poverty.R > output/logs/cbma/cbma_poverty.log 2>&1; stamp "12b done"
  Rscript --vanilla 13_cbma_summary.R > output/logs/cbma/cbma_summary.log 2>&1; stamp "13 done"
  Rscript --vanilla 14_extra_figures.R > output/logs/link_14.log 2>&1; stamp "14 done"
  Rscript --vanilla 17_permutation.R 500 3 > output/logs/link_17.log 2>&1; stamp "17 top-up done"
  Rscript --vanilla 17_permutation.R summary >> output/logs/link_17.log 2>&1; stamp "17 summary done"
  Rscript --vanilla 18_sens_weights.R 1 > output/logs/link_18.log 2>&1; stamp "18 done"
}
chain_a &
chain_b &
wait
stamp "LINK INTEGRATION DONE"
