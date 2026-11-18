#!/bin/zsh
# R/run_cbma_k28.sh: the clusterBMA fits on the rule input set (decided 2026-11-15;
# 12_cbma.R with CBMA_SET=rule: k-means, kml, PAM, GBTM x k 2-8 x one per family / all
# children = 56 inputs per scale, K_BMA = 8). Per scale: Calinski-Harabasz, S_Dbw, equal,
# then averaged (needs the other three); plus the pooled all-scales run with equal weights
# (168 inputs). Two fits at a time (about 4-5 GB each on this 16 GB Mac), then
# 12b_cbma_poverty.R, 13_cbma_summary.R, 14_extra_figures.R and 15_class_gradient.R.
# Run from the project folder: caffeinate -i zsh R/run_cbma_k28.sh
# Before the first run, the all-methods fits (now the sensitivity analysis) were moved to
# output/cbma_allmethods/ and their summaries and figures to output/archive/pre_k28_2026-11-18/.
# Finished fits (output/cbma/<scale>_<weighting>.rds) are skipped, so it can be restarted.
# CBMA_VECLIB as in R/run_cbma.sh (unset: R's default BLAS).

export CBMA_SET=rule
RH=/Library/Frameworks/R.framework/Resources
mkdir -p output/logs/cbma
fit() {
  [[ -f output/cbma/$1_$2.rds ]] && { echo "skip $1 $2"; return; }
  if [[ -n "$CBMA_VECLIB" ]]; then
    /usr/bin/env DYLD_LIBRARY_PATH=$CBMA_VECLIB R_HOME=$RH R_ARCH= $RH/bin/exec/R --no-echo --no-restore \
      --vanilla --file=12_cbma.R --args $1 $2 > output/logs/cbma/cbma_$1_$2.log 2>&1
  else
    Rscript --vanilla 12_cbma.R $1 $2 > output/logs/cbma/cbma_$1_$2.log 2>&1
  fi
  echo "$(date '+%m-%d %H:%M')  $1 $2: $(tail -2 output/logs/cbma/cbma_$1_$2.log | tr '\n' ' ')"
}
echo "$(date '+%m-%d %H:%M')  start"
fit raw ch & fit t ch & wait
fit log ch & fit raw sdbw & wait
fit t sdbw & fit log sdbw & wait
fit raw equal & fit t equal & wait
fit log equal & fit all equal & wait
fit raw avg & fit t avg & wait
fit log avg
n=$(ls output/cbma/*.rds 2>/dev/null | wc -l | tr -d ' ')
[[ $n -eq 13 ]] || { echo "only $n of 13 fits; stopping before the summaries"; exit 1; }
Rscript --vanilla 12b_cbma_poverty.R > output/logs/cbma/cbma_poverty.log 2>&1 && echo "12b done" || exit 1
Rscript --vanilla 13_cbma_summary.R > output/logs/cbma/cbma_summary.log 2>&1 && echo "13 done" || exit 1
Rscript --vanilla 14_extra_figures.R > output/logs/extra_figures.log 2>&1 && echo "14 done" || exit 1
Rscript --vanilla 15_class_gradient.R > output/logs/class_gradient.log 2>&1 && echo "15 done" || exit 1
echo "$(date '+%m-%d %H:%M')  all done"
