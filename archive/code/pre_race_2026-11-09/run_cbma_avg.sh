#!/bin/zsh
# R/run_cbma_avg.sh: after R/run_cbma.sh, the averaged-weight consensus for each scale
# (12_cbma.R <scale> avg, which reads that scale's ch, sdbw and equal runs), the
# pooled all-scales equal-weight consensus and the k-balanced consensus for each scale
# (added 2026-11-08), two at a time; then 13_cbma_summary.R again.
# Run from the project folder: CBMA_VECLIB=<folder> caffeinate -i zsh R/run_cbma_avg.sh
# (CBMA_VECLIB as in R/run_cbma.sh). Added 2026-11-07.
# Finished fits are skipped, so it can be restarted.

RH=/Library/Frameworks/R.framework/Resources
fit() {
  [[ -f output/cbma/$1_$2.rds ]] && { echo "skip $1 $2"; return; }
  if [[ -n "$CBMA_VECLIB" ]]; then
    /usr/bin/env DYLD_LIBRARY_PATH=$CBMA_VECLIB R_HOME=$RH R_ARCH= $RH/bin/exec/R --no-echo --no-restore \
      --vanilla --file=12_cbma.R --args $1 $2 > output/logs/cbma_$1_$2.log 2>&1
  else
    Rscript --vanilla 12_cbma.R $1 $2 > output/logs/cbma_$1_$2.log 2>&1
  fi
  echo "$(date '+%m-%d %H:%M')  $1 $2: $(tail -2 output/logs/cbma_$1_$2.log | tr '\n' ' ')"
}
fit raw avg & fit t avg & wait
fit log avg & fit all equal & wait
fit raw kbal & fit t kbal & wait
fit log kbal
Rscript --vanilla 13_cbma_summary.R > output/logs/cbma_summary.log 2>&1 && echo "summary done"
