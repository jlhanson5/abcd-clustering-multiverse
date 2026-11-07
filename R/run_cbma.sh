#!/bin/zsh
# R/run_cbma.sh: the nine clusterBMA consensus fits (12_cbma.R), two at a time, then
# 13_cbma_summary.R. Run from the project folder: caffeinate -i zsh R/run_cbma.sh
# Two at a time keeps memory near 8 GB on a 16 GB Mac (about 3.7 GB per fit at n = 8,331).
# Uses Apple's vecLib BLAS for this process only, if CBMA_VECLIB points to a folder
# holding a libRblas.dylib symlink to R's libRblas.vecLib.dylib (results match the
# reference BLAS to rounding; see clusterBMA_port NOTES.md, Step 6). Otherwise R's default.
# Finished fits (output/cbma/<scale>_<weighting>.rds) are skipped, so it can be restarted.

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
pairs=("raw ch" "t ch" "log ch" "raw sdbw" "t sdbw" "log sdbw" "raw equal" "t equal" "log equal")
for ((i = 1; i <= ${#pairs}; i += 2)); do
  fit ${=pairs[i]} &
  (( i + 1 <= ${#pairs} )) && fit ${=pairs[i+1]} &
  wait
done
Rscript --vanilla 13_cbma_summary.R > output/logs/cbma_summary.log 2>&1 && echo "summary done"
