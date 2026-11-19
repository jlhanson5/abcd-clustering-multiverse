#!/bin/zsh
# R/run_cbma_k28.sh: the clusterBMA fits on the k = 2-8 input sets (12_cbma.R).
#   CBMA_SET=rule (default; primary, decided 2026-11-18): k-means and kml x k 2-8 x one per
#     family / all children = 28 inputs per scale, output/cbma/.
#   CBMA_SET=four (sensitivity): plus PAM and GBTM, 56 per scale (raw 54), output/cbma_four/.
# K_BMA = 8. Per scale: Calinski-Harabasz, S_Dbw, equal, then averaged (started only once
# that scale's other three exist), plus the pooled all-scales run with equal weights.
# At most two fits at a time (about 4-5 GB each on this 16 GB Mac), counting any 12_cbma.R
# already running (a fit already running is not started again). The second starts only
# when memory_pressure reports >= MIN_FREE % free and the compressor is under MAX_COMP MB
# on 3 checks 5 minutes apart, so other projects' jobs are counted. While another
# project's targets/crew job runs, or swap is growing, fits run one at a time (2026-11-19:
# two fits plus the crew workers thrashed the Mac). SKIP="raw sdbw,log
# avg" leaves fits out (here for a decision still open). For the rule set, 12b_cbma_poverty.R,
# 13_cbma_summary.R, 14_extra_figures.R and 15_class_gradient.R run afterwards.
# Run from the project folder: caffeinate -i zsh R/run_cbma_k28.sh
# The 2026-11-08 all-methods fits were moved to output/cbma_allmethods/ and their summaries
# and figures to output/archive/pre_k28_2026-11-18/ before the first run.
# Finished fits (<dest>/<scale>_<weighting>.rds) are skipped, so it can be restarted.
# CBMA_VECLIB as in R/run_cbma.sh (unset: R's default BLAS).

export CBMA_SET=${CBMA_SET:-rule}
case $CBMA_SET in rule) D=output/cbma ;; four) D=output/cbma_four ;; *) echo "CBMA_SET rule or four"; exit 1 ;; esac
MIN_FREE=${MIN_FREE:-50}; MAX_COMP=${MAX_COMP:-3000}
RH=/Library/Frameworks/R.framework/Resources
L=output/logs/cbma_$CBMA_SET; mkdir -p $L
fit() {
  if [[ -n "$CBMA_VECLIB" ]]; then
    /usr/bin/env DYLD_LIBRARY_PATH=$CBMA_VECLIB R_HOME=$RH R_ARCH= $RH/bin/exec/R --no-echo --no-restore \
      --vanilla --file=12_cbma.R --args $1 $2 > $L/cbma_$1_$2.log 2>&1
  else
    Rscript --vanilla 12_cbma.R $1 $2 > $L/cbma_$1_$2.log 2>&1
  fi
  echo "$(date '+%m-%d %H:%M')  $1 $2 done: $(tail -2 $L/cbma_$1_$2.log | tr '\n' ' ')"
}
headroom() {
  local free comp
  free=$(memory_pressure | awk '/free percentage/ {gsub("%",""); print $NF}')
  comp=$(top -l 1 | awk '/PhysMem/ {for (i=1;i<=NF;i++) if ($(i+1) ~ /compressor/) print $i}' | tr -d '(M')
  [[ $comp == *G ]] && comp=$(( ${comp%G} * 1024 ))
  local sw; sw=$(sysctl -n vm.swapusage | awk '{gsub("M","",$6); print int($6)}')
  local grow=$(( sw > ${LAST_SW:-999999} )); LAST_SW=$sw
  # A second fit only when no targets/crew job is running and swap is not growing.
  pgrep -f "targets::tar_make|crew::crew_worker" >/dev/null && return 1
  (( !grow && free >= MIN_FREE && comp < MAX_COMP ))
}
running() { pgrep -f "file=12_cbma.R" | wc -l | tr -d ' '; }
echo "$(date '+%m-%d %H:%M')  start, set $CBMA_SET${SKIP:+, skipping $SKIP}"
for job in "raw ch" "t ch" "log ch" "raw sdbw" "t sdbw" "log sdbw" "raw equal" "t equal" "log equal" \
           "all equal" "raw avg" "t avg" "log avg"; do
  s=${job% *}; w=${job#* }
  [[ -f $D/${s}_$w.rds ]] && { echo "skip $s $w (done)"; continue; }
  [[ ",$SKIP," == *",$job,"* ]] && { echo "skip $s $w (SKIP)"; continue; }
  pgrep -f "12_cbma.R --args $s $w" >/dev/null && { echo "skip $s $w (already running)"; continue; }
  if [[ $w == avg ]]; then
    until [[ -f $D/${s}_ch.rds && -f $D/${s}_sdbw.rds && -f $D/${s}_equal.rds ]]; do sleep 120; done
  fi
  # Wait for a slot: no fit running, or one running and headroom on 3 checks in a row.
  ok=0
  while true; do
    n=$(running)
    (( n == 0 )) && break
    if (( n == 1 )) && headroom; then (( ++ok >= 3 )) && break; else ok=0; fi
    sleep 300
  done
  echo "$(date '+%m-%d %H:%M')  launch $s $w"
  fit $s $w &
  sleep 300   # let the new fit reach its working size before judging headroom again
done
while (( $(running) > 0 )); do sleep 300; done
n=$(ls $D/*.rds 2>/dev/null | wc -l | tr -d ' ')
need=$(( 13 - $(echo ${SKIP//[^,]/} | wc -c) )); [[ -z $SKIP ]] && need=13
[[ $n -ge $need ]] || { echo "only $n of $need fits; stopping before the summaries"; exit 1; }
[[ -n $SKIP ]] && { echo "$(date '+%m-%d %H:%M')  fits done except $SKIP; summaries wait for those"; exit 0; }
[[ $CBMA_SET == rule ]] || { echo "$(date '+%m-%d %H:%M')  all fits done"; exit 0; }
Rscript --vanilla 12b_cbma_poverty.R > $L/cbma_poverty.log 2>&1 && echo "12b done" || exit 1
Rscript --vanilla 13_cbma_summary.R > $L/cbma_summary.log 2>&1 && echo "13 done" || exit 1
Rscript --vanilla 14_extra_figures.R > output/logs/extra_figures.log 2>&1 && echo "14 done" || exit 1
Rscript --vanilla 15_class_gradient.R > output/logs/class_gradient.log 2>&1 && echo "15 done" || exit 1
echo "$(date '+%m-%d %H:%M')  all done"
