#!/usr/bin/env bash
# run_rotland.sh - rotation length x planting phase x climate ("normal forest" landscape means)
# Safe to re-run: runs whose run.log contains "Finished" are skipped
B=$HOME/lpjg_lishui; E=$B/exp_base; GUESS2=$(cat $B/.guess2_bin); RL=$B/rotland; mkdir -p $RL
mkrot2() {  # mkrot2 name rotation first_manage_year scenario
  local d=$RL/$1 R=$2
  rm -rf $d && cp -r $E $d && rm -rf $d/output && mkdir $d/output
  local t1=$(awk -v r=$R 'BEGIN {printf "%.4f", 9.5 / r}') t2=$(awk -v r=$R 'BEGIN {printf "%.4f", 13.5 / r}')
  sed -i -E "s#^(\s*firstmanageyear\s+)[0-9]+#\1 $3#; s#^(\s*cutinterval\s+)[0-9]+#\1 $R#; s#^(\s*thintime\s+).*#\1 $t1 $t2 0 0 0#; s#^(\s*thinstr\s+).*#\1 0.28 0.24 0 0 0#; s#^(\s*thinstrunsel\s+).*#\1 0.95 0 0 0 0#; s#^(\s*thinselectdiam\s+).*#\1 1 1 0 0 0#" $d/lishui_fir.ins
  awk '/^[[:space:]]*pft[[:space:]]*"TeNE"/ {inb = 1}
       inb && /^[[:space:]]*\)/ {print "\troot_beta 0.96   ! Lishui drought calibration"; inb = 0}
       {print}' $d/lishui_global.ins > $d/tmp.ins && mv $d/tmp.ins $d/lishui_global.ins
  if [ "$4" != ctrl ]; then
    local S=$B/data/climate_scen/$4
    sed -i -E "s#^(\s*param\s+\"file_prec\"\s*\(\s*str\s+)\"[^\"]*\"#\1\"$S/pr.nc\"#; s#^(\s*param\s+\"file_temp\"\s*\(\s*str\s+)\"[^\"]*\"#\1\"$S/tas.nc\"#" $d/lishui_cf.ins
  fi
}
todo=()
for R in 16 21 26 31 41 500; do
  nk=4; [ $R = 500 ] && nk=1
  for k in $(seq 0 $((nk - 1))); do
    fy=$(awk -v r=$R -v k=$k 'BEGIN {printf "%d", 1901 + int(k * r / 4 + 0.5)}')
    for sc in ctrl warm2_9019 dry30_9019; do
      n=r${R}_o${k}_$sc
      grep -q Finished $RL/$n/run.log 2>/dev/null && continue
      mkrot2 $n $R $fy $sc; todo+=($n)
    done
  done
done
echo "runs to do: ${#todo[@]}, start $(date +%T)"
printf '%s\n' "${todo[@]}" | xargs -r -P 4 -I{} sh -c \
  'cd "$1/{}" && "$0" -input cf lishui_fir.ins > run.log 2>&1; echo "$(date +%T) done {}"' "$GUESS2" "$RL"
echo "all finished $(date +%T); failed runs: $(grep -L Finished $RL/*/run.log | wc -l)"
