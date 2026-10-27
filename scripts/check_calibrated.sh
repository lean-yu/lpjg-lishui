#!/usr/bin/env bash
# check_calibrated.sh - two quick checks with the calibrated fir parameters (observed climate):
#  (1) rotation length vs harvest and carbon (6 rotations x 4 planting phases)
#  (2) save/restart test: switch management in 2000 (proxy for 2025) from a saved model state
# Re-runnable: finished runs are skipped. Summary: scripts/summarize_checks.py
B=$HOME/lpjg_lishui; E=$B/exp_base; GUESS2=$(cat $B/.guess2_bin); C=$B/checks; mkdir -p $C
PF=$B/scripts/tene_calibrated.ins
# state_year counts simulation years from the first spin-up year, not calendar years (framework.cpp)
NSPIN=$(cat $E/lishui_global.ins $E/lishui_cf.ins $E/lishui_fir.ins | awk '$1 == "nyear_spinup" {v = $2} END {print v}')
[ -n "$NSPIN" ] || { echo "nyear_spinup not found"; exit 1; }
SY=$((NSPIN + 2000 - 1901)); echo "nyear_spinup $NSPIN -> state_year $SY (start of calendar year 2000)"

ins_tene() {  # append the calibrated lines of $PF at the end of the TeNE block of file $1
  awk -v pf="$PF" '
    BEGIN {while ((getline l < pf) > 0) if (l !~ /^[ \t]*(!|$)/) add = add "\t" l "   ! calibrated (tene_calibrated.ins)\n"}
    /^[[:space:]]*pft[[:space:]]*"TeNE"/ {inb = 1}
    inb && /^[[:space:]]*\)/ {printf "%s", add; inb = 0}
    {print}' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}
set_rot() {  # set_rot dir rotation : rotation length with thinning at stand age 9 and 13
  local t1 t2
  t1=$(awk -v r=$2 'BEGIN {printf "%.4f", 9.5 / r}'); t2=$(awk -v r=$2 'BEGIN {printf "%.4f", 13.5 / r}')
  sed -i -E "s#^(\s*cutinterval\s+)[0-9]+#\1 $2#; s#^(\s*thintime\s+).*#\1 $t1 $t2 0 0 0#" $1/lishui_fir.ins
}
mk() {  # mk dir rotation first_manage_year
  rm -rf $1 && cp -r $E $1 && rm -rf $1/output $1/state && mkdir $1/output
  sed -i -E "s#^(\s*firstmanageyear\s+)[0-9]+#\1 $3#; s#^(\s*thinstr\s+).*#\1 0.28 0.24 0 0 0#; s#^(\s*thinstrunsel\s+).*#\1 0.95 0 0 0 0#; s#^(\s*thinselectdiam\s+).*#\1 1 1 0 0 0#" $1/lishui_fir.ins
  set_rot $1 $2
  ins_tene $1/lishui_global.ins
}
runs() { printf '%s\n' "$@" | xargs -r -P 4 -I{} sh -c 'cd "$1/{}" && "$0" -input cf lishui_fir.ins > run.log 2>&1; echo "$(date +%T) done {}"' "$GUESS2" "$C"; }
done_ok() { grep -q Finished $C/$1/run.log 2>/dev/null; }

# (1) rotations, plus the state-saving run for (2)
todo=()
for R in 16 21 26 31 36 41; do for k in 0 1 2 3; do
  n=rot_r${R}_p$k; done_ok $n && continue
  mk $C/$n $R $(awk -v r=$R -v k=$k 'BEGIN {printf "%d", 1901 + int(k * r / 4 + 0.5)}'); todo+=($n)
done; done
if ! done_ok rs_save; then
  mk $C/rs_save 26 1901; mkdir -p $C/rs_save/state
  printf '\n! save the model state at the start of 2000 (proxy for 2025)\nstate_path "%s/state/"\nsave_state 1\nstate_year %d\n' $C/rs_save $SY >> $C/rs_save/lishui_fir.ins
  todo+=(rs_save)
fi
echo "stage 1: ${#todo[@]} runs, start $(date +%T)"; runs "${todo[@]}"

# (2) branches restarted from the saved state; rs_bau must reproduce rot_r26_p0 after 2000
branch() {  # branch name : copy the saving run, switch to restart
  rm -rf $C/rs_$1 && cp -r $C/rs_save $C/rs_$1 && rm -rf $C/rs_$1/output $C/rs_$1/state && mkdir $C/rs_$1/output
  sed -i -E 's#^(\s*save_state\s+)1#\1 0#' $C/rs_$1/lishui_fir.ins
  echo "restart 1" >> $C/rs_$1/lishui_fir.ins
}
todo=()
done_ok rs_bau  || { branch bau; todo+=(rs_bau); }
done_ok rs_ext  || { branch ext; set_rot $C/rs_ext 36; todo+=(rs_ext); }
done_ok rs_stop || { branch stop
  sed -i -E 's#^(\s*cutinterval\s+)[0-9]+#\1 1000#; s#^(\s*thinstr\s+).*#\1 0 0 0 0 0#; s#^(\s*thinstrunsel\s+).*#\1 0 0 0 0 0#' $C/rs_stop/lishui_fir.ins
  awk '/^[[:space:]]*st[[:space:]]*"Fir_CC"/ {inb = 1}
       inb && /^[[:space:]]*\)/ {print "\treestab \"all\"   ! natural regeneration after harvest stops"; inb = 0} {print}' \
      $C/rs_stop/lishui_fir.ins > $C/rs_stop/tmp && mv $C/rs_stop/tmp $C/rs_stop/lishui_fir.ins
  todo+=(rs_stop); }
echo "stage 2: ${#todo[@]} runs, start $(date +%T)"; runs "${todo[@]}"
echo "all finished $(date +%T); failed runs: $(grep -L Finished $C/*/run.log 2>/dev/null | wc -l)"
grep -L Finished $C/*/run.log 2>/dev/null | head -n 5
