#!/usr/bin/env bash
# get_future.sh - ISIMIP3b bias-adjusted daily climate (historical 1901-2014 + SSP 2015-2100) cut to Lishui,
#                 merged into one 1901-2100 file per GCM x SSP x variable, plus matching CO2 files.
# Usage:  GCMS="gfdl-esm4 ukesm1-0-ll" SSPS="ssp126 ssp585" bash get_future.sh
# Safe to re-run: finished pieces are skipped.
set -uo pipefail
BASE=$HOME/lpjg_lishui; RAW=$BASE/data/raw3b; OUT=$BASE/data/climate3b; CO2D=$BASE/data/co2
mkdir -p "$RAW" "$OUT" "$CO2D"
API=https://files.isimip.org/api/v2
ROOT=ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily
CROOT=ISIMIP3b/InputData/climate/atmosphere_composition/co2
BBOX='[118.5, 120.5, 27.0, 29.5]'          # west east south north, covers Lishui
GCMS=${GCMS:-"gfdl-esm4 ukesm1-0-ll"}
SSPS=${SSPS:-"ssp126 ssp585"}
VARS="tas pr rsds"

upper() { echo "$1" | tr 'a-z' 'A-Z'; }
member() { if [ "$1" = ukesm1-0-ll ]; then echo r1i1p1f2; else echo r1i1p1f1; fi; }
periods() {
  if [ "$1" = historical ]; then
    for s in 1901 1911 1921 1931 1941 1951 1961 1971 1981 1991 2001; do echo "${s}_$((s + 9))"; done; echo 2011_2014
  else
    echo 2015_2020; for s in 2021 2031 2041 2051 2061 2071 2081 2091; do echo "${s}_$((s + 9))"; done
  fi
}
parse() {  # status|done|total|job_url|file_url|errors
  python3 -c '
import sys, json
try:
    d = json.load(sys.stdin)
except Exception:
    print("noresponse|||||"); sys.exit()
m = d.get("meta") or {}
e = d.get("errors")
print("|".join(str(x) for x in (d.get("status", ""), m.get("created_files", "?"), m.get("total_files", "?"),
      d.get("job_url", ""), d.get("file_url", ""), json.dumps(e) if e else "")))'
}
body_for() {  # body_for gcm experiment variable
  local g=$1 x=$2 v=$3 m p="" per
  m=$(member "$g")
  for per in $(periods "$x"); do
    p+="\"$ROOT/$x/$(upper "$g")/${g}_${m}_w5e5_${x}_${v}_global_daily_${per}.nc\","
  done
  echo "{\"paths\":[${p%,}],\"operations\":[{\"operation\":\"select_bbox\",\"bbox\":$BBOX}]}"
}

# ---- 1. submit all cut-out jobs (historical once per GCM, plus each SSP) ----
declare -A JOB
for g in $GCMS; do for x in historical $SSPS; do for v in $VARS; do
  key=${g}_${x}_${v}
  [ -s "$RAW/$key.zip" ] && continue
  IFS='|' read -r st n t job_url file_url err <<< "$(curl -sS --retry 5 --max-time 120 \
      -H 'Content-Type: application/json' -d "$(body_for "$g" "$x" "$v")" "$API" | parse)"
  if [ -n "$err" ] || [ -z "$job_url" ]; then echo "[$key] submit failed: $st $err"; continue; fi
  JOB[$key]=$job_url; echo "[$key] submitted ($st)"
done; done; done

# ---- 2. wait, download ----
for key in "${!JOB[@]}"; do
  while :; do
    IFS='|' read -r st n t job_url file_url err <<< "$(curl -sS --retry 5 --max-time 120 "${JOB[$key]}" | parse)"
    [ "$st" = finished ] && [ -n "$file_url" ] && break
    [ "$st" = failed ] && { echo "[$key] server failed: $err"; continue 2; }
    echo "[$key] $(date +%T) status=$st $n/$t"; sleep 60
  done
  if curl -sSL --retry 5 -C - -o "$RAW/$key.zip.part" "$file_url"; then mv "$RAW/$key.zip.part" "$RAW/$key.zip"; echo "[$key] downloaded"
  else echo "[$key] download failed, re-run to resume"; fi
done

# ---- 3. merge historical + SSP into one 1901-2100 file per GCM x SSP x variable ----
for g in $GCMS; do for s in $SSPS; do for v in $VARS; do
  out=$OUT/${g}_${s}_${v}_1901_2100.nc
  [ -s "$out" ] && continue
  if [ ! -s "$RAW/${g}_historical_${v}.zip" ] || [ ! -s "$RAW/${g}_${s}_${v}.zip" ]; then echo "[$g $s $v] pieces missing"; continue; fi
  tmp=$RAW/tmp_${g}_${s}_${v}; rm -rf "$tmp"; mkdir -p "$tmp"
  unzip -q -o "$RAW/${g}_historical_${v}.zip" -d "$tmp/h"; unzip -q -o "$RAW/${g}_${s}_${v}.zip" -d "$tmp/f"
  mapfile -t ncs < <(find "$tmp" -name '*.nc' | sort)
  if cdo -s -O mergetime "${ncs[@]}" "$out"; then
    echo "[$g $s $v] merged ${#ncs[@]} files -> $(basename "$out"), years $(cdo -s showyear "$out" | awk '{print $1 "-" $NF}')"
  else echo "[$g $s $v] merge failed"; fi
  rm -rf "$tmp"
done; done; done

# ---- 4. CO2: historical 1850-2014 + SSP 2015-2100 ----
getco2() {  # getco2 experiment period
  local f=$CO2D/co2_$1_annual_$2.txt
  [ -s "$f" ] || curl -sSfL --retry 5 -o "$f" "https://files.isimip.org/$CROOT/$1/co2_$1_annual_$2.txt" || { rm -f "$f"; return 1; }
}
getco2 historical 1850_2014 || echo "CO2 historical download failed: find it on data.isimip.org (ISIMIP3b, co2, historical)"
for s in $SSPS; do
  getco2 "$s" 2015_2100 || { echo "CO2 $s download failed: find it on data.isimip.org (ISIMIP3b, co2, $s)"; continue; }
  cat "$CO2D/co2_historical_annual_1850_2014.txt" "$CO2D/co2_${s}_annual_2015_2100.txt" 2>/dev/null \
    | awk 'NF >= 2 && $1 + 0 == $1 {printf "%d %.3f\n", $1, $2}' > "$CO2D/co2_${s}_1850_2100_lpjg.txt"
  echo "CO2 $s: $(head -n 1 "$CO2D/co2_${s}_1850_2100_lpjg.txt" | cut -d' ' -f1)-$(tail -n 1 "$CO2D/co2_${s}_1850_2100_lpjg.txt" | cut -d' ' -f1)"
done
echo "done $(date +%T)"; ls -lh "$OUT"
