#!/usr/bin/env bash
# 下载丽水范围的 ISIMIP3a GSWP3-W5E5 逐日气候（在 ISIMIP 服务器端裁剪）和 CO2 年值
set -uo pipefail
BASE=$HOME/lpjg_lishui
RAW=$BASE/data/raw; CLIM=$BASE/data/climate; CO2D=$BASE/data/co2
mkdir -p "$RAW" "$CLIM" "$CO2D"

API=https://files.isimip.org/api/v2
DIR=ISIMIP3a/InputData/climate/atmosphere/obsclim/global/daily/historical/GSWP3-W5E5
BBOX='[118.5, 120.5, 27.0, 29.5]'   # 西 东 南 北，覆盖丽水市域
VARS=${VARS:-"tas pr rsds"}          # 以后要用 BLAZE 火模型：VARS="hurs sfcwind" bash get_data.sh

# 解析 API 返回的 JSON，输出：状态|已处理|总数|job_url|file_url|错误
parse() {
  python3 -c '
import sys, json
try:
    d = json.load(sys.stdin)
except Exception:
    print("noresponse|||||"); sys.exit()
m = d.get("meta") or {}
e = d.get("errors")
print("|".join(str(x) for x in (d.get("status", ""), m.get("created_files", "?"),
      m.get("total_files", "?"), d.get("job_url", ""), d.get("file_url", ""),
      json.dumps(e, ensure_ascii=False) if e else "")))'
}

# 某变量 1901-2019 共 12 个十年文件的裁剪请求
body_for() {
  local v=$1 p="" s e
  for s in 1901 1911 1921 1931 1941 1951 1961 1971 1981 1991 2001 2011; do
    e=$((s + 9)); [[ $s == 2011 ]] && e=2019
    p+="\"$DIR/gswp3-w5e5_obsclim_${v}_global_daily_${s}_${e}.nc\","
  done
  echo "{\"paths\":[${p%,}],\"operations\":[{\"operation\":\"select_bbox\",\"bbox\":$BBOX}]}"
}

# ---------- 1. 一次性提交所有变量的裁剪任务（重复提交同一请求只会返回已有任务）----------
declare -A JOB
for v in $VARS; do
  [[ -s $CLIM/${v}_lishui_daily_1901_2019.nc || -s $RAW/$v.zip ]] && continue
  IFS='|' read -r st n t job_url file_url err <<< "$(curl -sS --retry 5 --max-time 120 \
      -H 'Content-Type: application/json' -d "$(body_for "$v")" "$API" | parse)"
  if [[ -n $err || -z $job_url ]]; then echo "[$v] 提交失败：$st $err"; exit 1; fi
  JOB[$v]=$job_url
  echo "[$v] 已提交，状态 $st：$job_url"
done

# ---------- 2. 等待服务器处理，下载，合并成一个文件 ----------
for v in $VARS; do
  out=$CLIM/${v}_lishui_daily_1901_2019.nc
  zip=$RAW/$v.zip
  if [[ -s $out ]]; then echo "[$v] 已存在，跳过"; continue; fi
  if [[ ! -s $zip ]]; then
    while :; do
      IFS='|' read -r st n t job_url file_url err <<< "$(curl -sS --retry 5 --max-time 120 "${JOB[$v]}" | parse)"
      [[ $st == finished && -n $file_url ]] && break
      [[ $st == failed ]] && { echo "[$v] 服务器处理失败：$err"; exit 1; }
      echo "[$v] $(date +%T) 状态=$st 已处理 $n/$t 个文件，30 秒后再查"
      sleep 30
    done
    echo "[$v] 下载 $file_url"
    if curl -sSL --retry 5 -C - -o "$zip.part" "$file_url"; then mv "$zip.part" "$zip"
    else echo "[$v] 下载失败，重新运行本脚本会续传"; exit 1; fi
  fi
  rm -rf "${RAW:?}/$v"; mkdir -p "$RAW/$v"
  unzip -q -o "$zip" -d "$RAW/$v"
  mapfile -t ncs < <(find "$RAW/$v" -name '*.nc' | sort)
  echo "[$v] 合并 ${#ncs[@]} 个文件（应为 12 个）"
  cdo -s -O mergetime "${ncs[@]}" "$out" || { echo "[$v] 合并失败"; exit 1; }
  ncdump -h "$out" | grep -E "(${v}:units|${v}:standard_name|time:calendar)"
  cdo -s showyear "$out" | awk '{print "  年份", $1, "-", $NF, "（共", NF, "年）"}'
  cdo -s griddes "$out" | grep -E '^(xsize|ysize|xfirst|yfirst)'
done

# ---------- 3. CO2 年值（1850-2021）----------
co2raw=$CO2D/co2_obsclim_annual_1850_2021.txt
if [[ ! -s $co2raw ]]; then
  curl -sSfL --retry 5 -o "$co2raw" "https://files.isimip.org/ISIMIP3a/InputData/climate/atmosphere_composition/co2/obsclim/co2_obsclim_annual_1850_2021.txt" \
  || { u=$(curl -sSf "https://data.isimip.org/api/v1/files/61a859b6-8de3-465a-936a-dfd837f42e1e/" \
           | python3 -c 'import sys, json; print(json.load(sys.stdin)["file_url"])') \
       && curl -sSfL --retry 5 -o "$co2raw" "$u"; } \
  || echo "CO2 下载失败"
fi
if [[ -s $co2raw ]]; then
  awk 'NF >= 2 && $1 + 0 == $1 {printf "%d %.3f\n", $1, $2}' "$co2raw" > "$CO2D/co2_1850_2021_lpjg.txt"
  awk 'NR == 1 {a = $1} END {print "CO2 年份", a, "-", $1, "（共", NR, "年）"}' "$CO2D/co2_1850_2021_lpjg.txt"
fi
echo "全部完成：$(date +%T)。数据在 $BASE/data/"
