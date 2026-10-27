#!/usr/bin/env bash
# run_batch.sh：密度 × 经营 × 根系深度 × 气候情景 的批量试验
# 可以重复运行：已经跑完（run.log 里有 Finished）的会自动跳过
B=$HOME/lpjg_lishui; E=$B/exp_base; GUESS2=$(cat $B/.guess2_bin); X=$B/batch
DENS="1667 2500 3333 5000"          # 初植密度（株/ha）
MGMT="A B"                          # A 不间伐；B 规程间伐（第 9 年 28%、第 13 年 24%）
BETAS="0.95 0.96 0.97"              # 杉木根系深度参数，0.96 为中心值
SCEN="ctrl"; for y in 1993 2000 2008; do SCEN="$SCEN dry50_$y dry50x2_$y hot3_$y comp_$y"; done
mkdir -p $X/runs

mkrun() {  # mkrun 目录 密度 经营 根系参数 情景
  local d=$1
  rm -rf $d && cp -r $E $d && rm -rf $d/output && mkdir $d/output
  sed -i -E "s#plantdensity \"[0-9]+\"#plantdensity \"$2\"#" $d/lishui_fir.ins
  if [ "$3" = B ]; then
    sed -i -E 's#^(\s*thintime\s+).*#\1 0.35 0.5 0 0 0#; s#^(\s*thinstr\s+).*#\1 0.28 0.24 0 0 0#; s#^(\s*thinstrunsel\s+).*#\1 0.95 0 0 0 0#; s#^(\s*thinselectdiam\s+).*#\1 1 1 0 0 0#' $d/lishui_fir.ins
  fi
  awk -v line="root_beta $4" '/^[[:space:]]*pft[[:space:]]*"TeNE"/ {inb = 1}
       inb && /^[[:space:]]*\)/ {print "\t" line "   ! Lishui drought calibration"; inb = 0}
       {print}' $d/lishui_global.ins > $d/tmp.ins && mv $d/tmp.ins $d/lishui_global.ins
  if [ "$5" != ctrl ]; then
    local S=$B/data/climate_scen/$5
    sed -i -E "s#^(\s*param\s+\"file_prec\"\s*\(\s*str\s+)\"[^\"]*\"#\1\"$S/pr.nc\"#; s#^(\s*param\s+\"file_temp\"\s*\(\s*str\s+)\"[^\"]*\"#\1\"$S/tas.nc\"#" $d/lishui_cf.ins
  fi
}

echo "run,density,mgmt,beta,scenario" > $X/manifest.csv
todo=()
for den in $DENS; do for m in $MGMT; do for b in $BETAS; do for s in $SCEN; do
  name=d${den}_${m}_b${b#0.}_$s
  echo "$name,$den,$m,$b,$s" >> $X/manifest.csv
  if ! grep -q Finished $X/runs/$name/run.log 2>/dev/null; then
    mkrun $X/runs/$name $den $m $b $s; todo+=($name)
  fi
done; done; done; done
echo "共 $(( $(wc -l < $X/manifest.csv) - 1 )) 个运行，本次要跑 ${#todo[@]} 个，开始于 $(date +%T)"
printf '%s\n' "${todo[@]}" | xargs -r -P 4 -I{} sh -c \
  'cd "$1/{}" && "$0" -input cf lishui_fir.ins > run.log 2>&1; echo "$(date +%T) 完成 {}"' "$GUESS2" "$X/runs"
echo "全部结束于 $(date +%T)"
bad=$(grep -L Finished $X/runs/*/run.log 2>/dev/null | wc -l); echo "未正常结束的运行：$bad 个"
grep -L Finished $X/runs/*/run.log 2>/dev/null | head -n 5
