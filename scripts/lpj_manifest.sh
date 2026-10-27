#!/usr/bin/env bash
# lpj_manifest.sh —— 只读诊断，不修改、不运行任何模型
# 汇总：系统资源、guess 可执行文件与版本、脚本和 .ins 模板、定制代码改动、
#       重启(save_state/restart)设置、runs_main 各类运行的完整性与耗时
#
# 用法（在腾讯云 Linux 上运行，不要在 Mac 本地运行）：
#   cd 到你执行 tree 的那一层（能看到 lpjg_lishui/ 和 src/ 的目录）
#   bash lpj_manifest.sh            # 或：bash lpj_manifest.sh /path/to/root
# 结果写入当前目录 lpj_manifest.txt。已自动遮蔽 secret/password/token/AKID 等字段，
# 上传前仍请快速浏览一遍。

ROOT="${1:-.}"
OUT="$(pwd)/lpj_manifest.txt"
if [ -d "$ROOT/lpjg_lishui" ]; then P="$ROOT/lpjg_lishui"; else P="$ROOT"; fi
RM="$P/runs_main"
SRC="$ROOT/src/guess"
[ -d "$SRC" ] || SRC="$(dirname "$P")/src/guess"

redact() {
  sed -E -e '/-----BEGIN [A-Z ]*PRIVATE KEY-----/,/-----END [A-Z ]*PRIVATE KEY-----/c\***PRIVATE KEY REDACTED***' \
         -e 's/((secret|passw(or)?d|token|api[_-]?key|access[_-]?key)[A-Za-z0-9_]*["\x27]?[[:space:]]*[=:][[:space:]]*)[^[:space:],]+/\1***REDACTED***/Ig' \
         -e 's/AKID[0-9A-Za-z]{8,}/AKID***REDACTED***/g'
}
sec()  { printf '\n==================== %s ====================\n' "$*"; }
listf() {  # 列出文件（不含目录）：listf 目录 深度 行数上限
  [ -d "$1" ] || { echo "(不存在) $1"; return; }
  find "$1" -maxdepth "${2:-2}" -type f -not -path '*/__pycache__/*' \
       -printf '%TY-%Tm-%Td %TH:%TM  %10s  %p\n' 2>/dev/null | sort -k4 | head -n "${3:-100}"
}
showf() {  # 打印文本文件（有行数上限）
  [ -f "$1" ] || return 0
  echo "---------- $1 ($(wc -l < "$1") 行)"
  head -n "${2:-300}" "$1"
}
KEYS='import|pft|pstemp|longevity|greff|respcoeff|alphaa|npatch|nyear_spinup|spinup|ifnlim|ifcentury|save_state|restart|state_path|state_year|harvest|cut|thin|firstmanage|switch|rotation|pulse|prf|elev|dT|co2|ndep|outputdirectory|file_'

{
sec "系统"
date; uname -a; echo "nproc=$(nproc)"; free -g; df -h "$ROOT" | tail -1
lscpu 2>/dev/null | grep -E 'Model name|Thread\(s\)|Core\(s\)|Socket'
command -v parallel >/dev/null && echo "GNU parallel: 有" || echo "GNU parallel: 无"
command -v tmux >/dev/null && echo "tmux: 有" || echo "tmux: 无"

sec "guess 可执行文件与源码版本"
find "$ROOT" -maxdepth 6 -type f -name 'guess' -perm -u+x -printf '%TY-%Tm-%Td %TH:%TM  %p\n' 2>/dev/null
if [ -d "$SRC" ]; then
  ( cd "$SRC" && { svnversion 2>/dev/null; svn info 2>/dev/null | head -12; git log -1 --oneline 2>/dev/null; } )
  echo "--- svn/git 本地改动文件（定制分支改了哪些文件）"
  ( cd "$SRC" && { svn status 2>/dev/null | grep -v '^?' | head -60; git status --short 2>/dev/null | head -60; } )
fi

sec "关键目录的文件清单"
for d in scripts pulse batch model_r13130 docs export trans_test \
         checks/rs_save checks/rs_bau checks/rs_stop checks/mt_bau rot/r500_ctrl r13130_nat r13130_fir; do
  echo "--- $P/$d"; listf "$P/$d" 2 60
done
for d in scripts run logs; do echo "--- $ROOT/lishui/$d"; listf "$ROOT/lishui/$d" 2 40; done
echo "--- lishui/cloud 与 lishui/mac（只列文件名，不读内容）"
find "$ROOT/lishui/cloud" "$ROOT/lishui/mac" -maxdepth 2 -type f -printf '%p\n' 2>/dev/null | head -40
echo "--- data（前 80 项）"
find "$P/data" -maxdepth 2 -printf '%y %10s  %p\n' 2>/dev/null | head -80
echo "--- 是否已有存档状态文件（state）"
find "$P" -maxdepth 4 \( -iname '*state*' -o -iname '*.bin' \) -printf '%10s  %p\n' 2>/dev/null | head -30

sec "脚本内容（scripts / pulse / batch 下的 .py .sh .R，每个最多 250 行）"
for f in $(find "$P/scripts" "$P/pulse" "$P/batch" -maxdepth 2 -type f \
             \( -name '*.py' -o -name '*.sh' -o -name '*.R' \) 2>/dev/null | sort | head -25); do
  showf "$f" 250
done

sec "settings 下 .ins 的关键参数行"
for f in $(find "$P/model_r13130/settings" -maxdepth 2 -type f -name '*.ins' 2>/dev/null | sort); do
  echo "---------- $f"; grep -nE "$KEYS" "$f" | head -80
done

sec "示例运行目录（对照 / 2045 热脉冲 / 一个 800 m 的 STOP 安慰剂）"
C="$RM/e500_bau26_t25_gfdl-esm4_ssp585_p1"
H="${C}_hot3y2045"
L="$RM/e800_stop_t30_ukesm1-0-ll_ssp585_p3_plcy2045"
for r in "$C" "$H" "$L"; do echo "--- $r"; listf "$r" 2 40; done
for f in "$C"/*.ins; do showf "$f" 400; done
echo "--- 对照与脉冲 .ins 的差异"
diff <(cat "$C"/*.ins 2>/dev/null) <(cat "$H"/*.ins 2>/dev/null)
echo "--- guess.log 头尾（对照）"
for g in "$C"/guess.log "$C"/*.log; do [ -f "$g" ] && { head -5 "$g"; echo ...; tail -8 "$g"; break; }; done

sec "rot/r500_ctrl 与 checks/rs_save 的 .ins（确认是否为观测气候不采伐长序列、重启设置）"
for f in "$P"/rot/r500_ctrl/*.ins "$P"/checks/rs_save/*.ins; do showf "$f" 200; done
echo "--- rs_save 与 rs_stop 的 .ins 差异"
diff <(cat "$P"/checks/rs_save/*.ins 2>/dev/null) <(cat "$P"/checks/rs_stop/*.ins 2>/dev/null)

sec "源码中的关键改动位置（grep）"
if [ -d "$SRC" ]; then
  grep -rnE 'pstemp_(min|low|high|max)' "$SRC/modules" "$SRC/framework" 2>/dev/null | head -40
  grep -rnE 'pulse|firstmanageyear|switch_year' "$SRC/modules" "$SRC/framework" 2>/dev/null | head -60
  grep -rnE 'save_state|state_year|serialize' "$SRC/framework" 2>/dev/null | head -30
  grep -rnE 'longevity|greff_min' "$SRC/modules" 2>/dev/null | head -20
fi
echo "--- model_r13130/code 下文件"; listf "$P/model_r13130/code" 3 60
echo "--- docs/patches 内容（前 1500 行）"
for f in $(find "$P/docs/patches" -maxdepth 2 -type f 2>/dev/null | sort); do showf "$f" 1500; done

sec "runs_main：完整性、磁盘与耗时估计"
du -sh "$RM" 2>/dev/null
python3 - "$RM" <<'PY'
import os, re, sys, statistics, collections, datetime
RM = sys.argv[1]
if not os.path.isdir(RM):
    print("runs_main 不存在"); sys.exit()
rows = []
for d in sorted(os.listdir(RM)):
    p = os.path.join(RM, d)
    if not os.path.isdir(p):
        continue
    files = []
    for root, _, fs in os.walk(p):
        for f in fs:
            files.append(os.path.join(root, f))
    outs = [f for f in files if not f.endswith('.ins')]
    mt = [os.path.getmtime(f) for f in outs]
    m = re.search(r'_(hot3y2045|hot3y2075|plcy2045)$', d)
    kind = m.group(1) if m else 'ctrl'
    parts = d.split('_')
    rows.append(dict(d=d, kind=kind, elev=parts[0], mgmt=parts[1] if len(parts) > 1 else '',
                     n=len(files), size=sum(os.path.getsize(f) for f in files),
                     t0=min(mt) if mt else None, t1=max(mt) if mt else None))
fmt = lambda t: datetime.datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M')
ref = [r for r in rows if r['kind'] == 'ctrl' and r['n'] > 0]
typ_n = statistics.median([r['n'] for r in ref]) if ref else 0
typ_s = statistics.median([r['size'] for r in ref]) if ref else 0
print(f"对照目录典型文件数 {typ_n}，典型大小 {typ_s/1e6:.1f} MB")
for k in ['ctrl', 'hot3y2045', 'hot3y2075', 'plcy2045']:
    rr = [r for r in rows if r['kind'] == k]
    ok = [r for r in rr if r['n'] >= typ_n and r['size'] >= 0.8 * typ_s]
    print(f"\n[{k}] 目录 {len(rr)} 个，看起来完整 {len(ok)} 个")
    c = collections.Counter((r['elev'], r['mgmt']) for r in ok)
    print("  完整数按 海拔×经营:", dict(sorted(c.items())))
    bad = sorted([r for r in rr if r not in ok], key=lambda r: r['size'])[:8]
    for r in bad:
        print(f"  不完整示例 {r['d']}  files={r['n']}  size={r['size']/1e6:.2f}MB")
    t = sorted(r['t1'] for r in rr if r['t1'])
    if len(t) > 1:
        gaps = [(t[i+1] - t[i]) / 3600 for i in range(len(t) - 1)]
        big = [i for i, g in enumerate(gaps) if g > 1]
        segs, start = [], 0
        for i in big + [len(t) - 1]:
            segs.append((start, i)); start = i + 1
        print(f"  完成时间 {fmt(t[0])} → {fmt(t[-1])}，分 {len(segs)} 段（间隔>1小时即断段）")
        for a, b in segs:
            n = b - a + 1
            if n > 5:
                per4 = 4 * (t[b] - t[a]) / (n - 1) / 60
                print(f"    段 {fmt(t[a])}–{fmt(t[b])}: {n} 条，按 4 并行折算单条约 {per4:.1f} 分钟")
    spans = [(r['t1'] - r['t0']) / 60 for r in rr if r['t0'] and r['t1'] and r['t1'] > r['t0']]
    if spans:
        print(f"  单目录输出文件首末修改时间差（粗略单条耗时）中位数 {statistics.median(spans):.1f} 分钟")
PY
} 2>&1 | redact > "$OUT"

echo "已写入 $OUT（$(wc -l < "$OUT") 行，$(du -h "$OUT" | cut -f1)）"
