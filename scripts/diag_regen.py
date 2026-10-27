#!/usr/bin/env python3
"""diag_regen.py - read-only diagnosis of the stop_regen.py restarts. Changes nothing; takes a few seconds.

A  Which year's climate do the restarts get? Every monthly output that exists in both the restart
   (runs_regen/<run_id>_same) and the uninterrupted run (runs_main/<run_id>) is compared at each offset L:
   restart year y against main-run year y - L. mpet.out depends on the climate only, so for it the true offset gives
   an error of ~0. Independent check: correlation with the Jul-Sep temperature of the climate files
   (export/clim_annual.csv). Also: save run vs main run up to 2024, and the years each output covers.
B  LPJ-GUESS source: the run command (run_main.run), restart handling in framework.cpp and cfinput.cpp,
   Standpft::serialize (are the per-stand establishment flags part of the saved state?), every use of reestab and
   where the functions that set it are called.
C  Instruction files of one STOP run: stand and management blocks, the 2025 switch, PFTs; regen vs same diff.
   (printed before B, which also looks up the switch parameters found in C)
Usage:  python diag_regen.py      (the paste block also writes the output to regen/diag_regen.log)
"""
import csv
import difflib
import inspect
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_main as rm  # noqa: E402

BASE = Path(rm.BASE)
MAIN, REGEN = BASE / "runs_main", BASE / "runs_regen"
START, END, FIRST_CLIM = 2025, 2100, 1901
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
SKIP_DIRS = {"output", "state", "node_modules", "__pycache__", "lpjg_venv", "venv", "export", "data", "miniconda3",
             "anaconda3", "miniforge3", "mambaforge", "R", "snap", "go", "site-packages", "dist-packages", "proc", "sys"}
SUMMARY = []
SWITCH_PARAMS = set()


def head(t):
    print(f"\n==== {t} " + "=" * max(4, 100 - len(t)), flush=True)


def stop_rows():
    rows = list(csv.DictReader(open(BASE / "design.csv")))
    return [r for r in rows if r["pulse"] == "none" and r["mgmt"] == "stop" and r["pstemp_high"] == "25"]


def pick(rows, gcm=None, ssp=None):
    """Prefer 500 m, phase 0."""
    c = [r for r in rows if (gcm is None or r["gcm"] == gcm) and (ssp is None or r["ssp"] == ssp)]
    best = [r for r in c if str(r["elev_m"]) == "500" and str(r["phase"]).lstrip("p") == "0"]
    return (best or c or [None])[0]


# ---------------------------------------------------------------- A: climate year after the restart
def read_out(d, f):
    import pandas as pd
    p = d / "output" / f
    if not p.is_file():
        return None
    t = pd.read_csv(p, sep=r"\s+")
    if "Year" not in t.columns:
        return None
    t = t.drop(columns=[c for c in ("Lon", "Lat") if c in t.columns])
    return t.groupby("Year").mean(numeric_only=True)


def rel_err(a, b):
    import numpy as np
    den = float(np.abs(b).mean())
    return float(np.abs(a - b).mean()) / den if den > 0 else float("nan")


def lag_scan(R, M, lags=range(-20, 171)):
    """{L: relative error of restart year y against main-run year y - L, y in START..END (>= 30 years overlap)}"""
    out = {}
    ry = [y for y in R.index if START <= y <= END]
    for L in lags:
        ys = [y for y in ry if (y - L) in M.index]
        if len(ys) >= 30:
            out[L] = rel_err(R.loc[ys].to_numpy(float), M.loc[[y - L for y in ys]].to_numpy(float))
    return out


def self30(X):
    ys = [y for y in X.index if START + 30 <= y <= END and (y - 30) in X.index]
    if len(ys) < 10:
        return float("nan")
    return rel_err(X.loc[ys].to_numpy(float), X.loc[[y - 30 for y in ys]].to_numpy(float))


def best_match(R, M, y):
    import numpy as np
    a = R.loc[y].to_numpy(float)
    e = np.abs(M.to_numpy(float) - a).mean(axis=1)
    j = int(e.argmin())
    return int(M.index[j]), rel_err(a, M.iloc[j].to_numpy(float))


def corr_lag(series, clim, lags=range(-20, 171)):
    """Lag with the largest |r| between a restart series (START..END) and a climate series (by year)."""
    import numpy as np
    best = (None, 0.0, 0)
    for L in lags:
        ys = [y for y in series.index if START <= y <= END and (y - L) in clim.index]
        if len(ys) < 30:
            continue
        a = np.diff(series.loc[ys].to_numpy(float))        # year-to-year changes: removes slow stand trends
        b = np.diff(clim.loc[[y - L for y in ys]].to_numpy(float))
        if a.std() == 0 or b.std() == 0:
            continue
        r = float(np.corrcoef(a, b)[0, 1])
        if abs(r) > abs(best[1]):
            best = (L, r, len(ys))
    return best


def part_a(rows):
    import pandas as pd
    head("A  which year's climate do the restarted runs get?")
    clim = None
    cp = BASE / "export" / "clim_annual.csv"
    if cp.is_file():
        c = pd.read_csv(cp)
        if {"gcm", "ssp", "elev_m", "Year", "tas_jas"} <= set(c.columns):
            clim = c
    mpet_lags = []
    for gcm in sorted({r["gcm"] for r in rows}):
        for ssp in sorted({r["ssp"] for r in rows}):
            r = pick(rows, gcm, ssp)
            if r is None:
                continue
            rid = r["run_id"]
            dm, dv, ds, dr = MAIN / rid, REGEN / f"{rid}_save", REGEN / f"{rid}_same", REGEN / f"{rid}_regen"
            print(f"\n-- {rid}")
            cov = []
            for lab, d in (("main", dm), ("save", dv), ("same", ds), ("regen", dr)):
                x = read_out(d, "cpool.out")
                cov.append(f"{lab} {int(x.index.min())}-{int(x.index.max())}" if x is not None else f"{lab} -")
            print("   cpool.out years: " + ", ".join(cov))
            m, v, s = read_out(dm, "cpool.out"), read_out(dv, "cpool.out"), read_out(ds, "cpool.out")
            col = "Total-slow" if m is not None and "Total-slow" in m.columns else None
            if col and v is not None:
                ys = [y for y in v.index if y in m.index and y < START]
                dmax = float((v.loc[ys, col] - m.loc[ys, col]).abs().max()) if ys else float("nan")
                print(f"   save vs main up to {START - 1}: max |EcoC diff| {dmax:.6f} kg C m-2 ({len(ys)} years)")
                SUMMARY.append(("save_vs_main", dmax))
            if col and s is not None and START in s.index:
                print(f"   EcoC (Total-slow): main {START - 1} {m.loc[START - 1, col]:.3f}, main {START} "
                      f"{m.loc[START, col]:.3f}, same {START} {s.loc[START, col]:.3f}; "
                      f"in {END}: main {m.loc[END, col]:.3f}, same {s.loc[END, col]:.3f}")
            files = sorted(p.name for p in (ds / "output").glob("m*.out")) if (ds / "output").is_dir() else []
            files = sorted(files, key=lambda f: (f != "mpet.out", f))[:10]
            for f in files:
                R, M = read_out(ds, f), read_out(dm, f)
                if R is None or M is None or not set(MONTHS) <= set(R.columns) or not set(MONTHS) <= set(M.columns):
                    continue
                R, M = R[MONTHS], M[MONTHS]
                sc = lag_scan(R, M)
                if not sc:
                    print(f"   {f:18s} no overlap for the lag scan (main years {int(M.index.min())}-"
                          f"{int(M.index.max())})")
                    continue
                Lb = min(sc, key=sc.get)
                SUMMARY.append(("lag", (f, Lb, sc[Lb])))
                k = START - FIRST_CLIM
                print(f"   {f:18s} best offset L={Lb:4d} (restart {START} ~ main {START - Lb}) rel.err {sc[Lb]:.4f};"
                      f" at L=0 {sc.get(0, float('nan')):.4f}; at L={k} {sc.get(k, float('nan')):.4f};"
                      f" 30-yr repeat: restart {self30(R):.4f}, main {self30(M):.4f}")
                if f == "mpet.out":
                    mpet_lags.append((rid, Lb, sc[Lb], sc.get(0, float("nan"))))
                    pairs = [best_match(R, M, y) + (y,) for y in (START, START + 1, START + 29, START + 30, END)
                             if y in R.index]
                    print("   mpet.out, best-matching main-run year: "
                          + "; ".join(f"{y} -> {my} ({e:.4f})" for my, e, y in pairs))
            if clim is not None:
                cc = clim[(clim.gcm == gcm) & (clim.ssp == ssp) & (clim.elev_m.astype(str) == str(r["elev_m"]))]
                cc = cc.groupby("Year").tas_jas.mean()
                for f, cols in (("mpet.out", ["Jul", "Aug", "Sep"]), ("mnpp.out", ["Jul", "Aug", "Sep"])):
                    R = read_out(ds, f)
                    if R is None or not set(cols) <= set(R.columns) or cc.empty:
                        continue
                    L, rr, n = corr_lag(R[cols].sum(axis=1), cc)
                    L0 = corr_lag(R[cols].sum(axis=1), cc, lags=[0])
                    print(f"   climate files: {f} Jul-Sep vs tas_jas, strongest at L={L} (r={rr:+.2f}, n={n});"
                          f" at L=0 r={L0[1]:+.2f}")
    if mpet_lags:
        Ls = sorted({L for _, L, _, _ in mpet_lags})
        err = max(e for _, _, e, _ in mpet_lags)
        e0 = min(e for _, _, _, e in mpet_lags)
        SUMMARY.append(("mpet", (Ls, err, e0, len(mpet_lags))))
    else:
        SUMMARY.append(("mpet", None))


# ---------------------------------------------------------------- B: model source
def walk_find(top, name, maxdepth=6, limit=5, budget=25.0):
    import time
    top, out, t0 = Path(top), [], time.time()
    if not top.is_dir() or top == Path("/"):
        return out
    for dp, dns, fns in os.walk(top):
        if time.time() - t0 > budget:
            print(f"   (search for {name} under {top} stopped after {budget:.0f} s)")
            break
        depth = len(Path(dp).relative_to(top).parts)
        dns[:] = [d for d in dns if d not in SKIP_DIRS and not d.startswith(("runs_", ".")) and depth < maxdepth]
        if name in fns:
            out.append(Path(dp) / name)
            if len(out) >= limit:
                break
    return out


def grep(path, pat, ctx=0, cap=80, label=None):
    lines = path.read_text(errors="ignore").splitlines()
    rx = re.compile(pat)
    keep = sorted({j for i, ln in enumerate(lines) if rx.search(ln)
                   for j in range(max(0, i - ctx), min(len(lines), i + ctx + 1))})
    print(f"-- {label or path.name}: /{pat}/ ({len(keep)} lines)")
    prev = None
    for n, j in enumerate(keep):
        if n >= cap:
            print(f"   [cut at {cap} lines]")
            break
        if prev is not None and j != prev + 1:
            print("   ...")
        print(f"{j + 1:6d}: {lines[j]}")
        prev = j
    return lines


def body(path, sig, cap=260):
    """Print the definition of the function whose header matches sig (brace counting)."""
    lines = path.read_text(errors="ignore").splitlines()
    rx = re.compile(sig)
    for i, ln in enumerate(lines):
        if rx.search(ln) and not ln.strip().startswith(("//", "*")) and not ln.rstrip().endswith(";"):
            print(f"-- {path.name}: definition of /{sig}/")
            depth, started = 0, False
            for j in range(i, min(len(lines), i + cap)):
                print(f"{j + 1:6d}: {lines[j]}")
                code = lines[j].split("//")[0]
                depth += code.count("{") - code.count("}")
                started = started or "{" in code
                if started and depth <= 0:
                    return "\n".join(lines[i:j + 1])
            print(f"   [cut at {cap} lines]")
            return "\n".join(lines[i:i + cap])
    print(f"-- {path.name}: no definition matching /{sig}/")
    return ""


FUNC_RX = re.compile(r"^[A-Za-z_][\w:<>,\*&\s]*?\b([A-Za-z_~][\w~]*(?:::[A-Za-z_~][\w~]*)*)\s*\(")


def enclosing(lines, i):
    for j in range(i, -1, -1):
        ln = lines[j]
        if ln[:1] and not ln[:1].isspace() and ln[:1] not in "{}#/*" and "(" in ln and not ln.rstrip().endswith(";"):
            m = FUNC_RX.match(ln)
            if m and m.group(1) not in ("if", "for", "while", "switch", "return"):
                return m.group(1)
    return "?"


def part_b():
    head("B  LPJ-GUESS source: restart, climate input, establishment flags")
    bins = []
    for k, v in sorted(vars(rm).items()):
        if isinstance(v, (str, Path)) and "guess" in str(v).lower() and "/" in str(v):
            p = Path(str(v)).expanduser()
            bins.append(p)
            print(f"run_main.{k} = {p}  ({'file' if p.is_file() else 'not a file'})")
    w = shutil.which("guess")
    if w:
        bins.append(Path(w))
        print(f"guess on PATH: {w}")
    if hasattr(rm, "run"):
        try:
            src = inspect.getsource(rm.run).splitlines()
            print("-- run_main.run:")
            print("\n".join(src[:60]) + ("\n   [cut]" if len(src) > 60 else ""))
        except (OSError, TypeError) as e:
            print(f"run_main.run: source not available ({e})")
    hits = []
    for b in bins:
        if b.exists():
            ups = ([b.resolve()] if b.is_dir() else []) + list(b.resolve().parents)[:4]
            for up in ups:
                hits = walk_find(up, "cfinput.cpp", maxdepth=3, limit=3)
                if hits:
                    break
        if hits:
            break
    if not hits:
        hits = walk_find(Path.home(), "cfinput.cpp", maxdepth=6)
    for extra in (Path("/opt"), Path("/usr/local/src")):
        if not hits:
            hits = walk_find(extra, "cfinput.cpp", maxdepth=4)
    print("cfinput.cpp found: " + (", ".join(map(str, hits)) if hits else "none"))
    if not hits:
        SUMMARY.append(("src", None))
        return
    cf = hits[0]
    root = cf.parent.parent if cf.parent.name == "modules" else cf.parent
    git = ""
    if (root / ".git").exists():
        try:
            git = subprocess.run(["git", "-C", str(root), "log", "-1", "--format=%h %ci %s"], capture_output=True,
                                 text=True, timeout=20).stdout.strip()
            dirty = subprocess.run(["git", "-C", str(root), "status", "--short", "--untracked-files=no"],
                                   capture_output=True, text=True, timeout=20).stdout.strip().splitlines()
            git += f"; locally modified files: {len(dirty)} {[d.split()[-1] for d in dirty[:8]]}"
        except Exception as e:  # noqa: BLE001
            git = f"git not readable ({e})"
    print(f"source root: {root}  {git}")
    SUMMARY.append(("src", f"{root} {git.split(';')[0]}"))

    def one(name):
        f = walk_find(root, name, maxdepth=3, limit=2)
        return f[0] if f else None

    fw, gc, gh = one("framework.cpp"), one("guess.cpp"), one("guess.h")
    if fw:
        grep(fw, r"\brestart\b|state_year|save_state|state_path", ctx=3, cap=110)
    grep(cf, r"restart|state_year|nyear_spinup|calendar_year|first_calendar|current_day|timestep|historic",
         cap=80, label="cfinput.cpp")
    SUMMARY.append(("cf_restart", bool(re.search(r"restart|state_year", cf.read_text(errors="ignore")))))
    gb = body(cf, r"CFInput::getclimate\s*\(", cap=280)
    idx = sorted({v for v in re.findall(r"\[\s*([A-Za-z_]\w*)", gb)
                  if v not in ("date", "day", "i", "j", "m", "d", "k", "n", "pft", "id")})
    if idx:
        grep(cf, r"\b(" + "|".join(map(re.escape, idx)) + r")\b", cap=60, label=f"cfinput.cpp, index variables {idx}")
    called = [f for f in dict.fromkeys(re.findall(r"\b([A-Za-z_]\w*)\s*\(", gb))
              if re.search(r"CFInput::" + f + r"\s*\(", cf.read_text(errors="ignore")) and f != "getclimate"]
    for f in called[:3]:
        body(cf, r"CFInput::" + f + r"\s*\(", cap=130)
    body(cf, r"CFInput::getgridcell\s*\(", cap=130)
    if gh:
        grep(gh, r"\bbool\s+(active|plant|reestab|irrigated)\b", ctx=1, cap=30, label="guess.h")
    ser = body(gc, r"Standpft::serialize\s*\(", cap=60) if gc else ""
    SUMMARY.append(("ser_reestab", ("reestab" in ser) if ser else None))

    srcs = [p for p in root.rglob("*") if p.suffix in (".cpp", ".h") and p.is_file()
            and not any(part.startswith(("build", ".git")) for part in p.relative_to(root).parts)]
    print(f"\n-- every use of reestab in {len(srcs)} source files (with the enclosing function)")
    funcs, n = {}, 0
    for p in sorted(srcs):
        lines = p.read_text(errors="ignore").splitlines()
        for i, ln in enumerate(lines):
            if "reestab" in ln:
                fn = enclosing(lines, i)
                blk = next((re.search(r"case\s+(BLOCK_\w+)", lines[j]).group(1) for j in range(i, max(-1, i - 400), -1)
                            if re.search(r"case\s+BLOCK_\w+", lines[j])), None)
                fn = fn + (f", {blk}" if blk else "")
                if n < 90:
                    print(f"   {p.relative_to(root)}:{i + 1}: [{fn}] {ln.strip()[:150]}")
                n += 1
                if re.search(r"reestab\s*=[^=]", ln):
                    funcs.setdefault(fn.split(",")[0].split("::")[-1], set()).add(f"{p.relative_to(root)}:{i + 1}")
    if n > 90:
        print(f"   [{n - 90} more not shown]")
    print(f"\n-- functions that assign reestab: {sorted(funcs)}; where they are called:")
    shown = 0
    for fn in sorted(funcs):
        if fn in ("?", "serialize", "operator"):
            continue
        rx = re.compile(r"\b" + re.escape(fn) + r"\s*\(")
        for p in sorted(srcs):
            lines = p.read_text(errors="ignore").splitlines()
            for i, ln in enumerate(lines):
                if rx.search(ln) and shown < 60 and not (ln[:1].strip() and ln.rstrip().endswith(("{", ")"))):
                    print(f"   {fn}: {p.relative_to(root)}:{i + 1}: [{enclosing(lines, i)}] {ln.strip()[:140]}")
                    shown += 1
    names = sorted({"naturalveg", "restrictpfts"} | SWITCH_PARAMS)
    print(f"\n-- naturalveg, restrictpfts and the switch parameters of part C in the source: {names}")
    for nm in names:
        rx, n = re.compile(r"\b" + re.escape(nm) + r"\b"), 0
        for p in sorted(srcs):
            lines = p.read_text(errors="ignore").splitlines()
            for i, ln in enumerate(lines):
                if rx.search(ln) and n < 12:
                    print(f"   {nm}: {p.relative_to(root)}:{i + 1}: [{enclosing(lines, i)}] {ln.strip()[:140]}")
                    n += 1


# ---------------------------------------------------------------- C: instruction files
KEYS = r"reestab|naturalveg|restrictpfts|plant|harvest|manage|rotation|cut|thin|clear|landcover|stinclude|selection"


def ins_blocks(lines, kinds=("st", "mt", "group", "pft")):
    rx = re.compile(r'^\s*(%s)\s+"([^"]+)"\s*\(' % "|".join(kinds))
    out, i = [], 0
    while i < len(lines):
        m = rx.match(lines[i])
        if not m:
            i += 1
            continue
        depth, j = 0, i
        while j < len(lines):
            code = lines[j].split("!")[0]
            depth += code.count("(") - code.count(")")
            if depth <= 0:
                break
            j += 1
        out.append((m.group(1), m.group(2), i, min(j, len(lines) - 1)))
        i = j + 1
    return out


def part_c(rows):
    head("C  instruction files of one STOP run")
    r = pick(rows, "gfdl-esm4", "ssp126") or pick(rows)
    rid = r["run_id"]
    d = MAIN / rid
    ins = sorted(d.glob("*.ins"))
    print(f"{d}: " + ", ".join(f"{p.name} ({len(p.read_text(errors='ignore').splitlines())} lines)" for p in ins))
    inside = {}
    pfts, printed = [], 0
    for p in ins:
        lines = p.read_text(errors="ignore").splitlines()
        for kind, name, a, b in ins_blocks(lines):
            for k in range(a, b + 1):
                inside.setdefault(p.name, set()).add(k)
            text = "\n".join(lines[a:b + 1])
            if kind == "pft":
                inc = re.findall(r"^\s*include\s+(\d)", text, flags=re.M)
                pfts.append(f"{name}: include {','.join(inc) or '(not set in block)'} ({p.name}:{a + 1})")
            elif kind in ("st", "mt"):
                for ln in lines[a + 1:b]:
                    w = ln.split("!")[0].split()
                    if w and re.search(r"manage|rot|year|cycle|switch|first|last|stop", w[0], flags=re.I):
                        SWITCH_PARAMS.add(w[0])
            if kind != "pft" and re.search(KEYS, text) and printed < 420:
                print(f"\n-- {p.name}:{a + 1}-{b + 1}  {kind} \"{name}\"")
                for k in range(a, min(b + 1, a + 90)):
                    print(f"{k + 1:6d}: {lines[k]}")
                    printed += 1
    print("\n-- pft blocks: " + ("; ".join(pfts) if pfts else "none found"))
    pat = (r"nyear_spinup|state_|restart|save_state|^\s*import|run_landcover|lcfrac|file_lu|stfrac|landcover|"
           r"manage|rotation|\b20[0-9]{2}\b|first|last|ifdisturb|distinterval|npatch|estab|reestab|naturalveg")
    print(f"\n-- other lines outside the blocks: /{pat}/")
    n = 0
    for p in ins:
        for k, ln in enumerate(p.read_text(errors="ignore").splitlines()):
            if k not in inside.get(p.name, set()) and re.search(pat, ln) and n < 140:
                print(f"   {p.name}:{k + 1}: {ln.rstrip()[:160]}")
                n += 1
                code = ln.split("!")[0]
                if re.search(r"\b20[0-9]{2}\b|manage|rot|switch|stop", code, flags=re.I):
                    m = re.match(r'\s*param\s+"(\w+)"', code) or re.match(r"\s*([A-Za-z_]\w*)\s", code)
                    if m and m.group(1) not in ("import", "param"):
                        SWITCH_PARAMS.add(m.group(1))
    cfi = d / "lishui_cf.ins"
    if cfi.is_file():
        lines = cfi.read_text(errors="ignore").splitlines()
        print(f"\n-- lishui_cf.ins ({len(lines)} lines, first 80)")
        for k, ln in enumerate(lines[:80]):
            print(f"{k + 1:6d}: {ln}")
    for lab, dd in (("main", d), ("regen", REGEN / f"{rid}_regen")):
        p = dd / "output" / "cmass.out"
        if p.is_file():
            print(f"\n-- cmass.out columns ({lab}): {open(p, errors='ignore').readline().split()}")
    a, b = REGEN / f"{rid}_same" / "lishui_fir.ins", REGEN / f"{rid}_regen" / "lishui_fir.ins"
    if a.is_file() and b.is_file():
        diff = list(difflib.unified_diff(a.read_text(errors="ignore").splitlines(), b.read_text(errors="ignore").splitlines(), "same", "regen",
                                         n=4, lineterm=""))
        print(f"\n-- lishui_fir.ins, same vs regen restart ({len(diff)} diff lines)")
        print("\n".join(diff[:50]))
        tail = a.read_text(errors="ignore").splitlines()[-10:]
        print("\n-- last lines of the same restart's lishui_fir.ins")
        print("\n".join(tail))
    rl = REGEN / f"{rid}_regen" / "run.log"
    if rl.is_file():
        t = rl.read_text(errors="ignore").splitlines()
        print(f"\n-- {rl.relative_to(BASE)}: first 15 and last 10 of {len(t)} lines")
        print("\n".join(t[:15] + ["   ..."] + t[-10:]))


def summary():
    head("SUMMARY")
    s = dict(SUMMARY)
    m = s.get("mpet")
    if m:
        Ls, err, e0, n = m
        yrs = ", ".join(str(START - L) for L in Ls)
        print(f"A  mpet.out (climate only): restart year {START} matches main-run year {yrs} "
              f"(offset {Ls}, {n} runs, rel. error <= {err:.4f}; at offset 0 >= {e0:.4f})")
    else:
        print("A  mpet.out not found in both runs: see the other monthly files and the tas_jas check in part A")
    by = {}
    for k, v in SUMMARY:
        if k == "lag":
            by.setdefault(v[0], []).append(v[1:])
    if by:
        print("A  best offset by monthly file (restart year y ~ main-run year y - L): " + "; ".join(
            f"{f} L={sorted({L for L, _ in x})} in {len(x)} runs, rel.err <= {max(e for _, e in x):.4f}"
            for f, x in sorted(by.items())))
    if "save_vs_main" in s:
        v = max(x for k, x in SUMMARY if k == "save_vs_main")
        print(f"A  save runs vs main runs up to {START - 1}: max |EcoC diff| {v:.6f} kg C m-2")
    print(f"B  LPJ-GUESS source: {s.get('src') or 'not found'}")
    if s.get("src"):
        print(f"B  cfinput.cpp mentions restart/state_year: {'yes' if s.get('cf_restart') else 'no'}")
        r = s.get("ser_reestab")
        print("B  Standpft::serialize saves reestab: " + {True: "yes", False: "no", None: "definition not found"}[r])
    print("full output: ~/lpjg_lishui/regen/diag_regen.log  -> please upload this file")


if __name__ == "__main__":
    rows = stop_rows()
    print(f"diag_regen.py  base {BASE}; {len(rows)} STOP rows (pstemp_high 25)")
    for part in (lambda: part_a(rows), lambda: part_c(rows), part_b):
        try:
            part()
        except Exception as e:  # noqa: BLE001  one failing part must not hide the others
            import traceback
            print(f"\n[this part failed: {e!r}]")
            traceback.print_exc(file=sys.stdout)
    summary()
