#!/usr/bin/env python3
"""calib_growth.py - calibrate one TeNE parameter against Lishui/Zhejiang Chinese fir carbon-by-age targets.

Usage:  python calib_growth.py [PARAM] [V1 V2 ...] [NAME=VALUE ...]    default: respcoeff 0.3 0.5 0.7 0.9 1.2
        NAME=VALUE fixes another TeNE parameter in all runs of the sweep (written before PARAM), e.g.
        python calib_growth.py respcoeff 0.40 0.45 0.50 0.55 0.60 0.70 pstemp_high=30
Template run: rotland/r26_o0_ctrl (26-yr rotation from 1901, 2500 stems/ha, thinning at age 9 and 13,
root_beta 0.96, observed climate). Targets: docs/calib_targets.csv (Diao 2022, Li 2016).
N deposition: every run gets ndep_timeseries "rcp60" (= the pre-patch behaviour of the original calibration).
Output: calib/<TAG>/<PARAM>_<value>/ (runs, each with run_info.txt) and docs/calib_<TAG>_summary.csv,
        TAG = PARAM plus the fixed parameters (e.g. respcoeff_pstemp_high30); without them TAG = PARAM.
"""
import datetime
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path.home() / "lpjg_lishui"
TEMPLATE = BASE / "rotland" / "r26_o0_ctrl"
GUESS = (BASE / ".guess2_bin").read_text().strip()
ROT, FIRST, LAST_YEAR = 26, 1901, 2019
GROUPS = {"young": (1, 10), "middle": (11, 20), "near_mature": (21, 25)}   # stand-age windows (LY/T 2908)


def targets():
    """Mean and range of the two sources per age group (kg C m-2), and soil+litter per source."""
    t = pd.read_csv(BASE / "docs" / "calib_targets.csv")
    tree = t[t.variable.isin(["living_biomass_C", "tree_layer_C_planted"])]
    groups = {g: (tree[tree.age_group == g].value.mean() / 10, tree[tree.age_group == g].value.min() / 10,
                  tree[tree.age_group == g].value.max() / 10) for g in GROUPS}
    soil = {s: t[(t.source == s) & t.variable.str.startswith(("soil_C", "litter_C"))].value.sum() / 10
            for s in ("Diao2022", "Li2016")}
    return groups, soil


def make_run(param, value, fixed, tag):
    d = BASE / "calib" / tag / f"{param}_{value}"
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(TEMPLATE, d, ignore=shutil.ignore_patterns("output", "run.log", "state"))
    (d / "output").mkdir()
    g = d / "lishui_global.ins"
    out, inb = [], False
    for ln in g.read_text().splitlines():
        s = ln.strip()
        if s.startswith("pft") and '"TeNE"' in s:
            inb = True
        if inb and s.startswith(")"):
            out += [f"\t{k} {v}   ! fixed in this sweep (calib_growth.py)" for k, v in fixed]
            out.append(f"\t{param} {value}   ! growth calibration (calib_growth.py)")
            inb = False
        out.append(ln)
    g.write_text("\n".join(out) + "\n")
    cf = d / "lishui_cf.ins"
    cf.write_text("\n".join(ln for ln in cf.read_text().splitlines() if not ln.strip().startswith("ndep_timeseries"))
                  + '\nndep_timeseries "rcp60"\n')
    commit = subprocess.run(["git", "-C", str(BASE), "rev-parse", "--short", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    (d / "run_info.txt").write_text(f"param={param}\nvalue={value}\nfixed={' '.join(f'{k}={v}' for k, v in fixed)}\n"
                                    f"template={TEMPLATE}\ngit={commit}\n"
                                    f"created={datetime.datetime.now().isoformat(timespec='seconds')}\n")
    return d


def run(d):
    with open(d / "run.log", "w") as log:
        subprocess.run([GUESS, "-input", "cf", "lishui_fir.ins"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    return "Finished" in (d / "run.log").read_text(errors="ignore")


def series(d, f, col):
    p = d / "output" / f
    if not p.exists():
        return pd.Series(dtype=float)
    df = pd.read_csv(p, sep=r"\s+")
    if col not in df.columns:
        return pd.Series(0.0, index=sorted(df.Year.unique()))
    return df.groupby("Year")[col].mean()


def evaluate(d):
    fir = series(d, "cmass.out", "TeNE")                     # tree (fir) carbon only, excludes understory grass
    res = {}
    starts = range(FIRST + ROT, LAST_YEAR - ROT + 2, ROT)   # complete rotations after the first one
    for g, (a0, a1) in GROUPS.items():
        res[g] = float(np.nanmean([fir.reindex([s + a for a in range(a0, a1 + 1)]).mean() for s in starts]))
    tot, veg = series(d, "cpool.out", "Total"), series(d, "cpool.out", "VegC")
    prod = series(d, "cpool.out", "HarvSlowC").reindex(tot.index).fillna(0)
    res["soil_litter"] = float((tot - veg - prod).loc[1990:2019].mean())
    npp, gpp = series(d, "anpp.out", "Total"), series(d, "agpp.out", "Total")
    res["NPP"] = float(npp.loc[1990:2019].mean())
    res["NPP_GPP"] = float((npp / gpp).loc[1990:2019].mean()) if len(gpp) else np.nan
    res["harvest_yr"] = float(series(d, "cflux.out", "Harvest").loc[1960:2019].mean())
    return res


if __name__ == "__main__":
    fixed = [tuple(a.split("=", 1)) for a in sys.argv[1:] if "=" in a]
    plain = [a for a in sys.argv[1:] if "=" not in a]
    param = plain[0] if plain else "respcoeff"
    values = plain[1:] or ["0.3", "0.5", "0.7", "0.9", "1.2"]
    tag = param + "".join(f"_{k}{v}" for k, v in fixed)
    if fixed:
        print("Fixed in all runs: " + ", ".join(f"{k} {v}" for k, v in fixed))
    tg, soil = targets()
    print("Targets (kg C m-2): " + ", ".join(f"{g} {m:.2f} ({lo:.2f}-{hi:.2f})" for g, (m, lo, hi) in tg.items())
          + f"; soil+litter {soil['Diao2022']:.2f} (Diao) to {soil['Li2016']:.2f} (Li)")
    dirs = [make_run(param, v, fixed, tag) for v in values]
    with ThreadPoolExecutor(max_workers=4) as ex:
        ok = list(ex.map(run, dirs))
    rows = []
    for v, d, good in zip(values, dirs, ok):
        if not good:
            print(f"{param}={v}: run failed, see {d / 'run.log'}")
            continue
        r = {param: float(v), **evaluate(d)}
        r["rmse"] = float(np.sqrt(np.mean([(r[g] - tg[g][0]) ** 2 for g in GROUPS])))
        rows.append(r)
    if not rows:
        sys.exit("no successful runs")
    df = pd.DataFrame(rows).set_index(param)
    out = BASE / "docs" / f"calib_{tag}_summary.csv"
    df.round(4).to_csv(out)
    pd.set_option("display.width", 200)
    print(df.round(3).to_string())
    best = df.rmse.idxmin()
    print(f"\nBest {param} = {best} (RMSE {df.rmse.min():.2f} kg C m-2, NPP/GPP {df.loc[best, 'NPP_GPP']:.2f}); "
          f"summary written to {out}")
