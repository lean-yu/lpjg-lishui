#!/usr/bin/env python3
"""stop_regen.py - sensitivity of the STOP scenario to natural regeneration.

In the main runs a stopped stand keeps the managed stand type's reestab "none": after 2024 it only ages, nothing
establishes. Here the STOP runs (pstemp_high 25; all bands, GCMs, SSPs and phases = 48) are restarted from a model
state saved at the start of 2025, with reestab "all" (natural regeneration of every PFT), as in checks/rs_stop.
Every restart has a twin that changes nothing ("same"), restarted from the same state, so
  regeneration effect = regen - same
and "same" is compared with the uninterrupted run (restart identity check). Both are done with the SSP's rising CO2
and with CO2 held at its 2024 value (the files of co2fix.py).
Runs (runs_regen/): <run_id>_save (spin-up to 2024, saves the state); <run_id>_same, _regen, _same_fix, _regen_fix
(2025-2100, restarted). About 15-20 min with 4 workers.
Output (regen/): pairs.csv (one row per run and CO2 variant), summary.csv, annual_regen.csv, fig_regen.png.
Usage:  python stop_regen.py               run what is not finished yet, then analyse
        python stop_regen.py --analyze     analyse only
"""
import argparse
import csv
import datetime
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_main as rm  # noqa: E402  same build rules as the main runs
import co2fix  # noqa: E402  fixed-CO2 files

BASE = rm.BASE
MAIN, FIXD = BASE / "runs_main", BASE / "runs_co2fix"
OUT = BASE / "runs_regen"
RES = BASE / "regen"
START = 2025                                    # first restarted year
VARIANTS = {"same": (False, False), "regen": (True, False), "same_fix": (False, True), "regen_fix": (True, True)}
W = 0.38                                        # commercial share, as in analyze_pulse.py


def done(d):
    log = d / "run.log"
    return log.exists() and "Finished" in log.read_text(errors="ignore")


def selected_rows():
    rows = list(csv.DictReader(open(BASE / "design.csv")))
    return [r for r in rows if r["pulse"] == "none" and r["mgmt"] == "stop" and r["pstemp_high"] == "25"]


def nyear_spinup(d):
    """Last nyear_spinup in import order (global, cf, fir), as check_calibrated.sh."""
    v = None
    for f in ("lishui_global.ins", "lishui_cf.ins", "lishui_fir.ins"):
        for ln in (d / f).read_text().splitlines():
            t = ln.split()
            if len(t) >= 2 and t[0] == "nyear_spinup":
                v = int(t[1])
    if v is None:
        raise SystemExit(f"nyear_spinup not found in {d}")
    return v


def add_reestab_all(fir):
    """reestab "all" at the end of the st "Fir_CC" block; overrides the managed group's reestab "none"."""
    out, inb = [], False
    for ln in fir.splitlines():
        if re.match(r'\s*st\s*"Fir_CC"', ln):
            inb = True
        if inb and re.match(r"\s*\)", ln):
            out.append('\treestab "all"   ! natural regeneration after harvest stops (stop_regen.py)')
            inb = False
        out.append(ln)
    if inb or not any("stop_regen.py" in ln for ln in out):
        raise SystemExit('st "Fir_CC" block not found')
    return "\n".join(out) + "\n"


def build_save(row, commit):
    """The STOP run up to 2024, saving the state at the start of START. None if finished earlier."""
    d = rm.build(dict(row, run_id=row["run_id"] + "_save"), commit)
    if d is None:
        return None
    sy = nyear_spinup(d) + START - 1901          # state_year counts years from the first spin-up year
    (d / "state").mkdir()
    with open(d / "lishui_fir.ins", "a") as fh:
        fh.write(f'\n! save the model state at the start of {START} (stop_regen.py)\nstate_path "{d}/state/"\n'
                 f"save_state 1\nstate_year {sy}\n")
    return d


def build_restart(save_dir, name, regen, fixed, co2_file):
    d = OUT / name
    if done(d):
        return None
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(save_dir, d, ignore=shutil.ignore_patterns("output", "state", "run.log", "guess.log"))
    (d / "output").mkdir()
    fir = re.sub(r"(?m)^(\s*save_state\s+)1", r"\g<1>0", (d / "lishui_fir.ins").read_text()) + "restart 1\n"
    (d / "lishui_fir.ins").write_text(add_reestab_all(fir) if regen else fir)
    if fixed:
        cf = d / "lishui_cf.ins"
        cf.write_text(rm.set_param(cf.read_text(), "file_co2", str(co2_file)))
    with open(d / "run_info.txt", "a") as fh:
        fh.write(f"variant={name}: restart from {save_dir.name} at {START}; reestab "
                 f"{'all' if regen else 'none (unchanged)'}; CO2 {'held at 2024' if fixed else 'SSP'}\n")
    return d


def stage(label, todo, workers):
    print(f"{label}: {len(todo)} to run, start {datetime.datetime.now():%H:%M:%S}", flush=True)
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(rm.run, d) for d in todo]
        for i, f in enumerate(as_completed(futs), 1):
            name, ok = f.result()
            print(f"{datetime.datetime.now():%H:%M:%S} {i}/{len(todo)} {name} {'ok' if ok else 'FAILED'}", flush=True)
            if not ok:
                failed.append(name)
    print(f"{label}: finished {datetime.datetime.now():%H:%M:%S}; failed: {len(failed)} {failed[:5]}", flush=True)


def run_all(rows, workers):
    co2 = {s: co2fix.make_co2(s)[0] for s in sorted({r["ssp"] for r in rows})}
    rm.OUT = OUT                                 # run_main.build puts each run under rm.OUT / run_id
    commit = rm.git_commit()
    todo = []
    for r in rows:
        d = build_save(r, commit)
        if d is not None:
            todo.append(d)
            for v in VARIANTS:                   # a new state invalidates the old restarts
                shutil.rmtree(OUT / f"{r['run_id']}_{v}", ignore_errors=True)
    stage(f"save runs (to {START - 1})", todo, workers)
    todo = []
    for r in rows:
        sd = OUT / f"{r['run_id']}_save"
        if not (done(sd) and (sd / "state").is_dir() and any((sd / "state").iterdir())):
            print(f"no state from {sd.name}; its restarts are skipped")
            continue
        for v, (regen, fixed) in VARIANTS.items():
            d = build_restart(sd, f"{r['run_id']}_{v}", regen, fixed, co2[r["ssp"]])
            if d is not None:
                todo.append(d)
    stage(f"restarts ({START}-2100)", todo, workers)


def analyse(rows):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def read(d, f):
        p = d / "output" / f
        return pd.read_csv(p, sep=r"\s+").groupby("Year").mean(numeric_only=True)

    def series(d, mass=False):
        n, cf, cp = read(d, "anpp.out"), read(d, "cflux.out"), read(d, "cpool.out")
        x = pd.DataFrame({"NPP": n["Total"], "NEP": -cf["NEE-LU"], "VegC": cp["VegC"], "EcoC": cp["Total-slow"]})
        if mass:
            cm = read(d, "cmass.out")
            grass = [c for c in ("C3G", "C4G") if c in cm.columns]
            x["firC"] = cm["TeNE"]
            x["grassC"] = cm[grass].sum(axis=1) if grass else 0.0
            x["otherwoodyC"] = cm["Total"] - x.firC - x.grassC
        return x

    RES.mkdir(exist_ok=True)
    keys = ["run_id", "elev_m", "gcm", "ssp", "phase"]
    out, ann, cum, missing = [], [], [], []
    for r in rows:
        rid = r["run_id"]
        dirs = {v: OUT / f"{rid}_{v}" for v in VARIANTS}
        if not (done(MAIN / rid) and all(done(d) for d in dirs.values())):
            missing.append(rid)
            continue
        base = series(MAIN / rid)                 # identical to the save run up to 2024
        e24 = base.EcoC[START - 1]
        x = {v: series(d, mass=True) for v, d in dirs.items()}
        ref = {"": base, "_fix": series(FIXD / rid) if done(FIXD / rid) else None}
        fac = {k: r[k] for k in keys}
        for suf, co2 in (("", "SSP"), ("_fix", "fixed 2024")):
            s, g = x["same" + suf], x["regen" + suf]
            row = dict(fac, co2=co2)
            rr = ref[suf]
            row["restart_maxdiff"] = (float((s.EcoC - rr.EcoC).loc[START:2100].abs().max())
                                      if rr is not None else float("nan"))
            for lab, y in (("same", s), ("regen", g)):
                row[f"NEP_{lab}"] = y.NEP.loc[START:2100].sum()
                row[f"NEP_late_{lab}"] = y.NEP.loc[2081:2100].mean()
                row[f"dEcoC_{lab}"] = y.EcoC[2100] - e24
                row[f"VegC2100_{lab}"] = y.VegC[2100]
            for c in ("firC", "otherwoodyC", "grassC"):
                row[f"{c}2100_regen"] = g[c][2100]
            row["dNEP_regen"] = row["NEP_regen"] - row["NEP_same"]
            out.append(row)
            for lab, y in (("same", s), ("regen", g)):
                ann.append(y.loc[START:2100].reset_index().assign(variant=lab, **fac, co2=co2))
                cum.append(pd.DataFrame({"Year": range(START, 2101), "cumNEP": y.NEP.loc[START:2100].cumsum()
                                         .to_numpy()}).assign(variant=lab, co2=co2, gcm=r["gcm"], ssp=r["ssp"]))
    if not out:
        raise SystemExit(f"no finished sets ({len(missing)} missing)")
    p = pd.DataFrame(out)
    p.round(5).to_csv(RES / "pairs.csv", index=False)
    pd.concat(ann, ignore_index=True).round(5).to_csv(RES / "annual_regen.csv", index=False)
    pd.set_option("display.width", 200)
    print(f"\n{len(p) // 2} stop runs analysed (x 2 CO2 variants), {len(missing)} missing")
    md = p.groupby("co2").restart_maxdiff.max()
    print("restart identity (max |same - uninterrupted run|, ecosystem C 2025-2100, kg C m-2): "
          + ", ".join(f"{k} {v:.4f}" for k, v in md.items())
          + ("  -> restart reproduces the main runs" if (md < 1e-3).all() else
             "  -> restart is not exact; use regen - same (both restarted) for the effect"))

    g = p.groupby(["co2", "gcm", "ssp"])
    t = pd.DataFrame({"NEP_same": g.NEP_same.mean(), "NEP_regen": g.NEP_regen.mean(), "effect": g.dNEP_regen.mean(),
                      "late_NEP_same": g.NEP_late_same.mean(), "late_NEP_regen": g.NEP_late_regen.mean()})
    print("\n== STOP stands: cumulative NEP 2025-2100 and mean NEP 2081-2100 (kg C m-2; per yr for late NEP), "
          "mean over bands and phases ==")
    print(t.round(3).to_string())
    t.reset_index().round(4).to_csv(RES / "summary.csv", index=False)

    c = p.groupby(["co2", "ssp"])[["firC2100_regen", "otherwoodyC2100_regen", "grassC2100_regen",
                                   "VegC2100_same", "VegC2100_regen"]].mean()
    print("\n== vegetation in 2100 (kg C m-2): fir, other woody PFTs and grasses with regeneration; "
          "total vegetation without (same) and with (regen) ==")
    print(c.round(2).to_string())

    e = p.groupby("co2")[["dEcoC_same", "dEcoC_regen", "dNEP_regen"]].mean()
    e["landscape_BAU_or_EXT_change"] = (1 - W) * e.dNEP_regen
    print(f"\n== effect on the landscapes (kg C m-2, 2025-2100): STOP changes by dNEP_regen, BAU and EXT by "
          f"(1 - W) x that, W = {W} ==")
    print(e.round(2).to_string())

    q = pd.concat(cum, ignore_index=True).groupby(["gcm", "ssp", "co2", "variant", "Year"]).cumNEP.mean()
    gcms, ssps = sorted(p.gcm.unique()), sorted(p.ssp.unique())
    fig, axes = plt.subplots(len(ssps), len(gcms), figsize=(5 * len(gcms), 3.6 * len(ssps)), squeeze=False,
                             sharex=True, sharey=True)
    style = {("same", "SSP"): ("0.45", "-"), ("regen", "SSP"): ("tab:green", "-"),
             ("same", "fixed 2024"): ("0.45", "--"), ("regen", "fixed 2024"): ("tab:green", "--")}
    for i, s_ in enumerate(ssps):
        for j, g_ in enumerate(gcms):
            ax = axes[i, j]
            for (v, co2), (col, ls) in style.items():
                y = q.loc[(g_, s_, co2, v)]
                ax.plot(y.index, y.values, color=col, ls=ls, lw=1.5,
                        label=f"{'regeneration' if v == 'regen' else 'no regeneration'}, CO2 {co2}")
            ax.set_title(f"STOP, {g_}, {s_}", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"cumulative NEP since {START}\n(kg C m$^{{-2}}$)")
            if i == len(ssps) - 1:
                ax.set_xlabel("year")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(RES / "fig_regen.png", dpi=150)
    plt.close(fig)
    print(f"\nWritten: {RES}/pairs.csv, summary.csv, annual_regen.csv, fig_regen.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="STOP with natural regeneration from 2025 (save/restart sensitivity)")
    ap.add_argument("--analyze", action="store_true", help="analyse only, do not build or run")
    ap.add_argument("--workers", type=int, default=4, help="parallel runs (default 4)")
    a = ap.parse_args()
    rows = selected_rows()
    if not a.analyze:
        run_all(rows, a.workers)
    analyse(rows)
