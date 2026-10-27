#!/usr/bin/env python3
"""co2fix.py - CO2 attribution: the main runs again, with CO2 held at its 2024 value from 2025 on.

Each run is paired with the main run of the same run_id (runs_main/), which has the SSP's rising CO2. Both are
identical up to 2024 (checked below), so from 2025 on:
  CO2 effect      = main run - fixed-CO2 run (same SSP climate)
  climate effect  = fixed-CO2 SSP5-8.5 - fixed-CO2 SSP1-2.6 (their CO2 differs only by the 2024 values, printed)
  SSP difference  = climate effect + (CO2 effect in SSP5-8.5 - CO2 effect in SSP1-2.6)
Selection: pstemp_high 25 (the heat parameter changes cumulative NEP by ~5 % only), both GCMs and SSPs, all bands,
managements and phases = 144 runs, about 30 min with 4 workers. The runs are built by run_main.build (same template,
calibrated fir, management switch); only file_co2 then points to data/co2/co2_<ssp>_fix2024_1850_2100_lpjg.txt.
Output: runs_co2fix/<run_id>/; co2fix/pairs.csv (one row per pair), summary.csv, annual_co2fix.csv (fixed-CO2
runs, 2020-2100, same columns as export/annual_main.csv where available), fig_co2fix.png.
Usage:  python co2fix.py               build and run what is not finished yet, then analyse
        python co2fix.py --analyze     analyse only
        python co2fix.py --workers N   parallel runs (default 4)
"""
import argparse
import csv
import datetime
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_main as rm  # noqa: E402  same build rules as the main runs

BASE = rm.BASE
MAIN = BASE / "runs_main"
FIX = BASE / "runs_co2fix"
RES = BASE / "co2fix"
FIX_YEAR = 2024
PSTEMP = "25"
EARLY, LATE = (2025, 2044), (2081, 2100)


def co2_paths(ssp):
    src = BASE / "data" / "co2" / f"co2_{ssp}_1850_2100_lpjg.txt"
    return src, src.with_name(f"co2_{ssp}_fix{FIX_YEAR}_1850_2100_lpjg.txt")


def make_co2(ssp):
    """Copy of the SSP CO2 file with every year after FIX_YEAR set to the FIX_YEAR value. Returns (path, ppm 2024,
    ppm 2100 of the rising path)."""
    src, dst = co2_paths(ssp)
    data = [ln.split() for ln in src.read_text().splitlines() if ln.split() and ln.split()[0].isdigit()]
    vals = {int(y): float(v) for y, v, *_ in data}
    if FIX_YEAR not in vals or max(vals) < 2100:
        raise SystemExit(f"{src}: no value for {FIX_YEAR} or does not reach 2100")
    lines = [f"{int(y)} {vals[FIX_YEAR] if int(y) > FIX_YEAR else float(v):.3f}" for y, v, *_ in data]
    tmp = dst.with_name(dst.name + ".part")
    tmp.write_text("\n".join(lines) + "\n")
    tmp.rename(dst)
    return dst, vals[FIX_YEAR], vals[2100]


def selected_rows():
    rows = list(csv.DictReader(open(BASE / "design.csv")))
    return [r for r in rows if r["pulse"] == "none" and r["pstemp_high"] == PSTEMP]


def build_and_run(rows, workers):
    co2 = {}
    for s in sorted({r["ssp"] for r in rows}):
        co2[s], c24, c2100 = make_co2(s)
        print(f"CO2 {s}: held at {c24:.1f} ppm from {FIX_YEAR + 1} (rising path reaches {c2100:.1f} ppm in 2100)")
    rm.OUT = FIX                                   # run_main.build puts each run under rm.OUT / run_id
    commit = rm.git_commit()
    todo = []
    for r in rows:
        d = rm.build(r, commit)
        if d is None:                              # finished in an earlier call
            continue
        cf = d / "lishui_cf.ins"
        cf.write_text(rm.set_param(cf.read_text(), "file_co2", str(co2[r["ssp"]])))
        with open(d / "run_info.txt", "a") as fh:
            fh.write(f"co2=held at the {FIX_YEAR} value from {FIX_YEAR + 1} ({co2[r['ssp']].name})\n")
        todo.append(d)
    print(f"{len(rows)} fixed-CO2 runs, {len(todo)} to run (git {commit}), start {datetime.datetime.now():%H:%M:%S}",
          flush=True)
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(rm.run, d) for d in todo]
        for i, f in enumerate(as_completed(futs), 1):
            name, ok = f.result()
            print(f"{datetime.datetime.now():%H:%M:%S} {i}/{len(todo)} {name} {'ok' if ok else 'FAILED'}", flush=True)
            if not ok:
                failed.append(name)
    print(f"runs finished {datetime.datetime.now():%H:%M:%S}; failed: {len(failed)} {failed[:5]}")


def analyse(rows):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def finished(d):
        log = d / "run.log"
        return log.exists() and "Finished" in log.read_text(errors="ignore")

    def read(d, f):
        p = d / "output" / f
        return pd.read_csv(p, sep=r"\s+").groupby("Year").mean(numeric_only=True)

    def series(d):
        n, g, cf, cp = (read(d, f) for f in ("anpp.out", "agpp.out", "cflux.out", "cpool.out"))
        x = pd.DataFrame({"GPP": g["Total"], "NPP": n["Total"], "NEP": -cf["NEE-LU"], "VegC": cp["VegC"],
                          "HarvSlowC": cp["HarvSlowC"], "EcoC": cp["Total-slow"]})
        x["removals"] = cf["Harvest"] + cp["HarvSlowC"].diff().fillna(0) + cf["Slow_h"]
        return x

    def totals(x):
        return dict(NEP=x.NEP.loc[2025:2100].sum(), removals=x.removals.loc[2025:2100].sum(),
                    dEcoC=x.EcoC[2100] - x.EcoC[2024], dVegC=x.VegC[2100] - x.VegC[2024],
                    dProducts=x.HarvSlowC[2100] - x.HarvSlowC[2024],
                    NPP_early=x.NPP.loc[EARLY[0]:EARLY[1]].mean(), NPP_late=x.NPP.loc[LATE[0]:LATE[1]].mean())

    RES.mkdir(exist_ok=True)
    keys = ["run_id", "elev_m", "mgmt", "gcm", "ssp", "phase"]
    pairs, ann, cum, missing = [], [], [], []
    for r in rows:
        dm, df = MAIN / r["run_id"], FIX / r["run_id"]
        if not (finished(dm) and finished(df)):
            missing.append(r["run_id"])
            continue
        xm, xf = series(dm), series(df)
        fac = {k: r[k] for k in keys}
        row = dict(fac, pre_maxdiff=float((xf.EcoC - xm.EcoC).loc[:FIX_YEAR].abs().max()))
        row.update({f"{k}_main": v for k, v in totals(xm).items()})
        row.update({f"{k}_fix": v for k, v in totals(xf).items()})
        pairs.append(row)
        ann.append(xf.loc[2020:2100].reset_index().assign(**fac))
        cum.append(pd.DataFrame({"Year": range(2025, 2101), "main": xm.NEP.loc[2025:2100].cumsum().to_numpy(),
                                 "fix": xf.NEP.loc[2025:2100].cumsum().to_numpy()}).assign(**fac))
    if not pairs:
        raise SystemExit(f"no finished pairs ({len(missing)} missing)")
    p = pd.DataFrame(pairs)
    p["elev_m"], p["phase"] = p.elev_m.astype(int), p.phase.astype(int)
    p["CO2_effect_NEP"] = p.NEP_main - p.NEP_fix
    p.to_csv(RES / "pairs.csv", index=False)
    a = pd.concat(ann, ignore_index=True)
    a[keys + [c for c in a.columns if c not in keys]].round(5).to_csv(RES / "annual_co2fix.csv", index=False)
    pd.set_option("display.width", 200)
    print(f"\n{len(p)} pairs analysed, {len(missing)} missing (unfinished main or fixed-CO2 run)")
    bad = p[p.pre_maxdiff > 1e-6]
    print("pre-2025 check: " + (f"WARNING {len(bad)} pairs differ before {FIX_YEAR + 1}: {bad.run_id.head(3).tolist()}"
                                if len(bad) else f"all pairs identical up to {FIX_YEAR}"))

    g = p.groupby(["gcm", "ssp", "mgmt"])
    t = pd.DataFrame({"NEP_SSP_CO2": g.NEP_main.mean(), "NEP_fixed_CO2": g.NEP_fix.mean(),
                      "CO2_effect": g.CO2_effect_NEP.mean(),
                      "CO2_share_%": 100 * g.CO2_effect_NEP.mean() / g.NEP_main.mean()})
    print("\n== cumulative NEP 2025-2100 (kg C m-2), mean over bands and phases ==")
    print(t.round(2).to_string())

    w = p.groupby(["gcm", "mgmt", "ssp"])[["NEP_main", "NEP_fix"]].mean().unstack("ssp")
    if {"ssp126", "ssp585"} <= set(p.ssp):
        dec = pd.DataFrame({"total": w.NEP_main.ssp585 - w.NEP_main.ssp126,
                            "climate": w.NEP_fix.ssp585 - w.NEP_fix.ssp126})
        dec["CO2"] = dec.total - dec.climate
        print("\n== SSP5-8.5 minus SSP1-2.6, cumulative NEP 2025-2100 (kg C m-2) = climate + CO2 ==")
        print("   climate = difference of the fixed-CO2 runs; CO2 = the rest (CO2 effect in SSP5-8.5 minus in SSP1-2.6)")
        print(dec.round(2).to_string())
        dec.reset_index().round(4).to_csv(RES / "summary.csv", index=False)
    else:
        t.reset_index().round(4).to_csv(RES / "summary.csv", index=False)

    b = p[p.mgmt == "bau26"].groupby(["gcm", "ssp"])
    npp = pd.DataFrame({"SSP_CO2": 100 * (b.NPP_late_main.mean() / b.NPP_early_main.mean() - 1),
                        "fixed_CO2": 100 * (b.NPP_late_fix.mean() / b.NPP_early_fix.mean() - 1)})
    print(f"\n== bau26 NPP change, {LATE[0]}-{LATE[1]} vs {EARLY[0]}-{EARLY[1]} (%; rotation-stationary) ==")
    print(npp.round(1).to_string())

    e = p.groupby(["ssp", "elev_m", "mgmt"]).CO2_effect_NEP.mean().unstack("mgmt")
    print("\n== CO2 effect on cumulative NEP by band (kg C m-2), mean over GCMs and phases ==")
    print(e.round(2).to_string())

    s = p.groupby("mgmt")[["dEcoC_main", "dEcoC_fix", "removals_main", "removals_fix", "dProducts_main",
                           "dProducts_fix"]].mean()
    print("\n== where the CO2 effect goes, 2025-2100 (kg C m-2), mean over all ==")
    print(s.round(2).to_string())

    c = pd.concat(cum, ignore_index=True).groupby(["gcm", "ssp", "mgmt", "Year"])[["main", "fix"]].mean()
    gcms, ssps = sorted(p.gcm.unique()), sorted(p.ssp.unique())
    fig, axes = plt.subplots(len(ssps), len(gcms), figsize=(5 * len(gcms), 3.6 * len(ssps)), squeeze=False,
                             sharex=True, sharey=True)
    colors = {"bau26": "tab:blue", "ext36": "tab:orange", "stop": "tab:green"}
    for i, s_ in enumerate(ssps):
        for j, g_ in enumerate(gcms):
            ax = axes[i, j]
            for m, col in colors.items():
                if (g_, s_, m) not in c.index.droplevel("Year"):
                    continue
                q = c.loc[(g_, s_, m)]
                ax.plot(q.index, q["main"], color=col, lw=1.5, label=f"{m}, SSP CO2")
                ax.plot(q.index, q["fix"], color=col, lw=1.5, ls="--", label=f"{m}, CO2 fixed at {FIX_YEAR}")
            ax.set_title(f"{g_}, {s_}", fontsize=10)
            if j == 0:
                ax.set_ylabel("cumulative NEP since 2025\n(kg C m$^{-2}$)")
            if i == len(ssps) - 1:
                ax.set_xlabel("year")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(RES / "fig_co2fix.png", dpi=150)
    plt.close(fig)
    print(f"\nWritten: {RES}/pairs.csv, summary.csv, annual_co2fix.csv, fig_co2fix.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="CO2 attribution runs (CO2 held at its 2024 value) and their analysis")
    ap.add_argument("--analyze", action="store_true", help="analyse only, do not build or run")
    ap.add_argument("--workers", type=int, default=4, help="parallel runs (default 4)")
    a = ap.parse_args()
    rows = selected_rows()
    if not a.analyze:
        build_and_run(rows, a.workers)
    analyse(rows)
