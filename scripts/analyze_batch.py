#!/usr/bin/env python3
# analyze_batch.py - metrics from the density x management x root-depth x climate batch
# Output: batch/results.csv (one row per disturbed run), tables on screen, batch/fig_resistance.png
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

X = Path.home() / "lpjg_lishui" / "batch"
PLANT_YEAR = 1988
man = pd.read_csv(X / "manifest.csv")
_cache = {}


def series(run, fname, col="Total"):
    """Yearly series from one output file; col=None sums all value columns."""
    key = (run, fname, col)
    if key not in _cache:
        try:
            df = pd.read_csv(X / "runs" / run / "output" / fname, sep=r"\s+")
            vals = df.drop(columns=[c for c in ("Lon", "Lat", "Year") if c in df.columns])
            s = vals.select_dtypes("number").sum(axis=1) if col is None else df[col]
            _cache[key] = s.groupby(df["Year"]).mean()
        except Exception:  # missing file or unexpected format
            _cache[key] = pd.Series(dtype=float)
    return _cache[key]


def mean_years(s, years):
    v = s.reindex(list(years))
    return v.mean() if v.notna().all() else np.nan


# ---- 1. control runs: growth check ----
ctrl = man[man.scenario == "ctrl"].copy()
ctrl["NPP_1995_99"] = [mean_years(series(r, "anpp.out"), range(1995, 2000)) for r in ctrl.run]
ctrl["VegC_2013"] = [series(r, "cpool.out", "VegC").get(2013, np.nan) for r in ctrl.run]
ctrl["stems_2013"] = [series(r, "dens.out", "TeNE").get(2013, np.nan) * 10000 for r in ctrl.run]
print("== Control runs (no extra disturbance): growth by density, management and root depth ==")
print(ctrl.pivot_table(index=["mgmt", "density"], columns="beta",
                       values=["NPP_1995_99", "VegC_2013", "stems_2013"]).round(2).to_string(), "\n")

# ---- 2. disturbed runs: paired metrics and Lloret indices ----
rows = []
for r in man.itertuples():
    if r.scenario == "ctrl":
        continue
    kind, y = r.scenario.rsplit("_", 1)
    y = int(y)
    ev = [y, y + 1] if kind == "dry50x2" else [y]
    pre, post = range(ev[0] - 3, ev[0]), range(ev[-1] + 1, ev[-1] + 4)
    c = r.run[: -len(r.scenario)] + "ctrl"
    nd, nc = series(r.run, "anpp.out"), series(c, "anpp.out")
    vd, vc = series(r.run, "cpool.out", "VegC"), series(c, "cpool.out", "VegC")
    md, mc = series(r.run, "mort_greff_sts.out", None), series(c, "mort_greff_sts.out", None)
    during_d, during_c = mean_years(nd, ev), mean_years(nc, ev)
    pre_d, post_d, post_c = mean_years(nd, pre), mean_years(nd, post), mean_years(nc, post)
    win = list(range(ev[0], ev[-1] + 6))                 # event year(s) + 5 years
    mwin = list(range(ev[0], ev[-1] + 3))
    rows.append(dict(
        run=r.run, density=r.density, mgmt=r.mgmt, beta=r.beta, kind=kind, year=y, age=y - PLANT_YEAR,
        dNPP_event_pct=(during_d / during_c - 1) * 100,      # paired resistance: event NPP vs control
        dNPP_post_pct=(post_d / post_c - 1) * 100,           # 3 post-event years vs control
        cum_dNPP=(nd.reindex(win) - nc.reindex(win)).sum(min_count=1),       # kg C m-2
        dVegC_5y=vd.get(ev[-1] + 5, np.nan) - vc.get(ev[-1] + 5, np.nan),   # kg C m-2
        dMort_greff=md.reindex(mwin).sum(min_count=1) - mc.reindex(mwin).sum(min_count=1),
        Rt=during_d / pre_d, Rc=post_d / during_d, Rs=post_d / pre_d, RRs=(post_d - during_d) / pre_d,
    ))
res = pd.DataFrame(rows)
res.to_csv(X / "results.csv", index=False)
print(f"Written: {X / 'results.csv'} ({len(res)} rows)\n")
pd.set_option("display.width", 200)

names = {"dry50": "Jul-Sep rain -50%", "dry50x2": "rain -50% two years", "hot3": "Jul-Sep +3 C",
         "comp": "rain -50% and +3 C"}
for m, label in (("A", "no thinning"), ("B", "standard thinning")):
    print(f"######## {label}: event NPP change vs control (%), root_beta 0.96 (range over 0.95-0.97) ########")
    for k in ("dry50", "dry50x2", "hot3", "comp"):
        sub = res[(res.mgmt == m) & (res.kind == k)]
        if sub.empty:
            continue
        center = sub[np.isclose(sub.beta, 0.96)].pivot_table(index="density", columns="age", values="dNPP_event_pct")
        lo = sub.pivot_table(index="density", columns="age", values="dNPP_event_pct", aggfunc="min")
        hi = sub.pivot_table(index="density", columns="age", values="dNPP_event_pct", aggfunc="max")
        table = center.round(1).astype(str) + " (" + lo.round(1).astype(str) + "~" + hi.round(1).astype(str) + ")"
        table.columns = [f"age{a}" for a in table.columns]
        print(f"-- {names[k]} --")
        print(table.to_string(), "\n")

print("== Mean over all scenarios and ages, root_beta 0.96 ==")
summ = res[np.isclose(res.beta, 0.96)].groupby(["mgmt", "density"])[
    ["dNPP_event_pct", "dNPP_post_pct", "cum_dNPP", "dVegC_5y", "dMort_greff", "Rt", "Rs"]].mean().round(3)
print(summ.to_string(), "\n")

# ---- 3. figure ----
kinds = [k for k in ("dry50", "dry50x2", "hot3", "comp") if k in set(res.kind)]
fig, axes = plt.subplots(1, len(kinds), figsize=(4 * len(kinds), 3.6), sharey=True)
axes = np.atleast_1d(axes)
for ax, k in zip(axes, kinds):
    for m, ls in (("A", "-"), ("B", "--")):
        for age, colr in zip(sorted(res.age.unique()), ("tab:blue", "tab:orange", "tab:red")):
            sub = res[(res.kind == k) & (res.mgmt == m) & (res.age == age)]
            g = sub.groupby("density")["dNPP_event_pct"]
            cen = sub[np.isclose(sub.beta, 0.96)].groupby("density")["dNPP_event_pct"].mean()
            ax.plot(cen.index, cen.values, ls, color=colr, marker="o", label=f"{m} age {age}")
            ax.fill_between(g.min().index, g.min().values, g.max().values, color=colr, alpha=0.08)
    ax.set_title(names[k])
    ax.set_xlabel("planting density (stems/ha)")
    ax.axhline(0, color="grey", lw=0.5)
axes[0].set_ylabel("event NPP change vs control (%)")
axes[-1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(X / "fig_resistance.png", dpi=150)
print(f"Figure: {X / 'fig_resistance.png'}")
