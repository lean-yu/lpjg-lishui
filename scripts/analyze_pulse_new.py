#!/usr/bin/env python3
"""analyze_pulse.py - paired heat-pulse metrics: every pulse run against its control, the main run with the same
factors. Both runs are identical before the pulse year (deterministic model; checked below), so the differences
from the pulse year on are the pulse effect.

Per pulse run (Y = pulse year; differences = pulse minus control, kg C m-2 unless stated):
  age              stand age in Y: years since the control's last clearcut (dens.out TeNE = 0 with harvest)
  vegc_pre         control VegC in Y-1
  dNPP, dNPP_pct   NPP change in Y (absolute, %); dNPP_pct is NaN when the control NPP rounds to 0
  dNEP             ecosystem NEP change in Y, NEP = -(NEE-LU): without harvest and wood-product fluxes
  dNPP_post_pct    NPP change over Y+1..Y+3 (%)
  dVegC_5y         VegC change in Y+5
  dEcoC_10y        ecosystem carbon change (cpool Total-slow) in Y+10
  dMort_3y         change in tree_mort_sts.out (Fir_CC) summed over Y..Y+2
  rec_VegC         years until VegC is back within 1 % of the control, searched up to the control's next clearcut
                   (0 = never more than 1 % below; NaN = not recovered before that clearcut or 2100)
  dRemoval_2100    change in harvest removals (Harvest + product-pool inflow) summed over Y..2100
  dRemoval_10y     change in harvest removals summed over Y..Y+10
  dUptake_10y      change in NEP summed over Y..Y+10 = dEcoC_10y + dRemoval_10y: carbon the stand did not take up,
                   wherever it ended up. Use it to compare managements: where the control is clearcut or thinned
                   in the window, part of the loss leaves as lower removals instead of staying in dEcoC_10y
  harv_events_10y  years in Y..Y+10 in which the control exports harvest (Harvest flux > 0); a thinning that leaves
                   the cut trees on site (the age-9 thinning in the current setup) is not counted
  yrs_to_cc        years from Y to the control's next clearcut (NaN: none by 2100). A clearcut soon after the
                   pulse removes the stand before it can make up the loss, so the loss stays (see the timing table)
Caveats: the model has no physiological legacy (docs/LOG.md), so recovery shows mainly through mortality and
biomass; for very young stands use the absolute metrics (percentages are unstable when NPP is small). Outputs carry
3 decimals, so the placebo's dNPP is 0 or -0.001; its continuous metrics (dNEP, dUptake_10y) give the noise floor.

Output (pulse/): results.csv; phase_means.csv; landscape.csv (landscape BAU = W*bau26 + (1-W)*stop,
EXT = W*ext36 + (1-W)*stop, from phase means); age_vs_mgmt.csv (hot pulses only); background_jas.csv (cache of
the background pulse-season temperatures); fig_pulse_temp.png (pulse effect vs background temperature, colour =
stand age); fig_pulse_age.png (pulse effect vs stand age, adjusted for background climate and heat parameter).
Usage: python analyze_pulse.py            metrics from the run directories, then tables and figures
       python analyze_pulse.py --replot   tables and figures only, from an existing pulse/results.csv
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BASE = Path.home() / "lpjg_lishui"
RUNS = BASE / "runs_main"
OUT = BASE / "pulse"
W = 0.38                                   # commercial share of fir area (docs/sources.md); sensitivity 0.3-0.6
KEYS = ["elev_m", "mgmt", "pstemp_high", "gcm", "ssp", "phase"]
AGE_BINS = [0, 5, 10, 15, 20, 25, 30, 40, 60, 200]
LAND = ["dNPP", "dNEP", "dVegC_5y", "dEcoC_10y", "dRemoval_10y", "dUptake_10y", "dRemoval_2100"]
TEST_METRICS = ("dNPP", "dUptake_10y", "dEcoC_10y")
NEW_COLS = ("dRemoval_10y", "dUptake_10y", "harv_events_10y", "yrs_to_cc")
CELL = (119.75, 28.25)                     # lon, lat of the simulated cell
MGMT_STYLE = {"bau26": ("tab:blue", "o", -0.7), "ext36": ("tab:orange", "s", 0.0), "stop": ("tab:green", "^", 0.7)}
_cache = {}


def table(run, fname):
    if (run, fname) not in _cache:
        p = RUNS / run / "output" / fname
        _cache[run, fname] = (pd.read_csv(p, sep=r"\s+").groupby("Year").mean(numeric_only=True)
                              if p.exists() else None)
    return _cache[run, fname]


def finished(run):
    log = RUNS / run / "run.log"
    return log.exists() and "Finished" in log.read_text(errors="ignore")


def removals(run):
    cf, cp = table(run, "cflux.out"), table(run, "cpool.out")
    return cf["Harvest"] + cp["HarvSlowC"].diff().fillna(0) + cf["Slow_h"]


def clearcuts(run):
    dens, harv = table(run, "dens.out")["TeNE"], table(run, "cflux.out")["Harvest"]
    return [int(y) for y in dens.index if dens[y] < 1e-6 and harv.get(y, 0) > 0]


def recovery(vp, vc, y0, y_end):
    below = [y for y in range(y0, y_end) if y in vc.index and vp[y] < 0.99 * vc[y]]
    if not below:
        return 0
    return max(below) + 1 - y0 if max(below) + 1 < y_end else np.nan


def metrics(prun, crun, Y):
    npp_p, npp_c = table(prun, "anpp.out")["Total"], table(crun, "anpp.out")["Total"]
    nee_p, nee_c = table(prun, "cflux.out")["NEE-LU"], table(crun, "cflux.out")["NEE-LU"]
    veg_p, veg_c = table(prun, "cpool.out")["VegC"], table(crun, "cpool.out")["VegC"]
    eco_p, eco_c = table(prun, "cpool.out")["Total-slow"], table(crun, "cpool.out")["Total-slow"]
    mort_p, mort_c = table(prun, "tree_mort_sts.out")["Fir_CC"], table(crun, "tree_mort_sts.out")["Fir_CC"]
    harv_c = table(crun, "cflux.out")["Harvest"]
    d_rem = removals(prun) - removals(crun)
    cc = clearcuts(crun)
    past, future = [y for y in cc if y <= Y], [y for y in cc if y > Y]
    post = list(range(Y + 1, Y + 4))
    return dict(
        pre_maxdiff=float((eco_p - eco_c).loc[:Y - 1].abs().max()),
        age=Y - max(past) if past else np.nan,
        vegc_pre=veg_c[Y - 1],
        dNPP=npp_p[Y] - npp_c[Y],
        dNPP_pct=(npp_p[Y] / npp_c[Y] - 1) * 100 if npp_c[Y] > 0 else np.nan,
        dNEP=-(nee_p[Y] - nee_c[Y]),
        dNPP_post_pct=(npp_p[post].mean() / npp_c[post].mean() - 1) * 100,
        dVegC_5y=veg_p[Y + 5] - veg_c[Y + 5],
        dEcoC_10y=eco_p[Y + 10] - eco_c[Y + 10],
        dMort_3y=float((mort_p - mort_c).loc[Y:Y + 2].sum()),
        rec_VegC=recovery(veg_p, veg_c, Y, future[0] if future else 2101),
        dRemoval_2100=float(d_rem.loc[Y:].sum()),
        # new columns are appended, so the original column order of results.csv is unchanged
        dRemoval_10y=float(d_rem.loc[Y:Y + 10].sum()),
        dUptake_10y=float((nee_c - nee_p).loc[Y:Y + 10].sum()),
        harv_events_10y=int((harv_c.loc[Y:Y + 10] > 1e-6).sum()),
        yrs_to_cc=future[0] - Y if future else np.nan,
    )


def ols_r2(y, X):
    X = np.column_stack([np.ones(len(y)), X]).astype(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    ssr = float(((y - X @ beta) ** 2).sum())
    return 1 - ssr / float(((y - y.mean()) ** 2).sum()), ssr, np.linalg.matrix_rank(X)


def hot_only(res):
    return res[res.pulse.astype(str).str.startswith("hot")]


def age_vs_mgmt(res):
    """Does management explain the pulse effect beyond stand age? Hot pulses only: the placebo has no effect to
    explain and would inflate R2 through the pulse term. Only age bins with >= 2 managements are used."""
    d = hot_only(res)
    d = d.assign(age_bin=pd.cut(d.age, AGE_BINS, right=False).astype(str)).dropna(subset=["age"])
    d = d[d.groupby("age_bin").mgmt.transform("nunique") >= 2]
    rows = []
    if d.empty or d.mgmt.nunique() < 2:
        return pd.DataFrame(rows), d
    for m in TEST_METRICS:
        e = d.dropna(subset=[m])
        base = pd.get_dummies(e[["age_bin", "elev_m", "pstemp_high", "gcm", "ssp", "pulse"]].astype(str),
                              drop_first=True)
        full = pd.concat([base, pd.get_dummies(e.mgmt, prefix="mgmt", drop_first=True)], axis=1)
        y = e[m].to_numpy(float)
        r2a, ssa, ka = ols_r2(y, base.to_numpy())
        r2b, ssb, kb = ols_r2(y, full.to_numpy())
        dfe = len(y) - kb
        f = ((ssa - ssb) / max(kb - ka, 1)) / (ssb / dfe) if dfe > 0 and ssb > 0 else np.nan
        try:
            from scipy.stats import f as fdist
            p = float(fdist.sf(f, kb - ka, dfe))
        except ImportError:
            p = np.nan
        rows.append(dict(metric=m, pulses="hot only", n=len(y), age_bins=e.age_bin.nunique(), R2_age_and_controls=r2a,
                         R2_plus_mgmt=r2b, dR2_mgmt=r2b - r2a, F=f, p=p))
    return pd.DataFrame(rows), d


def landscape(res):
    ph = res.groupby(["elev_m", "pstemp_high", "gcm", "ssp", "pulse", "mgmt"])[LAND].mean().reset_index()
    wide = ph.pivot_table(index=["elev_m", "pstemp_high", "gcm", "ssp", "pulse"], columns="mgmt", values=LAND)
    out = []
    for name, managed in (("BAU", "bau26"), ("EXT", "ext36"), ("STOP", "stop")):
        part = pd.DataFrame({m: W * wide[m][managed] + (1 - W) * wide[m]["stop"] for m in LAND})
        out.append(part.assign(landscape=name).reset_index())
    return ph, pd.concat(out, ignore_index=True)


def months_of(spec):
    """'7-9' -> [7, 8, 9]"""
    a, b = (int(x) for x in str(spec).split("-"))
    return list(range(a, b + 1))


def tas_path(g, s, e):
    """Control temperature file of a gcm x ssp x band, as named by the runner."""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run_main import climate_paths
        return Path(climate_paths({"gcm": g, "ssp": s, "elev_m": e, "pulse": "none"})["file_temp"])
    except Exception:  # runner not importable: same naming rule as run_main.climate_paths
        return BASE / "data" / "climate3b" / "elev" / f"{g}_{s}_tas_e{e}_1901_2100.nc"


def background_temperature(res, design):
    """Mean temperature of the pulse months in the pulse year, in the control climate of each gcm x ssp x band: the
    background the pulse is added to. Sources in this order: cache pulse/background_jas.csv, export/clim_annual.csv
    (export_for_r.py, Jul-Sep only), the climate files (netCDF4). Returns None if none of them works."""
    spec = {r["pulse"]: r["pulse_months"] for r in design if r["pulse"] != "none"}
    months = sorted({tuple(months_of(v)) for v in spec.values()})
    if len(months) != 1:
        print(f"background temperature: pulses use different months {months}; skipped")
        return None
    months = list(months[0])
    need = res[["gcm", "ssp", "elev_m", "pulse_year"]].drop_duplicates().astype({"elev_m": int, "pulse_year": int})
    key = ["gcm", "ssp", "elev_m", "pulse_year"]
    cache = OUT / "background_jas.csv"
    tag = f"{min(months)}-{max(months)}"
    if cache.exists():
        c = pd.read_csv(cache)
        c = c[c.months == tag]
        got = need.merge(c, on=key, how="left")
        if got.tas_bg.notna().all():
            return got
    clim = BASE / "export" / "clim_annual.csv"
    got = None
    if clim.exists() and months == [7, 8, 9]:
        c = pd.read_csv(clim).rename(columns={"Year": "pulse_year", "tas_jas": "tas_bg"})
        got = need.merge(c[key + ["tas_bg"]], on=key, how="left")
        if got.tas_bg.isna().any():
            got = None
    if got is None:
        try:
            from netCDF4 import Dataset, num2date
        except ImportError:
            print("background temperature: no cache, no export/clim_annual.csv and no netCDF4; fig_pulse_temp skipped")
            return None
        vals = []
        for (g, s, e), grp in need.groupby(["gcm", "ssp", "elev_m"]):
            path = tas_path(g, s, e)
            with Dataset(path) as nc:
                t = nc.variables["time"]
                dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
                lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
                i, j = int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))
                v = nc.variables["tas"]
                x = np.ma.filled(np.ma.asarray(v[:, i, j]).astype(float), np.nan)
                units = str(getattr(v, "units", ""))
            if units.strip().lower() in ("k", "kelvin") or np.nanmean(x) > 100:
                x = x - 273.15
            ym = np.array([(d.year, d.month) for d in dates])
            for y in grp.pulse_year:
                sel = (ym[:, 0] == y) & np.isin(ym[:, 1], months)
                vals.append(dict(gcm=g, ssp=s, elev_m=e, pulse_year=y, tas_bg=float(np.nanmean(x[sel]))))
        got = need.merge(pd.DataFrame(vals), on=key, how="left")
    got.assign(months=tag).round({"tas_bg": 3}).to_csv(cache, index=False)
    return got


def r2_dummies(d, m, cols):
    X = pd.get_dummies(d[cols].astype(str), drop_first=True).to_numpy(float)
    return ols_r2(d[m].to_numpy(float), X)[0]


def adjusted(h, m):
    """Partial residuals: metric m with the effects of background climate (gcm x ssp x band x pulse year) and
    pstemp_high removed. One OLS with these factors and the stand (management x phase x pulse) as dummies; the
    removed part is centred, so the overall mean is kept. What is left is the stand (age, management) signal."""
    d = h.dropna(subset=[m]).copy()
    bg = d.gcm + "|" + d.ssp + "|" + d.elev_m.astype(str) + "|" + d.pulse_year.astype(str)
    stand = d.mgmt + "|" + d.phase.astype(str) + "|" + d.pulse.astype(str)
    Xb = pd.get_dummies(bg, drop_first=True).to_numpy(float)
    Xp = pd.get_dummies(d.pstemp_high.astype(str), drop_first=True).to_numpy(float)
    Xs = pd.get_dummies(stand, drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(d)), Xb, Xp, Xs])
    beta, *_ = np.linalg.lstsq(X, d[m].to_numpy(float), rcond=None)
    nb, nt = Xb.shape[1], Xp.shape[1]
    removed = Xb @ beta[1:1 + nb] + Xp @ beta[1 + nb:1 + nb + nt]
    d["adj"] = d[m].to_numpy(float) - removed + removed.mean()
    return d


def placebo_band(res, m="dUptake_10y"):
    pl = res[res.pulse.astype(str).str.startswith("plc")][m].dropna()
    return float(np.percentile(pl.abs(), 95)) if len(pl) else None


ROWS = (("dNPP_pct", "NPP change in pulse year (%)", "{:+.2f} % per K"),
        ("dUptake_10y", "uptake change, pulse year + 10 yr\n(kg C m-2)", "{:+.4f} kg C m-2 per K"))


def legend(fig, band, extra=()):
    handles = [Line2D([], [], color="0.4", marker=mk, ls="", label=mg) for mg, (_, mk, _) in MGMT_STYLE.items()]
    if band:
        handles.append(Patch(color="0.8", label=f"placebo noise, 95% of |x| ({band:.4f})"))
    fig.legend(handles=handles + list(extra), loc="lower center", ncol=len(handles) + len(extra), fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.02))


def fig_temperature(h, band, path):
    pulses = sorted(h.pulse.unique())
    fig, axes = plt.subplots(2, len(pulses), figsize=(5.2 * len(pulses) + 1.2, 7.5), squeeze=False, sharex=True)
    norm, cmap = Normalize(vmin=0, vmax=max(80, float(np.nanmax(h.age)))), plt.get_cmap("viridis")
    for j, pu in enumerate(pulses):
        sub = h[h.pulse == pu]
        for i, (m, lab, fmt) in enumerate(ROWS):
            ax = axes[i, j]
            if i == 1 and band:
                ax.axhspan(-band, band, color="0.8", zorder=0)
            for mg, (_, mk, _) in MGMT_STYLE.items():
                s = sub[sub.mgmt == mg].dropna(subset=[m])
                ax.scatter(s.tas_bg, s[m], c=s.age, cmap=cmap, norm=norm, marker=mk, s=18, alpha=0.8, linewidths=0)
            s = sub.dropna(subset=[m])
            if len(s) > 2 and s.tas_bg.nunique() > 1:
                b, a = np.polyfit(s.tas_bg, s[m], 1)
                xx = np.array([s.tas_bg.min(), s.tas_bg.max()])
                ax.plot(xx, a + b * xx, color="k", lw=1)
                ax.text(0.98, 0.97, "slope " + fmt.format(b) + f", r = {s.tas_bg.corr(s[m]):.2f}",
                        transform=ax.transAxes, ha="right", va="top", fontsize=8)
            ax.axhline(0, color="grey", lw=0.5)
            if j == 0:
                ax.set_ylabel(lab)
            ax.set_title(pu if i == 0 else "")
        axes[1, j].set_xlabel("Jul-Sep mean temperature of the control climate\nin the pulse year (deg C)")
    legend(fig, band)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes, shrink=0.85, label="stand age in the pulse year (yr)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_age(h, band, path):
    pulses = sorted(h.pulse.unique())
    fig, axes = plt.subplots(2, len(pulses), figsize=(5.2 * len(pulses), 7.5), squeeze=False, sharex=True)
    for j, pu in enumerate(pulses):
        for i, (m, lab, _) in enumerate(ROWS):
            ax = axes[i, j]
            if i == 1 and band:
                ax.axhspan(-band, band, color="0.8", zorder=0)
            d = adjusted(h, m)
            d = d[d.pulse == pu]
            for mg, (col, mk, off) in MGMT_STYLE.items():
                s = d[d.mgmt == mg]
                ax.scatter(s.age + off, s.adj, color=col, marker=mk, s=14, alpha=0.45, linewidths=0)
                mean = s.groupby("age").agg(adj=("adj", "mean"), soon=("yrs_to_cc", lambda x: (x <= 3).all()))
                edge = np.where((i == 1) & mean.soon.to_numpy(), "red", "k")
                ax.scatter(mean.index + off, mean.adj, color=col, marker=mk, s=60, edgecolors=edge,
                           linewidths=np.where(edge == "red", 2.0, 0.8), zorder=3)
            ax.axhline(0, color="grey", lw=0.5)
            if j == 0:
                ax.set_ylabel(lab + "\nadjusted for background climate")
            ax.set_title(pu if i == 0 else "")
        axes[1, j].set_xlabel("stand age in the pulse year (yr)\n(managements offset by 0.7 yr for visibility)")
    soon = Line2D([], [], color="w", marker="o", markeredgecolor="red", markeredgewidth=2, ls="",
                  label="mean of a stand clearcut <= 3 yr after the pulse")
    legend(fig, band, [soon])
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compute(design):
    controls = {tuple(r[k] for k in KEYS): r["run_id"] for r in design if r["pulse"] == "none"}
    rows, missing = [], []
    for r in design:
        if r["pulse"] == "none":
            continue
        crun = controls[tuple(r[k] for k in KEYS)]
        if not (finished(r["run_id"]) and finished(crun)):
            missing.append(r["run_id"])
            continue
        rows.append({**{k: r[k] for k in ["run_id", *KEYS, "pulse", "pulse_year"]}, "control": crun,
                     **metrics(r["run_id"], crun, int(r["pulse_year"]))})
    if not rows:
        raise SystemExit(f"no finished pulse/control pairs ({len(missing)} missing)")
    res = pd.DataFrame(rows)
    for k in ("elev_m", "pstemp_high", "phase", "pulse_year"):
        res[k] = res[k].astype(int)
    res.to_csv(OUT / "results.csv", index=False)
    print(f"{len(res)} pulse runs analysed, {len(missing)} missing (unfinished pulse or control run)")
    return res


def report(res, design):
    pd.set_option("display.width", 220)
    bad = res[res.pre_maxdiff > 1e-6]
    print("pre-pulse check: " + (f"WARNING {len(bad)} pairs differ before the pulse year: {bad.run_id.head(3).tolist()}"
                                 if len(bad) else "all pairs identical before the pulse year"))
    nan_pct = res[res.dNPP_pct.isna()]
    if len(nan_pct):
        print(f"note: {len(nan_pct)} runs have dNPP_pct = NaN (control NPP rounds to 0), e.g. "
              f"{nan_pct.run_id.head(2).tolist()}; they are left out of the % table and the NPP panels")

    for m, label in (("dNPP_pct", "NPP change in the pulse year (%)"),
                     ("dUptake_10y", "uptake change, pulse year + 10 years (kg C m-2) = dEcoC_10y + dRemoval_10y")):
        print(f"\n== {label}; mean over GCMs and phases ==")
        print(res.pivot_table(index=["pstemp_high", "elev_m", "mgmt"], columns=["ssp", "pulse"], values=m)
              .round(3).to_string())

    h = hot_only(res)
    acc = h.groupby("mgmt").agg(dEcoC_10y=("dEcoC_10y", "mean"), dRemoval_10y=("dRemoval_10y", "mean"),
                                dUptake_10y=("dUptake_10y", "mean"),
                                share_with_harvest_in_window=("harv_events_10y", lambda x: (x > 0).mean()))
    print("\n== where the 10-year loss ends up, hot pulses (kg C m-2): ecosystem vs lower removals ==")
    print(acc.round(4).to_string())
    pos = h[h.dEcoC_10y > 0]
    if len(pos):
        print(f"   dEcoC_10y > 0 in {len(pos)} runs, {int((pos.harv_events_10y > 0).sum())} of them with harvest "
              f"in the window")
    t = h.assign(next_clearcut=pd.cut(h.yrs_to_cc.fillna(999), [0, 3, 6, 10, 1000],
                                      labels=["1-3 yr", "4-6 yr", "7-10 yr", ">10 yr or none"]))
    print("\n== uptake loss by years from the pulse to the control's next clearcut, hot pulses ==")
    print(t.groupby("next_clearcut", observed=True).agg(
        n=("run_id", "size"), dNEP_pulse_year=("dNEP", "mean"), dUptake_10y=("dUptake_10y", "mean"),
        share_kept=("dUptake_10y", lambda x: float(np.nanmedian(x / t.loc[x.index, "dNEP"])))).round(3).to_string())
    print("   share_kept: median dUptake_10y / pulse-year dNEP; the rest is made up within the 10 years")

    ph, land = landscape(res)
    land.to_csv(OUT / "landscape.csv", index=False)
    ph.to_csv(OUT / "phase_means.csv", index=False)
    print(f"\n== landscape (phase means, W = {W}): uptake change, pulse year + 10 years (kg C m-2), mean over GCMs ==")
    print(land.pivot_table(index=["pstemp_high", "elev_m", "landscape"], columns=["ssp", "pulse"], values="dUptake_10y")
          .round(3).to_string())

    av, d = age_vs_mgmt(res)
    av.to_csv(OUT / "age_vs_mgmt.csv", index=False)
    print("\n== does management matter beyond stand age? hot pulses only (OLS with age bins, elevation, heat parameter,"
          " GCM, SSP, pulse; age bins shared by >= 2 managements) ==")
    if len(av):
        print(av.round(4).to_string(index=False))
        print("   runs per age bin and management in the test:")
        print("   " + pd.crosstab(d.age_bin, d.mgmt).to_string().replace("\n", "\n   "))
    else:
        print("   not enough overlap between managements")

    band = placebo_band(res)
    pl = res[res.pulse.astype(str).str.startswith("plc")]
    if len(pl):
        print("\n== placebo noise floor (95th percentile of |x|) vs median |x| of the hot pulses ==")
        for m in ("dNEP", "dVegC_5y", "dEcoC_10y", "dUptake_10y"):
            print(f"   {m:12s} placebo {np.nanpercentile(pl[m].abs(), 95):.4f} (max {pl[m].abs().max():.4f})"
                  f"   hot {h[m].abs().median():.4f}")

    if len(h):
        hh = h.dropna(subset=["dNPP_pct"])
        print("\n== share of the variance of the pulse-year NPP change (%), hot pulses ==")
        stand = hh.assign(stand=hh.mgmt + "|" + hh.phase.astype(str) + "|" + hh.pulse)
        stand["bg"] = stand.gcm + "|" + stand.ssp + "|" + stand.elev_m.astype(str) + "|" + stand.pulse_year.astype(str)
        for m in ("dNPP_pct", "dNPP"):
            r_st, r_bg = r2_dummies(stand, m, ["stand"]), r2_dummies(stand, m, ["bg"])
            r_all = r2_dummies(stand, m, ["stand", "bg", "pstemp_high"])
            print(f"   {m:9s} stand (management x phase) alone {r_st:.2f}; background (GCM x SSP x band x year) "
                  f"alone {r_bg:.2f}; both + pstemp_high {r_all:.2f}")

    figs = []
    bt = background_temperature(h, design) if len(h) else None
    if bt is not None:
        hb = h.merge(bt, on=["gcm", "ssp", "elev_m", "pulse_year"], how="left")
        cond = hb.groupby(["gcm", "ssp", "elev_m", "pulse_year"]).agg(dNPP_pct=("dNPP_pct", "mean"),
                                                                      tas_bg=("tas_bg", "first"))
        if len(cond) > 2:
            b, _ = np.polyfit(cond.tas_bg, cond.dNPP_pct, 1)
            print(f"\n== pulse-year NPP change vs background temperature of the pulse months ({len(cond)} GCM x SSP x "
                  f"band x year means): {b:+.2f} % per K, r = {cond.tas_bg.corr(cond.dNPP_pct):.2f} ==")
        fig_temperature(hb, band, OUT / "fig_pulse_temp.png")
        figs.append("fig_pulse_temp.png")
    if len(h):
        fig_age(h, band, OUT / "fig_pulse_age.png")
        figs.append("fig_pulse_age.png")
    print(f"\nWritten: {OUT}/results.csv, landscape.csv, phase_means.csv, age_vs_mgmt.csv"
          + (", background_jas.csv" if bt is not None else "") + (", " + ", ".join(figs) if figs else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="paired heat-pulse metrics, tables and figures")
    ap.add_argument("--replot", action="store_true", help="tables and figures only, from an existing pulse/results.csv")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    design = list(csv.DictReader(open(BASE / "design.csv")))
    if a.replot:
        res = pd.read_csv(OUT / "results.csv")
        lacking = [c for c in NEW_COLS if c not in res.columns]
        if lacking:
            raise SystemExit(f"pulse/results.csv lacks {lacking}: run once without --replot")
    else:
        res = compute(design)
    report(res, design)
