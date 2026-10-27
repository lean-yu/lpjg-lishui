#!/usr/bin/env python3
"""export_for_r.py - compact tables for the paper analysis in R (analysis/paper_analysis.R).

Writes export/ (plain CSV files) and export.tar.gz:
  design.csv            design rows (copy)
  annual_main.csv       main runs, 1990-2100, per year: GPP, NPP, Ra, Rh, NEP (= -(NEE-LU)), removals
                        (= Harvest + product-pool inflow), pools (VegC, TreeC, LitterC, SoilC, HarvSlowC,
                        EcoC = Total-slow), stems_ha, LAI, mort, stand age, clearcut flag
  clearcuts.csv         every clearcut of the main runs since 1950: carbon payback time (years until ecosystem
                        carbon is back at its pre-harvest level; status reached / cut_before / end_of_run) and
                        sink recovery time (years until the 3-year mean NEP is positive)
  pulse_metrics.csv     analyze_pulse metrics per pulse run plus mechanism and uptake columns:
                        gpp/npp of the control in the pulse year, dGPP_pct, dRa, dRh, JAS temperature of the
                        pulse season, dUptake_10y = sum of the NEP change over Y..Y+10 (= dEcoC_10y + dRemoval_10y)
  pulse_diff.csv        pulse minus control for years Y-1..Y+20
  clim_annual.csv       simulated cell, per GCM x SSP x band and year: annual and Jul-Sep tas (deg C), days with
                        daily mean tas > 30 deg C, annual and Jul-Sep precipitation (mm)
  obs_clim_annual.csv, obs_calib_annual.csv   observed climate and the observed-climate calibration runs
                        (respcoeff 0.40, pstemp_high 25 and 30), if present
  batch_results.csv, calib_*.csv               copies of earlier results
  checks.csv            reference numbers that paper_analysis.R recomputes and compares
Usage (server):  source ~/lpjg_venv/bin/activate && python scripts/export_for_r.py
"""
import csv
import re
import shutil
import sys
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_pulse import BASE, KEYS, RUNS, finished, metrics  # noqa: E402
from run_main import climate_paths  # noqa: E402

EXP = BASE / "export"
Y0, Y1 = 1990, 2100
W = 0.38                      # commercial share of fir area (docs/sources.md), for the reference check only
CELL = (119.75, 28.25)
_cache = {}


def read_out(d, fname):
    p = d / "output" / fname
    if not p.exists():
        return None
    return pd.read_csv(p, sep=r"\s+").groupby("Year").mean(numeric_only=True)


def annual(d):
    """Yearly series of one run directory (all simulated years)."""
    if d in _cache:
        return _cache[d]
    g, n, cf, cp, cm = (read_out(d, f) for f in ("agpp.out", "anpp.out", "cflux.out", "cpool.out", "cmass.out"))
    dens, lai, mort = (read_out(d, f) for f in ("dens.out", "lai.out", "tree_mort_sts.out"))
    x = pd.DataFrame({
        "GPP": g["Total"], "NPP": n["Total"], "Rh": cf["Soil"], "NEP": -cf["NEE-LU"],
        "harvest_imm": cf["Harvest"], "slow_h": cf["Slow_h"],
        "VegC": cp["VegC"], "TreeC": cm["TeNE"], "LitterC": cp["LitterC"], "SoilC": cp["SoilC"],
        "HarvSlowC": cp["HarvSlowC"], "EcoC": cp["Total-slow"], "stems_ha": dens["TeNE"] * 1e4,
        "LAI": lai["TeNE"] if lai is not None else np.nan,
        "mort": mort["Fir_CC"] if mort is not None else np.nan})
    x["Ra"] = x.GPP - x.NPP
    x["removals"] = x.harvest_imm + x.HarvSlowC.diff().fillna(0) + x.slow_h
    cc = (x.stems_ha < 0.01) & (x.harvest_imm > 0)
    x["clearcut"] = cc.astype(int)
    last, age = np.nan, []
    for y, c in cc.items():
        last = y if c else last
        age.append(y - last)
    x["age"] = age
    x.index.name = "Year"
    _cache[d] = x
    return x


def clearcut_events(x):
    ccy = [int(y) for y in x.index[x.clearcut == 1]]
    nep3 = x.NEP.rolling(3).mean().shift(-2)            # mean NEP of years y..y+2
    out = []
    for i, c in enumerate(ccy):
        if c < 1950:
            continue
        nxt = ccy[i + 1] if i + 1 < len(ccy) else None
        end = nxt if nxt is not None else int(x.index.max()) + 1
        pre = x.EcoC.get(c - 1, np.nan)
        pay = next((y - c for y in range(c + 1, end) if x.EcoC[y] >= pre), None)
        rec = next((y - c for y in range(c + 1, end) if nep3.get(y, np.nan) > 0), None)
        out.append(dict(cc_year=c, age_at_cut=c - ccy[i - 1] if i > 0 else np.nan, eco_pre=pre,
                        veg_pre=x.VegC.get(c - 1, np.nan), removal=x.removals[c], next_cc=nxt,
                        payback_yr=pay,
                        payback_status="reached" if pay is not None else ("cut_before" if nxt else "end_of_run"),
                        sink_rec_yr=rec))
    return out


def cell_daily(path, var):
    from netCDF4 import Dataset, num2date
    with Dataset(path) as nc:
        t = nc.variables["time"]
        dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
        lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
        i, j = int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))
        v = nc.variables[var]
        vals = np.ma.filled(np.ma.asarray(v[:, i, j]).astype(float), np.nan)
        units = str(getattr(v, "units", ""))
    ym = np.array([(d.year, d.month) for d in dates])
    if var == "tas" and (units.strip().lower() in ("k", "kelvin") or np.nanmean(vals) > 100):
        vals = vals - 273.15
    if var == "pr" and ("s-1" in units or np.nanmax(vals) < 0.1):
        vals = vals * 86400.0                          # kg m-2 s-1 -> mm per day
    return pd.DataFrame({"Year": ym[:, 0], "month": ym[:, 1], var: vals})


def yearly_climate(tas_path, pr_path):
    t, p = cell_daily(tas_path, "tas"), cell_daily(pr_path, "pr")
    jas_t, jas_p = t[t.month.between(7, 9)], p[p.month.between(7, 9)]
    return pd.DataFrame({
        "tas_ann": t.groupby("Year").tas.mean(), "tas_jas": jas_t.groupby("Year").tas.mean(),
        "hot30_days": t.assign(h=t.tas > 30).groupby("Year").h.sum(),
        "pr_ann": p.groupby("Year").pr.sum(), "pr_jas": jas_p.groupby("Year").pr.sum()}).reset_index()


def main():
    EXP.mkdir(exist_ok=True)
    design = list(csv.DictReader(open(BASE / "design.csv")))
    shutil.copy(BASE / "design.csv", EXP / "design.csv")
    main_rows = [r for r in design if r["pulse"] == "none"]

    # --- main runs: annual series and clearcut events ---
    ann, ccs = [], []
    for r in main_rows:
        if not finished(r["run_id"]):
            print(f"skip unfinished {r['run_id']}")
            continue
        x = annual(RUNS / r["run_id"])
        fac = {k: r[k] for k in ["run_id", *KEYS]}
        ann.append(x.loc[Y0:Y1].reset_index().assign(**fac))
        ccs += [dict(fac, **e) for e in clearcut_events(x)]
    ann = pd.concat(ann, ignore_index=True)
    ann = ann[["run_id", *KEYS] + [c for c in ann.columns if c not in ("run_id", *KEYS)]]
    ann.round(5).to_csv(EXP / "annual_main.csv", index=False)
    ccs = pd.DataFrame(ccs)
    ccs.to_csv(EXP / "clearcuts.csv", index=False)

    # --- climate of the simulated cell ---
    clim = []
    for g, s, e in sorted({(r["gcm"], r["ssp"], r["elev_m"]) for r in main_rows}):
        cp = climate_paths({"gcm": g, "ssp": s, "elev_m": e})
        try:
            clim.append(yearly_climate(cp["file_temp"], cp["file_prec"]).assign(gcm=g, ssp=s, elev_m=int(e)))
        except Exception as err:  # keep exporting even without netCDF4 or files
            warnings.warn(f"climate {g} {s} {e}: {err}")
    clim = pd.concat(clim, ignore_index=True) if clim else pd.DataFrame(
        columns=["Year", "tas_ann", "tas_jas", "hot30_days", "pr_ann", "pr_jas", "gcm", "ssp", "elev_m"])
    clim.round(3).to_csv(EXP / "clim_annual.csv", index=False)
    jas = {(row.gcm, row.ssp, str(row.elev_m), int(row.Year)): row.tas_jas for row in clim.itertuples()}

    # --- pulse runs: analyze_pulse metrics + mechanism and uptake columns, and difference time series ---
    controls = {tuple(r[k] for k in KEYS): r["run_id"] for r in main_rows}
    pm, pdiff = [], []
    for r in design:
        if r["pulse"] == "none":
            continue
        crun = controls[tuple(r[k] for k in KEYS)]
        if not (finished(r["run_id"]) and finished(crun)):
            continue
        Y = int(r["pulse_year"])
        xp, xc = annual(RUNS / r["run_id"]), annual(RUNS / crun)
        d = xp - xc
        m = metrics(r["run_id"], crun, Y)
        m.update(
            gpp_ctrl=xc.GPP[Y], npp_ctrl=xc.NPP[Y], npp_gpp_ctrl=xc.NPP[Y] / xc.GPP[Y] if xc.GPP[Y] > 0 else np.nan,
            dGPP_pct=(xp.GPP[Y] / xc.GPP[Y] - 1) * 100 if xc.GPP[Y] > 0 else np.nan, dRa=d.Ra[Y], dRh=d.Rh[Y],
            dUptake_10y=d.NEP.loc[Y:Y + 10].sum(), dRemoval_10y=d.removals.loc[Y:Y + 10].sum(),
            tas_jas=jas.get((r["gcm"], r["ssp"], r["elev_m"], Y), np.nan) + float(r["pulse_dT"]))
        pm.append({**{k: r[k] for k in ["run_id", *KEYS, "pulse", "pulse_year", "pulse_dT"]}, "control": crun, **m})
        for rel in range(-1, 21):
            y = Y + rel
            if y in d.index:
                pdiff.append({**{k: r[k] for k in ["run_id", *KEYS, "pulse"]}, "rel_year": rel, "Year": y,
                              "npp_ctrl": xc.NPP[y], **{f"d{v}": d[v][y] for v in
                                                         ("NPP", "GPP", "Ra", "Rh", "NEP", "VegC", "EcoC", "removals")}})
    pm, pdiff = pd.DataFrame(pm), pd.DataFrame(pdiff)
    pm.round(5).to_csv(EXP / "pulse_metrics.csv", index=False)
    pdiff.round(5).to_csv(EXP / "pulse_diff.csv", index=False)
    bal = (pm.dEcoC_10y - (pm.dUptake_10y - pm.dRemoval_10y)).abs().max() if len(pm) else np.nan

    # --- observed climate and observed-climate calibration runs ---
    cfins = (BASE / "exp_base" / "lishui_cf.ins").read_text()
    obs = {k: re.search(rf'param "{k}"\s*\(str "([^"]+)"\)', cfins) for k in ("file_temp", "file_prec")}
    try:
        yearly_climate(Path(obs["file_temp"].group(1)), Path(obs["file_prec"].group(1))).round(3).to_csv(
            EXP / "obs_clim_annual.csv", index=False)
    except Exception as err:
        warnings.warn(f"observed climate: {err}")
    oc = []
    for d in sorted((BASE / "calib").glob("respcoeff_pstemp_high*/respcoeff_0.4*")):
        if (d / "run.log").exists() and "Finished" in (d / "run.log").read_text(errors="ignore"):
            t = int(re.search(r"pstemp_high(\d+)", str(d)).group(1))
            oc.append(annual(d).reset_index().assign(pstemp_high=t, run=str(d.relative_to(BASE))))
    if oc:
        pd.concat(oc, ignore_index=True).round(5).to_csv(EXP / "obs_calib_annual.csv", index=False)

    # --- copies of earlier results ---
    copies = {"batch_results.csv": BASE / "batch" / "results.csv", "calib_targets.csv": BASE / "docs" / "calib_targets.csv"}
    copies.update({p.name: p for p in (BASE / "docs").glob("calib_*summary*.csv")})
    for name, src in copies.items():
        if src.exists():
            shutil.copy(src, EXP / name)

    # --- reference numbers for paper_analysis.R ---
    a = ann[ann.Year.isin([2024, 2100])].pivot_table(index=["run_id", *KEYS], columns="Year", values="EcoC").reset_index()
    a["dEco"] = a[2100] - a[2024]
    ph = a.groupby(["elev_m", "pstemp_high", "gcm", "ssp", "mgmt"]).dEco.mean().unstack("mgmt")
    hot = pm[pm.pulse.str.startswith("hot")] if len(pm) else pm
    post = ccs[(ccs.mgmt == "bau26") & (ccs.cc_year >= 2025) & (ccs.payback_status == "reached")] if len(ccs) else ccs
    checks = {"n_main_runs": ann.run_id.nunique(), "n_pulse_pairs": len(pm),
              "dEco_2025_2100_landscape_BAU_mean": float((W * ph["bau26"] + (1 - W) * ph["stop"]).mean()),
              "payback_bau26_post2025_mean": float(post.payback_yr.mean()) if len(post) else np.nan,
              "hot_dUptake_10y_mean": float(hot.dUptake_10y.mean()) if len(hot) else np.nan,
              "hot_dNPP_pct_mean": float(hot.dNPP_pct.mean()) if len(hot) else np.nan}
    pd.Series(checks, name="value").rename_axis("name").reset_index().to_csv(EXP / "checks.csv", index=False)

    with tarfile.open(BASE / "export.tar.gz", "w:gz") as tar:
        tar.add(EXP, arcname="export")
    print(f"main runs {checks['n_main_runs']}, clearcut events {len(ccs)}, pulse pairs {len(pm)}, "
          f"climate rows {len(clim)}; carbon balance check (dEcoC_10y vs uptake - removals) max |diff| {bal:.4f}")
    for k, v in checks.items():
        print(f"  check {k} = {v}")
    for f in sorted(EXP.iterdir()):
        print(f"  {f.name:40s} {f.stat().st_size / 1e6:6.2f} MB")
    print(f"archive: {BASE / 'export.tar.gz'} ({(BASE / 'export.tar.gz').stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
