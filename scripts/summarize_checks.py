#!/usr/bin/env python3
"""summarize_checks.py - summary of checks/: rotation length vs harvest and carbon (calibrated model),
and the save/restart management-switch test (rs_bau must reproduce the uninterrupted run rot_r26_p0)."""
from pathlib import Path

import numpy as np
import pandas as pd

C = Path.home() / "lpjg_lishui" / "checks"


def s(run, f, col):
    p = C / run / "output" / f
    if not p.exists() or p.stat().st_size == 0:
        return pd.Series(dtype=float)
    try:
        df = pd.read_csv(p, sep=r"\s+")
    except pd.errors.EmptyDataError:
        return pd.Series(dtype=float)
    return df.groupby("Year")[col].mean() if col in df.columns else pd.Series(dtype=float)


rows = []
for R in (16, 21, 26, 31, 36, 41):
    m = []
    for k in range(4):
        r = f"rot_r{R}_p{k}"
        tot, veg = s(r, "cpool.out", "Total"), s(r, "cpool.out", "VegC")
        m.append(dict(TotC=tot.loc[1960:2019].mean(), VegC=veg.loc[1960:2019].mean(),
                      Harvest_yr=s(r, "cflux.out", "Harvest").loc[1960:2019].mean(),
                      NEE=s(r, "cflux.out", "NEE").loc[1960:2019].mean(),
                      NPP=s(r, "anpp.out", "TeNE").loc[1990:2019].mean()))
    rows.append({"rotation": R, **pd.DataFrame(m).mean().to_dict()})
pd.set_option("display.width", 200)
print("== Calibrated model: landscape means over 4 planting phases (kg C m-2; fluxes kg C m-2 yr-1) ==")
print(pd.DataFrame(rows).round(3).to_string(index=False))

runs = ["rot_r26_p0", "rs_save", "rs_bau", "rs_ext", "rs_stop"]
years = [1999, 2001, 2004, 2007, 2010, 2013, 2016, 2019]
print("\n== Save/restart test: vegetation carbon (kg C m-2); rot_r26_p0 = uninterrupted reference ==")
print(pd.DataFrame({r: s(r, "cpool.out", "VegC").reindex(years) for r in runs}, index=years).round(3).to_string())
ref, bau = s("rot_r26_p0", "cpool.out", "Total"), s("rs_bau", "cpool.out", "Total")
common = ref.index.intersection(bau.index)
common = common[common >= 2000]
if len(common):
    print(f"\nmax |rs_bau - reference| total carbon 2000-2019: {np.abs(bau[common] - ref[common]).max():.4f} "
          f"kg C m-2 over {len(common)} years")
else:
    print("\nno overlapping years after 2000 between rs_bau and the reference")
bl = {r: s(r, "cmass.out", "TeBE").reindex([2019]).fillna(0).iloc[0]
      + s(r, "cmass.out", "TeBS").reindex([2019]).fillna(0).iloc[0] for r in runs[2:]}
print("broadleaf carbon in 2019 (kg C m-2): " + ", ".join(f"{k} {v:.3f}" for k, v in bl.items()))
print("harvest 2000-2019 (kg C m-2): "
      + ", ".join(f"{r} {s(r, 'cflux.out', 'Harvest').loc[2000:2019].sum():.3f}" for r in runs[2:]))
