#!/usr/bin/env python3
# summarize_rotland.py - landscape means over planting phases for each rotation and climate
from pathlib import Path
import pandas as pd

RL = Path.home() / "lpjg_lishui" / "rotland"


def s(run, f, col):
    df = pd.read_csv(RL / run / "output" / f, sep=r"\s+")
    return df.groupby("Year")[col].mean()


def metrics(run):
    tot, veg, prod = s(run, "cpool.out", "Total"), s(run, "cpool.out", "VegC"), s(run, "cpool.out", "HarvSlowC")
    return dict(NPP=s(run, "anpp.out", "TeNE").loc[1990:2019].mean(),
                TotC=tot.loc[1960:2019].mean(), VegC=veg.loc[1960:2019].mean(), Prod=prod.loc[1960:2019].mean(),
                SoilLit=(tot - veg - prod).loc[1960:2019].mean(),
                NEE=s(run, "cflux.out", "NEE").loc[1960:2019].mean(),
                Harv=s(run, "cflux.out", "Harvest").loc[1960:2019].mean(),
                TotC2019=tot.loc[2019])


rows = []
for R in (16, 21, 26, 31, 41, 500):
    ks = [0] if R == 500 else [0, 1, 2, 3]
    m = {sc: pd.DataFrame([metrics(f"r{R}_o{k}_{sc}") for k in ks]).mean()
         for sc in ("ctrl", "warm2_9019", "dry30_9019")}
    c, w, d = m["ctrl"], m["warm2_9019"], m["dry30_9019"]
    rows.append({"rotation": "none" if R == 500 else R, "NPP": c.NPP, "TotC": c.TotC, "VegC": c.VegC,
                 "SoilLitC": c.SoilLit, "Products": c.Prod, "NEE": c.NEE, "Harvest_yr": c.Harv,
                 "warm_dNPP%": 100 * (w.NPP / c.NPP - 1), "dry_dNPP%": 100 * (d.NPP / c.NPP - 1),
                 "warm_dTotC2019": w.TotC2019 - c.TotC2019, "dry_dTotC2019": d.TotC2019 - c.TotC2019})
pd.set_option("display.width", 200)
print("Landscape means over 4 planting phases (stocks kg C m-2, fluxes kg C m-2 yr-1; x10 = t C/ha)")
print(pd.DataFrame(rows).round(3).to_string(index=False))
