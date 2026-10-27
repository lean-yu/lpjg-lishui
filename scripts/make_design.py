#!/usr/bin/env python3
"""make_design.py - build the experiment matrix (design.csv) for the Lishui elevation x management study.

Edit the settings block, re-run, and commit design.csv together with the scripts.
Elevation bands are represented by a lapse-rate temperature offset relative to the
W5E5 grid-cell elevation (the elevation the 0.5-degree climate data represent).

Rows 1-288: main runs (pulse "none"). Then, for every main run and every entry of PULSES, a pulse run whose
climate is the same except for one extreme season (make_pulse_clim.py builds those files). The main run with
the same factors is the pulse run's control: the two are identical before the pulse year.
Each heat-parameter level has its own respcoeff, calibrated to the same growth targets (calib_growth.py).
"""
import csv
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# ---------------- settings ----------------
BASE = Path.home() / "lpjg_lishui"
CELL = (119.75, 28.25)               # lon, lat of the representative 0.5-degree cell (Lishui city)
LAPSE = 0.6                          # deg C per 100 m (standard environmental lapse rate)
ELEV = [200, 500, 800]               # m; valley / main Chinese fir belt / upper limit (docs/sources.md)
MGMT = ["bau26", "ext36", "stop"]    # from the switch year: 26-yr rotation / 36-yr (LY/T 2908) / no harvest
PSTEMP_HIGH = [25, 30]               # TeNE upper end of the photosynthesis temperature optimum (deg C)
RESPCOEFF = {25: 0.40, 30: 0.40}     # TeNE respcoeff calibrated for each pstemp_high level (calib_growth.py);
                                     # fill in the pstemp_high 30 value before generating the design
GCMS = ["gfdl-esm4", "ukesm1-0-ll"]
SSPS = ["ssp126", "ssp585"]
ROTATION_BEFORE_2025 = 26            # current practice, identical for all runs before the switch
PHASES = 4                           # staggered planting years -> landscape ("normal forest") mean
SWITCH_YEAR = 2024                   # startman2; set from the switch-year test (docs/LOG.md)
PULSES = [("hot3", 3.0, 1.0, 2045),  # (name, tas change in K, precipitation factor, calendar year)
          ("hot3", 3.0, 1.0, 2075),
          ("plc", 0.01, 1.0, 2045)]   # placebo: noise floor
PULSE_MONTHS = "7-9"                 # Jul-Sep, as in the earlier batch scenarios (make_scenarios.py)
# ------------------------------------------

missing = [p for p in PSTEMP_HIGH if RESPCOEFF.get(p) is None]
if missing:
    sys.exit(f"RESPCOEFF has no calibrated value for pstemp_high {missing}: enter it in the settings block of make_design.py")


def first_year(k):
    return 1901 + int(k * ROTATION_BEFORE_2025 / PHASES + 0.5)


def events(first, mgmt):
    """{year: "cc" or "thin"} of one stand: age counts from the last clearcut, thinning at age 9 and 13, and the
    rotation of the new management applies from SWITCH_YEAR on (matches the model output of the smoke runs)."""
    new_rot = {"bau26": 26, "ext36": 36, "stop": None}[mgmt]
    ev, last = {first: "cc"}, first
    for y in range(first + 1, 2101):
        rot = ROTATION_BEFORE_2025 if y < SWITCH_YEAR else new_rot
        age = y - last
        if rot and age >= rot:
            ev[y], last = "cc", y
        elif rot and age in (9, 13):
            ev[y] = "thin"
    return ev


def stand_age(first, mgmt, year):
    return year - max(y for y, kind in events(first, mgmt).items() if kind == "cc" and y <= year)


with Dataset(BASE / "data" / "orog_W5E5v2.0.nc") as nc:
    lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
    z = np.squeeze(np.asarray(nc.variables["orog"][:]))
    zcell = float(z[int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))])

rows = []
for e in ELEV:
    dT = round(-LAPSE * (e - zcell) / 100.0, 2)
    for m in MGMT:
        for p in PSTEMP_HIGH:
            for g in GCMS:
                for s in SSPS:
                    for k in range(PHASES):
                        rows.append(dict(
                            run_id=f"e{e}_{m}_t{p}_{g}_{s}_p{k}", elev_m=e, dT=dT, mgmt=m, pstemp_high=p,
                            respcoeff=RESPCOEFF[p], gcm=g, ssp=s, phase=k, firstmanageyear=first_year(k),
                            switch_year=SWITCH_YEAR, pulse="none", pulse_year="", pulse_dT=0.0, pulse_prf=1.0,
                            pulse_months=""))
main = list(rows)
for name, dt, prf, year in PULSES:
    tag = f"{name}y{year}"
    for r in main:
        rows.append(dict(r, run_id=f"{r['run_id']}_{tag}", pulse=tag, pulse_year=year, pulse_dT=dt,
                         pulse_prf=prf, pulse_months=PULSE_MONTHS))

out = BASE / "design.csv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

print(f"W5E5 elevation of cell {CELL}: {zcell:.0f} m")
for e in ELEV:
    print(f"  band {e:5d} m -> temperature offset {-LAPSE * (e - zcell) / 100.0:+.2f} deg C")
print(f"respcoeff by pstemp_high: {RESPCOEFF}; switch year {SWITCH_YEAR}")
for name, dt, prf, year in PULSES:
    stands = [(m, k) for m in MGMT for k in range(PHASES)]
    clash = [f"{m} p{k}" for m, k in stands if year in events(first_year(k), m)]
    ages = [stand_age(first_year(k), m, year) for m, k in stands]
    print(f"pulse {name}y{year}: tas {dt:+g} K, pr x{prf:g}, months {PULSE_MONTHS}; stand age {min(ages)}-{max(ages)} yr; "
          + (f"WARNING: harvest or thinning in the pulse year for {clash}" if clash else "no harvest or thinning that year"))
print(f"{len(main)} main runs + {len(rows) - len(main)} pulse runs = {len(rows)} rows written to {out}")
