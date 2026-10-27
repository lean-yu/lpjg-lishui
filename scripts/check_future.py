#!/usr/bin/env python3
"""check_future.py - integrity check of the downloaded ISIMIP3b climate and CO2 files.

For each merged file data/climate3b/<gcm>_<ssp>_<var>_1901_2100.nc: time coverage, missing years,
gaps or duplicates in the daily axis, missing values at the Lishui cell and in the box, value ranges,
and the climate signal. Also counts the pieces in each downloaded zip and checks the CO2 files.
With --fix: removes incomplete zips, faulty merged files (and the zips they were built from) and
incomplete CO2 files, so that re-running get_future.sh downloads and merges them again.
"""
import sys
import zipfile
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BASE = Path.home() / "lpjg_lishui"
RAW, OUT, CO2D = BASE / "data" / "raw3b", BASE / "data" / "climate3b", BASE / "data" / "co2"
CELL = (119.75, 28.25)
UNIT_DAYS = {"days": 1.0, "hours": 1 / 24, "minutes": 1 / 1440, "seconds": 1 / 86400}
CONVERT = {"tas": (lambda x: x - 273.15, "degC"), "pr": (lambda x: x * 86400 * 365.25, "mm/yr"),
           "rsds": (lambda x: x, "W/m2")}


def to_float(a):
    a = np.ma.asarray(a).astype(float).filled(np.nan)
    a[np.abs(a) > 1e19] = np.nan
    return a


def ranges(years):
    """[1931, 1932, ..., 1940, 2051] -> '1931-1940, 2051'"""
    out, years = [], sorted(years)
    for y in years:
        if out and y == out[-1][1] + 1:
            out[-1][1] = y
        else:
            out.append([y, y])
    return ", ".join(f"{a}-{b}" if a != b else f"{a}" for a, b in out) or "none"


def check_nc(path):
    var = path.name.split("_")[2]
    with Dataset(path) as nc:
        t = nc.variables["time"]
        units, cal = t.units, getattr(t, "calendar", "standard")
        tv = np.asarray(t[:], dtype=float)
        years = np.array([d.year for d in num2date(tv, units, cal)])
        step = np.diff(tv) * UNIT_DAYS.get(units.split()[0], 1.0)
        lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
        j, i = int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))
        cell = to_float(nc.variables[var][:, j, i])
        box = to_float(nc.variables[var][0, :, :])
    f, unit = CONVERT.get(var, (lambda x: x, "?"))
    mean = lambda y0, y1: f(np.nanmean(cell[(years >= y0) & (years <= y1)])) if ((years >= y0) & (years <= y1)).any() else np.nan
    print(f"{path.name}: {years[0]}-{years[-1]}, {len(tv)} days, calendar {cal}")
    print(f"   missing years: {ranges(set(range(1901, 2101)) - set(years.tolist()))}; "
          f"gaps {int((step > 1.01).sum())}, duplicates/backwards {int((step < 0.99).sum())}")
    print(f"   Lishui cell: {int(np.isnan(cell).sum())} missing values (years: {ranges(set(years[np.isnan(cell)].tolist()))}); "
          f"box: {int(np.isnan(box).sum())} of {box.size} cells missing on day 1")
    print(f"   cell mean [{unit}]: 1981-2010 {mean(1981, 2010):.1f}, 2071-2100 {mean(2071, 2100):.1f}; "
          f"join check 2005-2014 {mean(2005, 2014):.1f} vs 2015-2024 {mean(2015, 2024):.1f}")
    return len(set(range(1901, 2101)) - set(years.tolist())) == 0 and not (step > 1.01).any() \
        and not (step < 0.99).any() and not np.isnan(cell).any()


if __name__ == "__main__":
    fix = "--fix" in sys.argv
    bad_zips, bad_merged, bad_co2 = set(), set(), set()
    print("== downloaded pieces (expected: historical 12, ssp 9) ==")
    for z in sorted(RAW.glob("*.zip")):
        try:
            n = sum(1 for m in zipfile.ZipFile(z).namelist() if m.endswith(".nc"))
        except zipfile.BadZipFile:
            n = -1
        exp = 12 if "_historical_" in z.name else 9
        print(f"   {z.name}: {n} files" + ("" if n == exp else f"   <-- expected {exp}"))
        if n != exp:
            bad_zips.add(z)
    print("\n== merged files ==")
    files = sorted(OUT.glob("*_1901_2100.nc"))
    if not files:
        print("   none found")
    for p in files:
        try:
            ok = check_nc(p)
        except Exception as e:  # report and continue
            print(f"{p.name}: cannot read ({e})")
            ok = False
        if not ok:
            bad_merged.add(p)
    print("\n== CO2 files ==")
    for p in sorted(CO2D.glob("co2_*")):
        lines = p.read_text(errors="ignore").splitlines()
        data = [ln.split() for ln in lines if ln.strip() and ln.split()[0].isdigit()]
        yrs = [int(r[0]) for r in data]
        blank = sum(1 for ln in lines if not ln.strip())
        span = f"{min(yrs)}-{max(yrs)}" if yrs else "no data"
        gaps = ranges(set(range(min(yrs), max(yrs) + 1)) - set(yrs)) if yrs else "-"
        print(f"   {p.name}: {len(data)} years ({span}), missing years {gaps}, blank lines {blank}")
        if p.name.endswith("_1850_2100_lpjg.txt") and (not yrs or gaps != "none" or min(yrs) > 1850 or max(yrs) < 2100):
            bad_co2.add(p)
    log = BASE / "future.log"
    if log.exists():
        bad = [ln for ln in log.read_text(errors="ignore").splitlines()
               if any(k in ln for k in ("failed", "missing", "error", "Error"))]
        print("\n== problems reported in future.log ==")
        print("\n".join("   " + ln for ln in bad[-20:]) or "   none")
    if fix:
        print("\n== --fix ==")
        for p in bad_merged:
            g, s_, v = p.name.split("_")[:3]
            src = {RAW / f"{g}_historical_{v}.zip", RAW / f"{g}_{s_}_{v}.zip"}
            if not (src & bad_zips):          # zips look complete, so the pieces themselves are suspect
                bad_zips.update(z for z in src if z.exists())
        for z in sorted(bad_zips):
            if z.exists():
                z.unlink()
                print(f"   removed {z.name}")
        for p in files:
            g, s_, v = p.name.split("_")[:3]
            if p.exists() and (p in bad_merged or {RAW / f"{g}_historical_{v}.zip", RAW / f"{g}_{s_}_{v}.zip"} & bad_zips):
                p.unlink()
                print(f"   removed {p.name}")
        for p in bad_co2:
            p.unlink()
            print(f"   removed {p.name}")
        print("   now re-run:  cd ~/lpjg_lishui && nohup bash scripts/get_future.sh > future.log 2>&1 &")
