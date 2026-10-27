#!/usr/bin/env python3
"""make_pulse_clim.py - climate files for the pulse runs in design.csv.

For every design row with a pulse, the row's control forcing is copied and changed only in the pulse months of
the pulse year: tas + pulse_dT (K), pr x pulse_prf. File names come from run_main.climate_paths, so the runner
reads exactly these files (data/climate3b/pulse/). Each file is written under a temporary name and renamed when
complete; existing files are kept (delete data/climate3b/pulse to rebuild).
Check printed per file: pulse minus control, pulse-month mean of the representative cell, for the year before,
the pulse year and the year after (expected 0 / change / 0).
Usage:  python make_pulse_clim.py
"""
import csv
import shutil
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_main import BASE, climate_paths  # noqa: E402  (same file names as the runner)

CELL = (119.75, 28.25)  # lon, lat of the simulated cell


def months_of(spec):
    """'7-9' -> [7, 8, 9]"""
    a, b = (int(x) for x in spec.split("-"))
    return list(range(a, b + 1))


def seasons(nc, years, months):
    """{year: time slice of the pulse months}; the months of each year must form one contiguous block of days."""
    t = nc.variables["time"]
    dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
    ym = np.array([(d.year, d.month) for d in dates])
    out = {}
    for y in years:
        idx = np.flatnonzero((ym[:, 0] == y) & np.isin(ym[:, 1], months))
        if len(idx) == 0 or idx[-1] - idx[0] + 1 != len(idx):
            raise ValueError(f"months {months} of {y} not found as one block of days")
        out[y] = slice(int(idx[0]), int(idx[-1]) + 1)
    return out


def make(src, dst, var, year, months, add, factor):
    if dst.exists():
        return "kept"
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".part")
    shutil.copyfile(src, tmp)
    with Dataset(tmp, "r+") as nc:
        v = nc.variables[var]
        if v.dimensions[0] != "time":
            raise ValueError(f"{src}: {var} dimensions {v.dimensions}, time must come first")
        sl = seasons(nc, [year], months)[year]
        v[sl, ...] = v[sl, ...] * factor + add
    tmp.rename(dst)
    return "made"


def cell_season_means(path, var, years, months):
    """Pulse-month mean of the simulated cell for each year."""
    with Dataset(path) as nc:
        lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
        i, j = int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))
        return {y: float(np.mean(nc.variables[var][sl, i, j])) for y, sl in seasons(nc, years, months).items()}


if __name__ == "__main__":
    rows = list(csv.DictReader(open(BASE / "design.csv")))
    if "pulse" not in rows[0]:
        sys.exit("design.csv has no pulse columns: re-run scripts/make_design.py")
    jobs = {}
    for r in rows:
        if r["pulse"] == "none":
            continue
        base, new = climate_paths(dict(r, pulse="none")), climate_paths(r)
        year, months = int(r["pulse_year"]), months_of(r["pulse_months"])
        for key, var, add, factor in (("file_temp", "tas", float(r["pulse_dT"]), 1.0),
                                      ("file_prec", "pr", 0.0, float(r["pulse_prf"]))):
            if new[key] != base[key]:
                jobs[new[key]] = (base[key], var, year, months, add, factor)
    print(f"{len(jobs)} pulse climate files for {sum(r['pulse'] != 'none' for r in rows)} pulse runs")
    for dst, (src, var, year, months, add, factor) in sorted(jobs.items()):
        status = make(src, dst, var, year, months, add, factor)
        years = (year - 1, year, year + 1)
        a, b = cell_season_means(src, var, years, months), cell_season_means(dst, var, years, months)
        diffs = [f"{b[y] - a[y]:+.2f}" if var == "tas" else f"x{b[y] / a[y]:.2f}" for y in years]
        print(f"{status} {dst.name}: {var} months {min(months)}-{max(months)}, {year - 1}/{year}/{year + 1}: {' '.join(diffs)}")
