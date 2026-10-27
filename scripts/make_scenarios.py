#!/usr/bin/env python3
# make_scenarios.py：在对照气候基础上制造 7–9 月的干旱、热浪情景
# 用法：python make_scenarios.py [年份 ...]   默认 1993 2000 2008（1988 年造林后第 5、12、20 年）
import shutil
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

BASE = Path.home() / "lpjg_lishui"
CLIM = BASE / "data" / "climate"
OUT = BASE / "data" / "climate_scen"
SRC = {"pr": CLIM / "pr_lishui_daily_1901_2019.nc", "tas": CLIM / "tas_lishui_daily_1901_2019.nc"}
MONTHS = (7, 8, 9)
YEARS = [int(a) for a in sys.argv[1:]] or [1993, 2000, 2008]


def scenarios():
    s = {}
    for y in YEARS:
        s[f"dry50_{y}"] = dict(pr=0.5, dt=0.0, years=[y])            # 7–9 月降水减半
        s[f"dry50x2_{y}"] = dict(pr=0.5, dt=0.0, years=[y, y + 1])   # 连续两年减半
        s[f"hot3_{y}"] = dict(pr=1.0, dt=3.0, years=[y])             # 7–9 月升温 3 °C
        s[f"comp_{y}"] = dict(pr=0.5, dt=3.0, years=[y])             # 减半 + 升温
    return s


def season_index(path, years):
    with Dataset(path) as nc:
        t = nc.variables["time"]
        dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
        return np.where([(d.year in years) and (d.month in MONTHS) for d in dates])[0]


def write(var, dst, idx, factor=1.0, add=0.0):
    shutil.copyfile(SRC[var], dst)
    if factor != 1.0 or add != 0.0:
        with Dataset(dst, "r+") as nc:
            v = nc.variables[var]
            v[idx, ...] = v[idx, ...] * factor + add


def season_stat(path, var, idx, nyears):
    with Dataset(path) as nc:
        x = np.asarray(nc.variables[var][idx, ...], dtype=float).mean(axis=(1, 2))
    return x.sum() * 86400 / nyears if var == "pr" else x.mean() - 273.15


for name, s in scenarios().items():
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    idx = season_index(SRC["pr"], s["years"])
    write("pr", d / "pr.nc", idx, factor=s["pr"])
    write("tas", d / "tas.nc", idx, add=s["dt"])
    n = len(s["years"])
    p0, p1 = season_stat(SRC["pr"], "pr", idx, n), season_stat(d / "pr.nc", "pr", idx, n)
    t0, t1 = season_stat(SRC["tas"], "tas", idx, n), season_stat(d / "tas.nc", "tas", idx, n)
    print(f"{name:14s} 7-9月降水 {p0:6.0f} -> {p1:6.0f} mm/年   气温 {t0:5.1f} -> {t1:5.1f} °C")
print(f"完成：{OUT}")
