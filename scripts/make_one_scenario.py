#!/usr/bin/env python3
# 生成一个自定义情景：python make_one_scenario.py 名称 降水倍数 升温(K) 月份(逗号分隔) 年份(逗号分隔)
# 例：python make_one_scenario.py dry80_2000 0.2 0 7,8,9 2000
import shutil
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset, num2date

name, fac, dt = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
months = [int(m) for m in sys.argv[4].split(",")]
years = [int(y) for y in sys.argv[5].split(",")]
BASE = Path.home() / "lpjg_lishui"
SRC = {v: BASE / "data" / "climate" / f"{v}_lishui_daily_1901_2019.nc" for v in ("pr", "tas")}
d = BASE / "data" / "climate_scen" / name
d.mkdir(parents=True, exist_ok=True)
with Dataset(SRC["pr"]) as nc:
    t = nc.variables["time"]
    dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
idx = np.where([(x.year in years) and (x.month in months) for x in dates])[0]
for var, f, a in (("pr", fac, 0.0), ("tas", 1.0, dt)):
    shutil.copyfile(SRC[var], d / f"{var}.nc")
    if f != 1.0 or a != 0.0:
        with Dataset(d / f"{var}.nc", "r+") as nc:
            v = nc.variables[var]
            v[idx, ...] = v[idx, ...] * f + a
with Dataset(SRC["pr"]) as a, Dataset(d / "pr.nc") as b:
    p0 = np.asarray(a.variables["pr"][idx, ...], float).mean(axis=(1, 2)).sum() * 86400 / len(years)
    p1 = np.asarray(b.variables["pr"][idx, ...], float).mean(axis=(1, 2)).sum() * 86400 / len(years)
print(f"{name}: 月份 {months} 年份 {years} 降水 {p0:.0f} -> {p1:.0f} mm/年，升温 {dt:g} °C")
