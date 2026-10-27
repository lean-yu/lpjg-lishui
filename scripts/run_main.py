#!/usr/bin/env python3
"""run_main.py - main experiment: elevation band x management from 2025 x heat parameter x GCM x SSP x phase,
plus heat-pulse runs (same factors, one extreme season added to the climate).

Reads design.csv (scripts/make_design.py). Each run is one uninterrupted simulation 1901-2100:
  climate   ISIMIP3b GCM historical+SSP (tas shifted to the elevation band, pr and rsds unchanged);
            pulse runs read the modified copies made by make_pulse_clim.py (file names: climate_paths)
  CO2       co2_<ssp>_1850_2100_lpjg.txt;  N deposition: ndep_timeseries "rcp60" for all runs (docs/LOG.md)
  fir       calibrated parameters (scripts/tene_calibrated.ins), then pstemp_high and respcoeff from the design
            (each pstemp_high level has its own calibrated respcoeff; the last value in the block applies)
  stand     2500 stems/ha, 26-yr rotation with thinning at age 9 and 13, first management year = phase year;
            in the design's switch_year the stand switches (mt blocks, startman2) to bau26 (identical copy),
            ext36 or stop.
Usage:  python run_main.py [--only N] [--match REGEX] [--workers N] [--dry]
        --only: first N design rows; --match: rows whose run_id matches REGEX (e.g. '_p[0-3]$' = main runs only,
        'hot3' = pulse runs only); --workers: parallel runs (default 4); --dry: build without running
Re-runnable: runs whose run.log contains "Finished" are skipped. Output: runs_main/<run_id>/
"""
import argparse
import csv
import datetime
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

BASE = Path.home() / "lpjg_lishui"
TEMPLATE = BASE / "exp_base"
OUT = BASE / "runs_main"
CLIM = BASE / "data" / "climate3b"
PULSE_DIR = CLIM / "pulse"
GUESS = (BASE / ".guess2_bin").read_text().strip()
TENE_FILE = BASE / "scripts" / "tene_calibrated.ins"
ROTATION = 26
ST_NAME = "Fir_CC"
REQUIRED = ("respcoeff", "switch_year", "pulse", "pulse_dT", "pulse_prf")   # design.csv from make_design.py
ST_ONLY = ("stinclude", "landcover", "naturalveg", "reestab", "intercrop", "multicrop", "firstrotyear",
           "distinterval", "rot_wait_for_cc", "management1", "management2", "management3",
           "startman1", "startman2", "startman3")


def git_commit():
    return subprocess.run(["git", "-C", str(BASE), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip()


def set_param(text, name, path):
    """Replace the string of param "name" (str "...") in a cf instruction file; must occur exactly once."""
    pat = re.compile(rf'^(\s*param\s+"{name}"\s*\(\s*str\s+)"[^"]*"', re.M)
    new, n = pat.subn(rf'\g<1>"{path}"', text)
    if n != 1:
        raise ValueError(f'param "{name}" found {n} times')
    return new


def set_key(block, key, value):
    """Replace the value of key in a list of lines, or append it."""
    for i, ln in enumerate(block):
        if ln.split() and ln.split()[0] == key:
            block[i] = f"\t{key} {value}"
            return block
    return block + [f"\t{key} {value}"]


def rotation_lines(block, rotation):
    block = set_key(block, "cutinterval", rotation)
    return set_key(block, "thintime", f"{9.5 / rotation:.4f} {13.5 / rotation:.4f} 0 0 0")


def tene_lines(row):
    lines = [ln.strip() for ln in TENE_FILE.read_text().splitlines() if ln.strip() and not ln.strip().startswith("!")]
    return ([f"\t{ln}   ! calibrated (tene_calibrated.ins)" for ln in lines]
            + [f"\tpstemp_high {row['pstemp_high']}   ! design",
               f"\trespcoeff {row['respcoeff']}   ! design: calibrated for this pstemp_high"])


def climate_paths(row):
    """Climate files of a design row. Pulse runs read copies changed in the pulse season (make_pulse_clim.py)."""
    g, s, e, pulse = row["gcm"], row["ssp"], row["elev_m"], row.get("pulse", "none")
    p = {"file_temp": CLIM / "elev" / f"{g}_{s}_tas_e{e}_1901_2100.nc",
         "file_prec": CLIM / f"{g}_{s}_pr_1901_2100.nc",
         "file_insol": CLIM / f"{g}_{s}_rsds_1901_2100.nc"}
    if pulse != "none":
        if float(row["pulse_dT"]) != 0.0:
            p["file_temp"] = PULSE_DIR / f"{g}_{s}_tas_e{e}_{pulse}_1901_2100.nc"
        if float(row["pulse_prf"]) != 1.0:
            p["file_prec"] = PULSE_DIR / f"{g}_{s}_pr_{pulse}_1901_2100.nc"
    return p


def insert_in_block(text, header_regex, new_lines):
    """Insert lines before the closing parenthesis of the first block whose header matches."""
    out, inb, done = [], False, False
    for ln in text.splitlines():
        if not done and re.match(header_regex, ln):
            inb = True
        if inb and ln.strip().startswith(")"):
            out.extend(new_lines)
            inb, done = False, True
        out.append(ln)
    if not done:
        raise ValueError(f"block {header_regex} not found")
    return "\n".join(out) + "\n"


def stand_file(text, first_year, mgmt, switch_year):
    lines = text.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if re.match(rf'\s*st\s+"{ST_NAME}"', ln))
    i1 = next(i for i in range(i0 + 1, len(lines)) if lines[i].strip().startswith(")"))
    st = lines[i0 + 1:i1]
    st = set_key(st, "firstmanageyear", first_year)
    st = rotation_lines(st, ROTATION)
    st = set_key(st, "thinstr", "0.28 0.24 0 0 0")
    st = set_key(st, "thinstrunsel", "0.95 0 0 0 0")
    st = set_key(st, "thinselectdiam", "1 1 0 0 0")
    mgmt_lines = [ln for ln in st if ln.split() and not ln.strip().startswith("!")
                  and ln.split()[0] not in ST_ONLY and len(ln.split()) > 1]
    m1 = ["\tmtinclude 1"] + mgmt_lines
    m2 = set_key(set_key(["\tmtinclude 1"] + list(mgmt_lines), "cutfirstyear", 0), "cutfirstyear_nonsel", 0)
    if mgmt == "ext36":
        m2 = rotation_lines(m2, 36)
    elif mgmt == "stop":
        m2 = set_key(set_key(set_key(m2, "cutinterval", 1000), "thinstr", "0 0 0 0 0"), "thinstrunsel", "0 0 0 0 0")
    elif mgmt != "bau26":
        raise ValueError(f"unknown management {mgmt}")
    mt = ['mt "M1" (', *m1, ")", f'mt "M2_{mgmt}" (', *m2, ")", ""]
    st += ['\tmanagement1 "M1"', f'\tmanagement2 "M2_{mgmt}"', f"\tstartman2 {switch_year}", "\trot_wait_for_cc 0"]
    return "\n".join(lines[:i0] + mt + [lines[i0]] + st + lines[i1:]) + "\n"


def build(row, commit):
    d = OUT / row["run_id"]
    if (d / "run.log").exists() and "Finished" in (d / "run.log").read_text(errors="ignore"):
        return None
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(TEMPLATE, d, ignore=shutil.ignore_patterns("output", "run.log", "state"))
    (d / "output").mkdir()
    paths = {**climate_paths(row), "file_co2": BASE / "data" / "co2" / f"co2_{row['ssp']}_1850_2100_lpjg.txt"}
    for p in paths.values():
        if not p.exists():
            raise FileNotFoundError(f"{p} (pulse climate files: python scripts/make_pulse_clim.py)")
    cf = (d / "lishui_cf.ins").read_text()
    for k, p in paths.items():
        cf = set_param(cf, k, p)
    cf = "\n".join(ln for ln in cf.splitlines() if not ln.startswith("ndep_timeseries")) + '\nndep_timeseries "rcp60"\n'
    (d / "lishui_cf.ins").write_text(cf)
    gl = (d / "lishui_global.ins").read_text()
    (d / "lishui_global.ins").write_text(insert_in_block(gl, r'\s*pft\s+"TeNE"', tene_lines(row)))
    (d / "lishui_fir.ins").write_text(stand_file((d / "lishui_fir.ins").read_text(), row["firstmanageyear"], row["mgmt"],
                                                        row["switch_year"]))
    info = dict(row, git=commit, created=datetime.datetime.now().isoformat(timespec="seconds"))
    (d / "run_info.txt").write_text("".join(f"{k}={v}\n" for k, v in info.items()))
    return d


def run(d):
    with open(d / "run.log", "w") as log:
        subprocess.run([GUESS, "-input", "cf", "lishui_fir.ins"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    return d.name, "Finished" in (d / "run.log").read_text(errors="ignore")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int, default=0, help="use only the first N design rows")
    ap.add_argument("--match", default="", help="use only rows whose run_id matches this regular expression")
    ap.add_argument("--workers", type=int, default=4, help="number of runs in parallel")
    ap.add_argument("--dry", action="store_true", help="build run directories without running")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(BASE / "design.csv")))
    lacking = [c for c in REQUIRED if c not in rows[0]]
    if lacking:
        raise SystemExit(f"design.csv lacks columns {lacking}: re-run scripts/make_design.py")
    rows = rows[:a.only] if a.only else rows
    rows = [r for r in rows if re.search(a.match, r["run_id"])]
    commit = git_commit()
    todo = [d for d in (build(r, commit) for r in rows) if d is not None]
    print(f"{len(rows)} design rows, {len(todo)} to run (git {commit}), start {datetime.datetime.now():%H:%M:%S}", flush=True)
    if a.dry:
        raise SystemExit(0)
    failed = []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run, d) for d in todo]
        for i, f in enumerate(as_completed(futs), 1):
            name, ok = f.result()
            print(f"{datetime.datetime.now():%H:%M:%S} {i}/{len(todo)} {name} {'ok' if ok else 'FAILED'}", flush=True)
            if not ok:
                failed.append(name)
    print(f"all finished {datetime.datetime.now():%H:%M:%S}; failed: {len(failed)} {failed[:5]}")
