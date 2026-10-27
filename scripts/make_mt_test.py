#!/usr/bin/env python3
"""make_mt_test.py - test switching forest management in a given year with management types (mt blocks).

Builds, from checks/rot_r26_p0 (calibrated fir, 26-yr rotation, planted 1901):
  mt_ref   unchanged configuration (reference)
  mt_bau   switch in SWITCH_YEAR to an identical copy of the 26-yr management  -> must equal mt_ref
  mt_ext   switch in SWITCH_YEAR to a 36-yr rotation (thinning at age 9 and 13)
  mt_stop  switch in SWITCH_YEAR to no harvest
  nd_hist  unchanged configuration with ndep_timeseries "historic" (N deposition sensitivity)
All runs use ndep_timeseries "rcp60" unless stated. Runs are executed 5 in parallel.
Mechanism (management.cpp): a forest stand with >= 2 management types rotates to management m in
calendar year startman<m> when rot_wait_for_cc is 0.
"""
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = Path.home() / "lpjg_lishui"
C = BASE / "checks"
TEMPLATE = C / "rot_r26_p0"
GUESS = (BASE / ".guess2_bin").read_text().strip()
SWITCH_YEAR = 2000
ST_NAME = "Fir_CC"
ST_ONLY = ("stinclude", "landcover", "naturalveg", "reestab", "intercrop", "multicrop", "firstrotyear",
           "distinterval", "rot_wait_for_cc", "management1", "management2", "management3",
           "startman1", "startman2", "startman3")


def st_block(lines):
    """Index of the first and last line (closing parenthesis) of st "Fir_CC"."""
    i0 = next(i for i, ln in enumerate(lines) if re.match(rf'\s*st\s+"{ST_NAME}"', ln))
    i1 = next(i for i in range(i0 + 1, len(lines)) if lines[i].strip().startswith(")"))
    return i0, i1


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


def build(name, second=None, ndep="rcp60"):
    d = C / name
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(TEMPLATE, d, ignore=shutil.ignore_patterns("output", "run.log", "state"))
    (d / "output").mkdir()
    cf = d / "lishui_cf.ins"
    cf.write_text("\n".join(ln for ln in cf.read_text().splitlines() if not ln.startswith("ndep_timeseries"))
                  + f'\nndep_timeseries "{ndep}"\n')
    if second is None:
        return d
    fir = d / "lishui_fir.ins"
    lines = fir.read_text().splitlines()
    i0, i1 = st_block(lines)
    mgmt = [ln for ln in lines[i0 + 1:i1]
            if ln.split() and not ln.strip().startswith("!") and ln.split()[0] not in ST_ONLY
            and len(ln.split()) > 1]                       # keyword + value lines (group includes are skipped)
    m1 = ["\tmtinclude 1"] + mgmt
    m2 = ["\tmtinclude 1"] + mgmt
    m2 = set_key(set_key(m2, "cutfirstyear", 0), "cutfirstyear_nonsel", 0)   # never clearcut at the switch
    if second == "ext36":
        m2 = rotation_lines(m2, 36)
    elif second == "stop":
        m2 = set_key(set_key(set_key(m2, "cutinterval", 1000), "thinstr", "0 0 0 0 0"), "thinstrunsel", "0 0 0 0 0")
    mt_text = [f'mt "M1" (', *m1, ")", f'mt "M2_{second}" (', *m2, ")", ""]
    st_extra = ['\tmanagement1 "M1"', f'\tmanagement2 "M2_{second}"', f"\tstartman2 {SWITCH_YEAR}",
                "\trot_wait_for_cc 0"]
    new = lines[:i0] + mt_text + lines[i0:i1] + st_extra + lines[i1:]
    fir.write_text("\n".join(new) + "\n")
    return d


def run(d):
    with open(d / "run.log", "w") as log:
        subprocess.run([GUESS, "-input", "cf", "lishui_fir.ins"], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    return d.name, "Finished" in (d / "run.log").read_text(errors="ignore")


if __name__ == "__main__":
    dirs = [build("mt_ref"), build("mt_bau", "copy26"), build("mt_ext", "ext36"), build("mt_stop", "stop"),
            build("nd_hist", ndep="historic")]
    with ThreadPoolExecutor(max_workers=5) as ex:
        for name, ok in ex.map(run, dirs):
            print(f"{name}: {'finished' if ok else 'FAILED, see run.log'}")
