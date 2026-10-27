#!/usr/bin/env bash
# next_steps.sh - (1) refresh the R export, (2) CO2 attribution runs and their analysis, (3) redraw the pulse figures.
# Started in the background by the paste block; log: ~/lpjg_lishui/next_steps.log
cd "$HOME/lpjg_lishui" || exit 1
PY=/usr/bin/python3                                   # system Python: pandas, matplotlib, scipy
$PY -c "import pandas, matplotlib, scipy" 2>/dev/null || PY="$HOME/lpjg_venv/bin/python"
VPY="$HOME/lpjg_venv/bin/python"                      # venv: netCDF4 for export_for_r.py
echo "== 1/3 export_for_r.py   $(date +%T)"
$VPY scripts/export_for_r.py; echo "-- exit $?"
echo "== 2/3 co2fix.py: 144 fixed-CO2 runs (about 30 min), then the analysis   $(date +%T)"
$PY scripts/co2fix.py; echo "-- exit $?"
echo "== 3/3 analyze_pulse.py --replot (figures with units)   $(date +%T)"
$PY scripts/analyze_pulse.py --replot > pulse/replot.log 2>&1; echo "-- exit $? (tables in pulse/replot.log)"
echo "== all done   $(date +%T)"
