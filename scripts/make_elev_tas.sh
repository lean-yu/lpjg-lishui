#!/usr/bin/env bash
# make_elev_tas.sh - elevation-band temperature forcing: tas + lapse-rate offset (elev_m, dT from design.csv)
B=$HOME/lpjg_lishui; IN=$B/data/climate3b; OUT=$IN/elev; mkdir -p $OUT
tail -n +2 $B/design.csv | cut -d, -f2,3 | sort -u | while IFS=, read -r elev dT; do
  for f in $IN/*_tas_1901_2100.nc; do
    o=$OUT/$(basename "${f%_1901_2100.nc}")_e${elev}_1901_2100.nc
    [ -s "$o" ] || cdo -s -O addc,$dT "$f" "$o"
    echo "$(basename "$o"): offset $dT C"
  done
done
