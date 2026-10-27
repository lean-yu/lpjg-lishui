#!/usr/bin/env python3
import os
import glob
import csv

BASE = os.path.dirname(os.path.abspath(__file__))
OUT_FILES = {
    'cmass': 'cmass.out',
    'agpp': 'agpp.out',
    'anpp': 'anpp.out',
    'lai': 'lai.out',
    'cflux': 'cflux.out',
    'fpc': 'fpc.out',
    'nflux': 'nflux.out',
    'aaet': 'aaet.out',
}

def extract_last_total(path):
    """返回 (lon, lat, year, total) 或 None"""
    try:
        with open(path) as f:
            lines = [l for l in f if l.strip()]
        if len(lines) < 2:
            return None
        header = lines[0].split()
        # 找 Total 列索引
        if 'Total' in header:
            tidx = header.index('Total')
        else:
            tidx = len(header) - 1
        last = lines[-1].split()
        if len(last) <= tidx:
            return None
        lon = last[0]; lat = last[1]; year = last[2]; total = last[tidx]
        return (lon, lat, year, total)
    except Exception:
        return None

rows = []
output_dirs = []
for root, dirs, files in os.walk(BASE):
    if os.path.basename(root) == 'output':
        output_dirs.append(root)

for od in sorted(output_dirs):
    run_root = os.path.dirname(od)
    run_id = os.path.relpath(run_root, BASE)
    row = {'run_id': run_id}
    lon = lat = year = None
    for key, fname in OUT_FILES.items():
        fp = os.path.join(od, fname)
        if os.path.exists(fp):
            res = extract_last_total(fp)
            if res:
                lon, lat, year, total = res
                row[key + '_total'] = total
        else:
            row[key + '_total'] = ''
    row['lon'] = lon or ''
    row['lat'] = lat or ''
    row['year'] = year or ''
    rows.append(row)

fieldnames = ['run_id', 'lon', 'lat', 'year'] + [k + '_total' for k in OUT_FILES]
out_csv = os.path.join(BASE, 'summary_results.csv')
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows:
        w.writerow(r)

print(f'汇总完成：{len(rows)} 个运行目录 -> {out_csv}')
