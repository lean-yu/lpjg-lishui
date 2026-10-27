#!/usr/bin/env python3
"""analyze_pulse.py - paired heat-pulse metrics: every pulse run against its control, the main run with the same
factors. Both runs are identical before the pulse year (deterministic model; checked below), so the differences
from the pulse year on are the pulse effect.

Per pulse run (Y = pulse year; differences = pulse minus control, kg C m-2 unless stated):
  age              stand age in Y: years since the control's last clearcut (dens.out TeNE = 0 with harvest)
  vegc_pre         control VegC in Y-1
  dNPP, dNPP_pct   NPP change in Y (absolute, %); dNPP_pct is NaN when the control NPP rounds to 0
  dNEP             ecosystem NEP change in Y, NEP = -(NEE-LU): without harvest and wood-product fluxes
  dNPP_post_pct    NPP change over Y+1..Y+3 (%)
  dVegC_5y         VegC change in Y+5
  dEcoC_10y        ecosystem carbon change (cpool Total-slow) in Y+10
  dMort_3y         change in tree_mort_sts.out (Fir_CC) summed over Y..Y+2
  rec_VegC         years until VegC is back within 1 % of the control, searched up to the control's next clearcut
                   (0 = never more than 1 % below; NaN = not recovered before that clearcut or 2100)
  dRemoval_2100    change in harvest removals (Harvest + product-pool inflow) summed over Y..2100
  dRemoval_10y     change in harvest removals summed over Y..Y+10
  dUptake_10y      change in NEP summed over Y..Y+10 = dEcoC_10y + dRemoval_10y: carbon the stand did not take up,
                   wherever it ended up. Use it to compare managements: where the control is clearcut or thinned
                   in the window, part of the loss leaves as lower removals instead of staying in dEcoC_10y
  harv_events_10y  years in Y..Y+10 in which the control exports harvest (Harvest flux > 0); a thinning that leaves
                   the cut trees on site (the age-9 thinning in the current setup) is not counted
  yrs_to_cc        years from Y to the control's next clearcut (NaN: none by 2100). A clearcut soon after the
                   pulse removes the stand before it can make up the loss, so the loss stays (see the timing table)
Flux decomposition of the pulse effect (columns appended; the budget identities are checked in the report):
  NEP = NPP - Rh - other, NPP = GPP - Ra, with Rh = cflux.out Soil, other = Fire + Repr + Est + Seed (the cflux.out
  columns present) and GPP = agpp.out Total. In the pulse year Y:
  dGPP, dGPP_pct   GPP change (NaN when the runs have no agpp.out)
  dRa, dRa_pct     autotrophic respiration change, Ra = GPP - NPP (NaN without agpp.out)
  dRh, dRh_pct     heterotrophic respiration change
  dOther           change in the other terms
  dGPP_pm, dRa_pm, dRh_pm   the same three over the pulse months of Y only (mgpp.out, mra.out, mrh.out; NaN without
                   them): the direct response; the annual value minus this is the rest of the year
  Summed over Y..Y+10 (dUptake_10y = dNPP_10y - dRh_10y - dOther_10y): dGPP_10y, dRa_10y, dNPP_10y, dRh_10y,
  dOther_10y; summed over Y+1..Y+10 only (the legacy after the pulse year): dNPP_post10, dRh_post10
  dMort_10y        change in tree_mort_sts.out (Fir_CC) summed over Y..Y+10, in the units of that file (a carbon mass
                   compares directly with dLitterIn_10y: the mortality part of the extra litter)
  dLitterIn_10y    change in the carbon moved from vegetation to litter over Y..Y+10 (turnover, mortality and harvest
                   residues), from the vegetation budget: sum(NPP - Repr - Est - Seed) - dVegC_10y - dRemoval_10y
  dVegC_10y, dLitSoilC_10y   VegC and LitterC + SoilC change in Y+10; dEcoC_10y = dVegC_10y + dLitSoilC_10y, and
                   dLitSoilC_10y = dLitterIn_10y - dRh_10y (- fire): dead carbon still to be respired after Y+10
Relative metrics (denominators from the control before the pulse):
  nep_pre, npp_pre control NEP and NPP, mean of Y-3..Y-1 (years after the last clearcut only); ecoc_pre: Total-slow
                   in Y-1
  dNEP_pctNEP      dNEP as % of nep_pre, NaN when nep_pre < NEP_MIN (recently clearcut or near-equilibrium stands,
                   where the ratio is unstable); dUptake_10y_yrNEP: dUptake_10y / nep_pre, the 10-year loss in years
                   of the pre-pulse sink, same rule
  dUptake_10y_pctEcoC, dEcoC_10y_pctEcoC   as % of ecoc_pre; dVegC_5y_pctVegC: as % of vegc_pre
Caveats: the model has no physiological legacy (docs/LOG.md), so recovery shows mainly through mortality and
biomass; for very young stands use the absolute metrics (percentages are unstable when NPP is small). Outputs carry
3 decimals, so the placebo's dNPP is 0 or -0.001; its continuous metrics (dNEP, dUptake_10y) give the noise floor.

Management vs stand age (hot pulses): the stop stands are older and larger when the pulse comes, so the raw
difference between managements mixes management with the age/biomass effect. Four views:
  same-age table   (printed) metrics by age bin x management
  matched_age.csv  each managed run against the stop run with the same background (band, heat parameter, GCM, SSP,
                   pulse, pulse year) and the same age within AGE_TOL years: mean paired differences, next to the
                   unmatched difference
  mediation.csv    OLS on the common age support (age bins with >= 2 managements): total effect of management with
                   background controls only, direct effect with age bins (and log pre-pulse biomass) held fixed,
                   indirect = total - direct and the proportion mediated, i.e. the share of the management difference
                   that runs through the age structure; raw_diff_all_runs is the unadjusted difference over all runs
  pathways.csv     the 10-year budget by management and by age bin, drawn in fig_pulse_pathway.png

Output (pulse/): results.csv; phase_means.csv; landscape.csv (landscape BAU = W*bau26 + (1-W)*stop,
EXT = W*ext36 + (1-W)*stop, from phase means; dUptake_10y_pctEcoC from the weighted sums); age_vs_mgmt.csv (hot
pulses only); matched_age.csv; mediation.csv; pathways.csv; background_jas.csv (cache of the background pulse-season
temperatures); fig_pulse_temp.png (pulse effect vs background temperature, colour = stand age); fig_pulse_age.png
(pulse effect vs stand age, adjusted for background climate and heat parameter); fig_pulse_pathway.png (budget).
Usage: python analyze_pulse.py            metrics from the run directories, then tables and figures
       python analyze_pulse.py --replot   tables and figures only, from an existing pulse/results.csv
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

BASE = Path.home() / "lpjg_lishui"
RUNS = BASE / "runs_main"
OUT = BASE / "pulse"
W = 0.38                                   # commercial share of fir area (docs/sources.md); sensitivity 0.3-0.6
KEYS = ["elev_m", "mgmt", "pstemp_high", "gcm", "ssp", "phase"]
AGE_BINS = [0, 5, 10, 15, 20, 25, 30, 40, 60, 200]
AGE_LABELS = [f"{a}-{b - 1}" for a, b in zip(AGE_BINS[:-1], AGE_BINS[1:])]
LAND = ["dNPP", "dNEP", "dVegC_5y", "dEcoC_10y", "dRemoval_10y", "dUptake_10y", "dRemoval_2100", "ecoc_pre"]
TEST_METRICS = ("dNPP", "dUptake_10y", "dEcoC_10y")
NEW_COLS = ("dRemoval_10y", "dUptake_10y", "harv_events_10y", "yrs_to_cc")
DECOMP_COLS = ("dGPP", "dRh", "dNPP_post10", "dLitterIn_10y", "dUptake_10y_pctEcoC", "nep_pre")
CELL = (119.75, 28.25)                     # lon, lat of the simulated cell
MGMT_STYLE = {"bau26": ("tab:blue", "o", -0.7), "ext36": ("tab:orange", "s", 0.0), "stop": ("tab:green", "^", 0.7)}
NEP_MIN = 0.05                             # kg C m-2 yr-1: ratios to a smaller pre-pulse NEP are not reported
AGE_TOL = 2                                # years: age tolerance of the matched-age pairs
BG = ["elev_m", "pstemp_high", "gcm", "ssp", "pulse", "pulse_year"]   # the design without management and phase
REL = ("dNPP_pct", "dUptake_10y", "dUptake_10y_pctEcoC", "dEcoC_10y_pctEcoC")   # metrics of the age analyses
MONTHS = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split()
_cache = {}


def table(run, fname):
    if (run, fname) not in _cache:
        p = RUNS / run / "output" / fname
        _cache[run, fname] = (pd.read_csv(p, sep=r"\s+").groupby("Year").mean(numeric_only=True)
                              if p.exists() else None)
    return _cache[run, fname]


def series(run, fname, name):
    """one column as a Series by year; None when the file or the column is missing"""
    t = table(run, fname)
    return t[name] if t is not None and name in t.columns else None


def delta(prun, crun, fname, name):
    a, b = series(prun, fname, name), series(crun, fname, name)
    return None if a is None or b is None else a - b


def summed(parts):
    """sum of the Series that exist (None = missing); None when none of them exists"""
    have = [s for s in parts if s is not None]
    return sum(have[1:], have[0]) if have else None


def month_delta(prun, crun, fname, months):
    """pulse minus control, summed over the given months, from a monthly output (Jan..Dec columns)"""
    tp, tc = table(prun, fname), table(crun, fname)
    names = [MONTHS[m - 1] for m in months]
    if tp is None or tc is None or not set(names) <= set(tp.columns):
        return None
    return tp[names].sum(axis=1) - tc[names].sum(axis=1)


def at(s, y):
    return float(s[y]) if s is not None and y in s.index else np.nan


def wsum(s, a, b):
    return float(s.loc[a:b].sum()) if s is not None else np.nan


def pct(d, base):
    return float(d / base * 100) if np.isfinite(d) and np.isfinite(base) and base > 0 else np.nan


def finished(run):
    log = RUNS / run / "run.log"
    return log.exists() and "Finished" in log.read_text(errors="ignore")


def removals(run):
    cf, cp = table(run, "cflux.out"), table(run, "cpool.out")
    return cf["Harvest"] + cp["HarvSlowC"].diff().fillna(0) + cf["Slow_h"]


def clearcuts(run):
    dens, harv = table(run, "dens.out")["TeNE"], table(run, "cflux.out")["Harvest"]
    return [int(y) for y in dens.index if dens[y] < 1e-6 and harv.get(y, 0) > 0]


def recovery(vp, vc, y0, y_end):
    below = [y for y in range(y0, y_end) if y in vc.index and vp[y] < 0.99 * vc[y]]
    if not below:
        return 0
    return max(below) + 1 - y0 if max(below) + 1 < y_end else np.nan


def metrics(prun, crun, Y, months):
    npp_p, npp_c = table(prun, "anpp.out")["Total"], table(crun, "anpp.out")["Total"]
    nee_p, nee_c = table(prun, "cflux.out")["NEE-LU"], table(crun, "cflux.out")["NEE-LU"]
    veg_p, veg_c = table(prun, "cpool.out")["VegC"], table(crun, "cpool.out")["VegC"]
    eco_p, eco_c = table(prun, "cpool.out")["Total-slow"], table(crun, "cpool.out")["Total-slow"]
    mort_p, mort_c = table(prun, "tree_mort_sts.out")["Fir_CC"], table(crun, "tree_mort_sts.out")["Fir_CC"]
    harv_c = table(crun, "cflux.out")["Harvest"]
    d_rem = removals(prun) - removals(crun)
    cc = clearcuts(crun)
    past, future = [y for y in cc if y <= Y], [y for y in cc if y > Y]
    post = list(range(Y + 1, Y + 4))
    # flux decomposition: NEP = NPP - Rh - other, NPP = GPP - Ra; the vegetation budget gives the litter input
    d_npp, d_nep, d_mort = npp_p - npp_c, -(nee_p - nee_c), mort_p - mort_c
    d_gpp, gpp_c = delta(prun, crun, "agpp.out", "Total"), series(crun, "agpp.out", "Total")
    d_ra = None if d_gpp is None else d_gpp - d_npp
    d_rh, rh_c = delta(prun, crun, "cflux.out", "Soil"), series(crun, "cflux.out", "Soil")
    d_other = summed([delta(prun, crun, "cflux.out", c) for c in ("Fire", "Repr", "Est", "Seed")])
    d_vegout = summed([delta(prun, crun, "cflux.out", c) for c in ("Repr", "Est", "Seed")])  # NPP not kept in VegC
    d_ls = summed([delta(prun, crun, "cpool.out", c) for c in ("LitterC", "SoilC")])
    pm = {f: month_delta(prun, crun, f, months) for f in ("mgpp.out", "mra.out", "mrh.out")}
    pre = [y for y in range(Y - 3, Y) if (not past or y > max(past)) and y in npp_c.index]
    nep_pre = float(-nee_c[pre].mean()) if pre else np.nan
    ecoc_pre, vegc_pre = eco_c[Y - 1], veg_c[Y - 1]
    d_nep_Y, d_uptake = float(d_nep[Y]), float(d_nep.loc[Y:Y + 10].sum())
    d_eco, d_veg5 = eco_p[Y + 10] - eco_c[Y + 10], veg_p[Y + 5] - veg_c[Y + 5]
    d_veg10, d_rem10 = veg_p[Y + 10] - veg_c[Y + 10], float(d_rem.loc[Y:Y + 10].sum())
    return dict(
        pre_maxdiff=float((eco_p - eco_c).loc[:Y - 1].abs().max()),
        age=Y - max(past) if past else np.nan,
        vegc_pre=vegc_pre,
        dNPP=npp_p[Y] - npp_c[Y],
        dNPP_pct=(npp_p[Y] / npp_c[Y] - 1) * 100 if npp_c[Y] > 0 else np.nan,
        dNEP=d_nep_Y,
        dNPP_post_pct=(npp_p[post].mean() / npp_c[post].mean() - 1) * 100,
        dVegC_5y=d_veg5,
        dEcoC_10y=d_eco,
        dMort_3y=float(d_mort.loc[Y:Y + 2].sum()),
        rec_VegC=recovery(veg_p, veg_c, Y, future[0] if future else 2101),
        dRemoval_2100=float(d_rem.loc[Y:].sum()),
        # new columns are appended, so the original column order of results.csv is unchanged
        dRemoval_10y=d_rem10,
        dUptake_10y=d_uptake,
        harv_events_10y=int((harv_c.loc[Y:Y + 10] > 1e-6).sum()),
        yrs_to_cc=future[0] - Y if future else np.nan,
        # flux decomposition, pulse year
        dGPP=at(d_gpp, Y), dGPP_pct=pct(at(d_gpp, Y), at(gpp_c, Y)),
        dRa=at(d_ra, Y), dRa_pct=pct(at(d_ra, Y), at(gpp_c, Y) - npp_c[Y]),
        dRh=at(d_rh, Y), dRh_pct=pct(at(d_rh, Y), at(rh_c, Y)),
        dOther=at(d_other, Y),
        dGPP_pm=at(pm["mgpp.out"], Y), dRa_pm=at(pm["mra.out"], Y), dRh_pm=at(pm["mrh.out"], Y),
        # flux decomposition, pulse year .. Y+10 and the legacy years Y+1 .. Y+10
        dGPP_10y=wsum(d_gpp, Y, Y + 10), dRa_10y=wsum(d_ra, Y, Y + 10), dNPP_10y=wsum(d_npp, Y, Y + 10),
        dRh_10y=wsum(d_rh, Y, Y + 10), dOther_10y=wsum(d_other, Y, Y + 10),
        dNPP_post10=wsum(d_npp, Y + 1, Y + 10), dRh_post10=wsum(d_rh, Y + 1, Y + 10),
        # mortality, litter and the pools
        dMort_10y=wsum(d_mort, Y, Y + 10),
        dLitterIn_10y=wsum(d_npp, Y, Y + 10) - wsum(d_vegout, Y, Y + 10) - d_veg10 - d_rem10 if d_vegout is not None
        else wsum(d_npp, Y, Y + 10) - d_veg10 - d_rem10,
        dVegC_10y=d_veg10,
        dLitSoilC_10y=at(d_ls, Y + 10),
        # relative metrics
        nep_pre=nep_pre,
        npp_pre=float(npp_c[pre].mean()) if pre else np.nan,
        ecoc_pre=ecoc_pre,
        dNEP_pctNEP=d_nep_Y / nep_pre * 100 if nep_pre >= NEP_MIN else np.nan,
        dUptake_10y_yrNEP=d_uptake / nep_pre if nep_pre >= NEP_MIN else np.nan,
        dUptake_10y_pctEcoC=pct(d_uptake, ecoc_pre),
        dEcoC_10y_pctEcoC=pct(d_eco, ecoc_pre),
        dVegC_5y_pctVegC=pct(d_veg5, vegc_pre),
    )


def ols(y, X):
    """OLS with an intercept: coefficients, their standard errors, residual sum of squares and the rank"""
    X = np.column_stack([np.ones(len(y)), X]).astype(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    k, ssr = np.linalg.matrix_rank(X), float(r @ r)
    s2 = ssr / (len(y) - k) if len(y) > k else np.nan
    se = np.sqrt(np.clip(np.diag(np.linalg.pinv(X.T @ X)) * s2, 0, None))
    return beta, se, ssr, k


def ols_r2(y, X):
    beta, se, ssr, k = ols(y, X)
    return 1 - ssr / float(((y - y.mean()) ** 2).sum()), ssr, k


def hot_only(res):
    return res[res.pulse.astype(str).str.startswith("hot")]


def with_age_bins(res):
    """hot pulses with a known age and their age bin (ordered categories); second value: the rows in age bins shared
    by >= 2 managements, the common support where a same-age comparison is possible"""
    d = hot_only(res).dropna(subset=["age"]).copy()
    d["age_bin"] = pd.cut(d.age, AGE_BINS, right=False, labels=AGE_LABELS)
    return d, d[d.groupby("age_bin", observed=True).mgmt.transform("nunique") >= 2]


def age_vs_mgmt(res):
    """Does management explain the pulse effect beyond stand age? Hot pulses only: the placebo has no effect to
    explain and would inflate R2 through the pulse term. Only age bins with >= 2 managements are used."""
    _, d = with_age_bins(res)
    rows = []
    if d.empty or d.mgmt.nunique() < 2:
        return pd.DataFrame(rows), d
    for m in TEST_METRICS:
        e = d.dropna(subset=[m])
        base = pd.get_dummies(e[["age_bin", "elev_m", "pstemp_high", "gcm", "ssp", "pulse"]].astype(str),
                              drop_first=True)
        full = pd.concat([base, pd.get_dummies(e.mgmt, prefix="mgmt", drop_first=True)], axis=1)
        y = e[m].to_numpy(float)
        r2a, ssa, ka = ols_r2(y, base.to_numpy())
        r2b, ssb, kb = ols_r2(y, full.to_numpy())
        dfe = len(y) - kb
        f = ((ssa - ssb) / max(kb - ka, 1)) / (ssb / dfe) if dfe > 0 and ssb > 0 else np.nan
        try:
            from scipy.stats import f as fdist
            p = float(fdist.sf(f, kb - ka, dfe))
        except ImportError:
            p = np.nan
        rows.append(dict(metric=m, pulses="hot only", n=len(y), age_bins=e.age_bin.nunique(), R2_age_and_controls=r2a,
                         R2_plus_mgmt=r2b, dR2_mgmt=r2b - r2a, F=f, p=p))
    return pd.DataFrame(rows), d


def mgmt_dummies(d):
    """one 0/1 column per management other than stop, named like the contrast it estimates"""
    return pd.DataFrame({f"{m}-stop": (d.mgmt == m).astype(float).to_numpy()
                         for m in sorted(d.mgmt.unique()) if m != "stop"}, index=d.index)


def mediation(res):
    """management -> age structure -> pulse effect? On the common age support, for every metric in REL: the total
    effect of management (OLS with the background and pulse factors, BG, as controls), the direct effect once age
    bins (and log pre-pulse biomass) are also in the model, and indirect = total - direct, the part of the
    management difference that runs through the age structure (difference method; equals the product method in
    linear models). 'stop' is the reference. raw_diff_all_runs: the unadjusted difference over all hot runs."""
    allh, d = with_age_bins(res)
    rows = []
    if d.empty or d.mgmt.nunique() < 2 or "stop" not in set(d.mgmt):
        return pd.DataFrame(rows)
    for m in REL:
        e = d.dropna(subset=[m, "vegc_pre"])
        if e.mgmt.nunique() < 2:
            continue
        base = pd.get_dummies(e[BG].astype(str), drop_first=True).to_numpy(float)
        A = pd.get_dummies(e.age_bin.astype(str), drop_first=True).to_numpy(float)
        V = np.log(e.vegc_pre.clip(lower=1e-3).to_numpy(float))[:, None]
        M = mgmt_dummies(e)
        y = e[m].to_numpy(float)
        fit = {}
        models = (("total", [base]), ("age", [base, A]), ("biomass", [base, V]), ("age_biomass", [base, A, V]))
        for name, parts in models:
            beta, se, *_ = ols(y, np.column_stack(parts + [M.to_numpy()]))
            fit[name] = beta[-M.shape[1]:], se[-M.shape[1]:]
        raw = allh.groupby("mgmt")[m].mean()
        for j, c in enumerate(M.columns):
            tot, d_age, d_bio, d_ab = (fit[k][0][j] for k in ("total", "age", "biomass", "age_biomass"))
            rows.append(dict(metric=m, contrast=c, n=len(e), age_bins=e.age_bin.nunique(),
                             raw_diff_all_runs=raw[c.split("-")[0]] - raw["stop"],
                             total=tot, se_total=fit["total"][1][j],
                             direct_age=d_age, se_direct_age=fit["age"][1][j],
                             indirect_age=tot - d_age, prop_mediated_age=(tot - d_age) / tot if tot else np.nan,
                             direct_biomass=d_bio, prop_mediated_biomass=(tot - d_bio) / tot if tot else np.nan,
                             direct_age_biomass=d_ab, prop_mediated_age_biomass=(tot - d_ab) / tot if tot else np.nan))
    return pd.DataFrame(rows)


def matched_age(res, tol=AGE_TOL):
    """Same background (BG) and the same stand age within tol years: every managed hot-pulse run against the stop
    run closest in age. Mean paired differences per metric, next to the unmatched difference of the means."""
    h = hot_only(res).dropna(subset=["age"])
    stop = h[h.mgmt == "stop"]
    rows = []
    for mg in sorted(set(h.mgmt) - {"stop"}):
        j = h[h.mgmt == mg].merge(stop, on=BG, suffixes=("", "_stop"))
        j = j.assign(dage=(j.age - j.age_stop).abs())
        j = j[j.dage <= tol].sort_values("dage").drop_duplicates("run_id")
        for m in REL:
            dd = (j[m] - j[m + "_stop"]).dropna()
            rows.append(dict(contrast=f"{mg}-stop", metric=m, n_pairs=len(dd),
                             age_managed=j.age.mean() if len(j) else np.nan,
                             age_stop=j.age_stop.mean() if len(j) else np.nan,
                             mean_diff=dd.mean() if len(dd) else np.nan, sd_diff=dd.std() if len(dd) > 1 else np.nan,
                             share_managed_lower=(dd < 0).mean() if len(dd) else np.nan,
                             unmatched_diff=h[h.mgmt == mg][m].mean() - stop[m].mean(),
                             n_managed_all=int(h[h.mgmt == mg][m].notna().sum()),
                             n_stop_all=int(stop[m].notna().sum())))
    return pd.DataFrame(rows)


BUDGET = ["dGPP", "dRa", "dNPP", "dRh", "dOther", "dNEP", "dGPP_10y", "dRa_10y", "dNPP_10y", "dRh_10y", "dOther_10y",
          "dNPP_post10", "dRh_post10", "dUptake_10y", "dMort_3y", "dMort_10y", "dLitterIn_10y", "dVegC_10y",
          "dLitSoilC_10y", "dEcoC_10y", "dRemoval_10y"]


def budget(d):
    """mean 10-year budget of a set of hot-pulse runs: the contributions to dUptake_10y by year (pulse year Y, the
    legacy years after it) and by flux, then the mortality-litter-Rh chain"""
    g = d[BUDGET].mean()
    legacy = g.dUptake_10y - g.dNEP
    share = lambda a, b: a / b if np.isfinite(b) and b != 0 else np.nan
    return pd.Series({
        "n": len(d),
        "Y: dGPP": g.dGPP, "Y: -dRa": -g.dRa, "Y: dNPP": g.dNPP, "Y: -dRh": -g.dRh, "Y: -dOther": -g.dOther,
        "Y: dNEP": g.dNEP,
        "post: dGPP": g.dGPP_10y - g.dGPP, "post: -dRa": -(g.dRa_10y - g.dRa), "post: dNPP": g.dNPP_post10,
        "post: -dRh": -g.dRh_post10, "post: -dOther": -(g.dOther_10y - g.dOther), "post: dNEP (legacy)": legacy,
        "dUptake_10y": g.dUptake_10y, "legacy / dUptake_10y": share(legacy, g.dUptake_10y),
        "legacy via NPP": share(g.dNPP_post10, legacy), "legacy via Rh": share(-g.dRh_post10, legacy),
        "dMort_3y": g.dMort_3y, "dMort_10y": g.dMort_10y, "dLitterIn_10y": g.dLitterIn_10y, "dRh_10y": g.dRh_10y,
        "dLitSoilC_10y (to come)": g.dLitSoilC_10y, "dVegC_10y": g.dVegC_10y, "dRemoval_10y": g.dRemoval_10y,
        "dEcoC_10y": g.dEcoC_10y})


def pathways(res):
    """budget by management and by age bin (all managements), hot pulses; columns = groups"""
    allh, _ = with_age_bins(res)
    h = hot_only(res)
    cols = {f"mgmt {mg}": budget(sub) for mg, sub in h.groupby("mgmt")}
    cols.update({f"age {b}": budget(sub) for b, sub in allh.groupby("age_bin", observed=True)})
    return pd.DataFrame(cols)


def landscape(res):
    ph = res.groupby(["elev_m", "pstemp_high", "gcm", "ssp", "pulse", "mgmt"])[LAND].mean().reset_index()
    wide = ph.pivot_table(index=["elev_m", "pstemp_high", "gcm", "ssp", "pulse"], columns="mgmt", values=LAND)
    out = []
    for name, managed in (("BAU", "bau26"), ("EXT", "ext36"), ("STOP", "stop")):
        part = pd.DataFrame({m: W * wide[m][managed] + (1 - W) * wide[m]["stop"] for m in LAND})
        part["dUptake_10y_pctEcoC"] = part.dUptake_10y / part.ecoc_pre * 100   # landscape loss / landscape stock
        out.append(part.assign(landscape=name).reset_index())
    return ph, pd.concat(out, ignore_index=True)


def months_of(spec):
    """'7-9' -> [7, 8, 9]"""
    a, b = (int(x) for x in str(spec).split("-"))
    return list(range(a, b + 1))


def tas_path(g, s, e):
    """Control temperature file of a gcm x ssp x band, as named by the runner."""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run_main import climate_paths
        return Path(climate_paths({"gcm": g, "ssp": s, "elev_m": e, "pulse": "none"})["file_temp"])
    except Exception:  # runner not importable: same naming rule as run_main.climate_paths
        return BASE / "data" / "climate3b" / "elev" / f"{g}_{s}_tas_e{e}_1901_2100.nc"


def background_temperature(res, design):
    """Mean temperature of the pulse months in the pulse year, in the control climate of each gcm x ssp x band: the
    background the pulse is added to. Sources in this order: cache pulse/background_jas.csv, export/clim_annual.csv
    (export_for_r.py, Jul-Sep only), the climate files (netCDF4). Returns None if none of them works."""
    spec = {r["pulse"]: r["pulse_months"] for r in design if r["pulse"] != "none"}
    months = sorted({tuple(months_of(v)) for v in spec.values()})
    if len(months) != 1:
        print(f"background temperature: pulses use different months {months}; skipped")
        return None
    months = list(months[0])
    need = res[["gcm", "ssp", "elev_m", "pulse_year"]].drop_duplicates().astype({"elev_m": int, "pulse_year": int})
    key = ["gcm", "ssp", "elev_m", "pulse_year"]
    cache = OUT / "background_jas.csv"
    tag = f"{min(months)}-{max(months)}"
    if cache.exists():
        c = pd.read_csv(cache)
        c = c[c.months == tag]
        got = need.merge(c, on=key, how="left")
        if got.tas_bg.notna().all():
            return got
    clim = BASE / "export" / "clim_annual.csv"
    got = None
    if clim.exists() and months == [7, 8, 9]:
        c = pd.read_csv(clim).rename(columns={"Year": "pulse_year", "tas_jas": "tas_bg"})
        got = need.merge(c[key + ["tas_bg"]], on=key, how="left")
        if got.tas_bg.isna().any():
            got = None
    if got is None:
        try:
            from netCDF4 import Dataset, num2date
        except ImportError:
            print("background temperature: no cache, no export/clim_annual.csv and no netCDF4; fig_pulse_temp skipped")
            return None
        vals = []
        for (g, s, e), grp in need.groupby(["gcm", "ssp", "elev_m"]):
            path = tas_path(g, s, e)
            with Dataset(path) as nc:
                t = nc.variables["time"]
                dates = num2date(t[:], t.units, getattr(t, "calendar", "standard"))
                lat, lon = np.asarray(nc.variables["lat"][:]), np.asarray(nc.variables["lon"][:])
                i, j = int(np.argmin(abs(lat - CELL[1]))), int(np.argmin(abs(lon - CELL[0])))
                v = nc.variables["tas"]
                x = np.ma.filled(np.ma.asarray(v[:, i, j]).astype(float), np.nan)
                units = str(getattr(v, "units", ""))
            if units.strip().lower() in ("k", "kelvin") or np.nanmean(x) > 100:
                x = x - 273.15
            ym = np.array([(d.year, d.month) for d in dates])
            for y in grp.pulse_year:
                sel = (ym[:, 0] == y) & np.isin(ym[:, 1], months)
                vals.append(dict(gcm=g, ssp=s, elev_m=e, pulse_year=y, tas_bg=float(np.nanmean(x[sel]))))
        got = need.merge(pd.DataFrame(vals), on=key, how="left")
    got.assign(months=tag).round({"tas_bg": 3}).to_csv(cache, index=False)
    return got


def r2_dummies(d, m, cols):
    X = pd.get_dummies(d[cols].astype(str), drop_first=True).to_numpy(float)
    return ols_r2(d[m].to_numpy(float), X)[0]


def adjusted(h, m):
    """Partial residuals: metric m with the effects of background climate (gcm x ssp x band x pulse year) and
    pstemp_high removed. One OLS with these factors and the stand (management x phase x pulse) as dummies; the
    removed part is centred, so the overall mean is kept. What is left is the stand (age, management) signal."""
    d = h.dropna(subset=[m]).copy()
    bg = d.gcm + "|" + d.ssp + "|" + d.elev_m.astype(str) + "|" + d.pulse_year.astype(str)
    stand = d.mgmt + "|" + d.phase.astype(str) + "|" + d.pulse.astype(str)
    Xb = pd.get_dummies(bg, drop_first=True).to_numpy(float)
    Xp = pd.get_dummies(d.pstemp_high.astype(str), drop_first=True).to_numpy(float)
    Xs = pd.get_dummies(stand, drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(len(d)), Xb, Xp, Xs])
    beta, *_ = np.linalg.lstsq(X, d[m].to_numpy(float), rcond=None)
    nb, nt = Xb.shape[1], Xp.shape[1]
    removed = Xb @ beta[1:1 + nb] + Xp @ beta[1 + nb:1 + nb + nt]
    d["adj"] = d[m].to_numpy(float) - removed + removed.mean()
    return d


def placebo_band(res, m="dUptake_10y"):
    pl = res[res.pulse.astype(str).str.startswith("plc")][m].dropna()
    return float(np.percentile(pl.abs(), 95)) if len(pl) else None


ROWS = (("dNPP_pct", "NPP change in pulse year (%)", "{:+.2f} % per K"),
        ("dUptake_10y", "uptake change, pulse year + 10 yr\n(kg C m-2)", "{:+.4f} kg C m-2 per K"),
        ("dUptake_10y_pctEcoC", "uptake change, pulse year + 10 yr\n(% of pre-pulse ecosystem C)", "{:+.3f} % per K"))


def legend(fig, band, extra=()):
    handles = [Line2D([], [], color="0.4", marker=mk, ls="", label=mg) for mg, (_, mk, _) in MGMT_STYLE.items()]
    if any(band.values()):
        handles.append(Patch(color="0.8", label="placebo noise, 95% of |x| (value in the panel)"))
    fig.legend(handles=handles + list(extra), loc="lower center", ncol=len(handles) + len(extra), fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.02))


def shade(ax, b):
    if b:
        ax.axhspan(-b, b, color="0.8", zorder=0)
        ax.text(0.02, 0.03, f"placebo {b:.4g}", transform=ax.transAxes, fontsize=7, color="0.4", va="bottom")


def fig_temperature(h, band, path):
    pulses = sorted(h.pulse.unique())
    fig, axes = plt.subplots(len(ROWS), len(pulses), figsize=(5.2 * len(pulses) + 1.2, 3.6 * len(ROWS)),
                             squeeze=False, sharex=True)
    norm, cmap = Normalize(vmin=0, vmax=max(80, float(np.nanmax(h.age)))), plt.get_cmap("viridis")
    for j, pu in enumerate(pulses):
        sub = h[h.pulse == pu]
        for i, (m, lab, fmt) in enumerate(ROWS):
            ax = axes[i, j]
            if i >= 1:
                shade(ax, band.get(m))
            for mg, (_, mk, _) in MGMT_STYLE.items():
                s = sub[sub.mgmt == mg].dropna(subset=[m])
                ax.scatter(s.tas_bg, s[m], c=s.age, cmap=cmap, norm=norm, marker=mk, s=18, alpha=0.8, linewidths=0)
            s = sub.dropna(subset=[m])
            if len(s) > 2 and s.tas_bg.nunique() > 1:
                b, a = np.polyfit(s.tas_bg, s[m], 1)
                xx = np.array([s.tas_bg.min(), s.tas_bg.max()])
                ax.plot(xx, a + b * xx, color="k", lw=1)
                ax.text(0.98, 0.97, "slope " + fmt.format(b) + f", r = {s.tas_bg.corr(s[m]):.2f}",
                        transform=ax.transAxes, ha="right", va="top", fontsize=8)
            ax.axhline(0, color="grey", lw=0.5)
            if j == 0:
                ax.set_ylabel(lab)
            ax.set_title(pu if i == 0 else "")
        axes[-1, j].set_xlabel("Jul-Sep mean temperature of the control climate\nin the pulse year (deg C)")
    legend(fig, band)
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes, shrink=0.85, label="stand age in the pulse year (yr)")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_age(h, band, path):
    pulses = sorted(h.pulse.unique())
    fig, axes = plt.subplots(len(ROWS), len(pulses), figsize=(5.2 * len(pulses), 3.6 * len(ROWS)), squeeze=False,
                             sharex=True)
    for j, pu in enumerate(pulses):
        for i, (m, lab, _) in enumerate(ROWS):
            ax = axes[i, j]
            if i >= 1:
                shade(ax, band.get(m))
            d = adjusted(h, m)
            d = d[d.pulse == pu]
            for mg, (col, mk, off) in MGMT_STYLE.items():
                s = d[d.mgmt == mg]
                ax.scatter(s.age + off, s.adj, color=col, marker=mk, s=14, alpha=0.45, linewidths=0)
                mean = s.groupby("age").agg(adj=("adj", "mean"), soon=("yrs_to_cc", lambda x: (x <= 3).all()))
                edge = np.where((i >= 1) & mean.soon.to_numpy(), "red", "k")
                ax.scatter(mean.index + off, mean.adj, color=col, marker=mk, s=60, edgecolors=edge,
                           linewidths=np.where(edge == "red", 2.0, 0.8), zorder=3)
            ax.axhline(0, color="grey", lw=0.5)
            if j == 0:
                ax.set_ylabel(lab + "\nadjusted for background climate")
            ax.set_title(pu if i == 0 else "")
        axes[-1, j].set_xlabel("stand age in the pulse year (yr)\n(managements offset by 0.7 yr for visibility)")
    soon = Line2D([], [], color="w", marker="o", markeredgecolor="red", markeredgewidth=2, ls="",
                  label="mean of a stand clearcut <= 3 yr after the pulse")
    legend(fig, band, [soon])
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


COMPONENTS = (("Y: dGPP", "GPP, pulse year", "#1b9e77"), ("Y: -dRa", "-Ra, pulse year", "#a6d854"),
              ("Y: dNPP", "NPP, pulse year", "#1b9e77"), ("Y: -dRh", "-Rh, pulse year", "#d95f02"),
              ("post: dNPP", "NPP, years 1-10 after", "#7570b3"), ("post: -dRh", "-Rh, years 1-10 after", "#e7298a"),
              ("other", "other (fire, repr., est.)", "0.6"))


def fig_pathway(pw, path):
    """stacked bars of the mean contributions to dUptake_10y (negative = loss) by management and by age bin, hot
    pulses; the diamond is dUptake_10y itself, the circle the pulse-year dNEP"""
    pw = pw.copy()
    pw.loc["other"] = pw.loc["Y: -dOther"].fillna(0) + pw.loc["post: -dOther"].fillna(0)
    has_gpp = np.isfinite(pw.loc["Y: dGPP"].to_numpy(float)).any()
    skip = {"Y: dNPP"} if has_gpp else {"Y: dGPP", "Y: -dRa"}     # GPP and -Ra when agpp.out exists, else NPP
    comps = [c for c in COMPONENTS if c[0] not in skip]
    groups = [("by management", [c for c in pw.columns if c.startswith("mgmt ")]),
              ("by stand age in the pulse year", [c for c in pw.columns if c.startswith("age ")])]
    groups = [g for g in groups if g[1]]
    fig, axes = plt.subplots(1, len(groups), figsize=(2.2 + 1.1 * sum(len(g[1]) for g in groups), 4.8), squeeze=False,
                             sharey=True, gridspec_kw=dict(width_ratios=[max(len(g[1]), 1) for g in groups]))
    for ax, (title, cols) in zip(axes[0], groups):
        x = np.arange(len(cols))
        pos, neg = np.zeros(len(cols)), np.zeros(len(cols))
        for key, lab, col in comps:
            v = pw.loc[key, cols].to_numpy(float)
            v = np.where(np.isfinite(v), v, 0)
            up = v >= 0
            ax.bar(x, np.where(up, v, 0), bottom=pos, color=col, width=0.7)
            ax.bar(x, np.where(up, 0, v), bottom=neg, color=col, width=0.7)
            pos, neg = pos + np.where(up, v, 0), neg + np.where(up, 0, v)
        ax.scatter(x, pw.loc["dUptake_10y", cols].to_numpy(float), color="k", marker="D", s=28, zorder=4)
        ax.scatter(x, pw.loc["Y: dNEP", cols].to_numpy(float), facecolor="w", edgecolor="k", marker="o", s=28, zorder=4)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{c.split(' ', 1)[1]}\n(n={int(pw.loc['n', c])})" for c in cols], fontsize=8)
        ax.set_title(title, fontsize=10)
    axes[0, 0].set_ylabel("contribution to the 10-yr uptake change\n(kg C m-2; negative = loss)")
    handles = [Patch(color=col, label=lab) for _, lab, col in comps]
    handles += [Line2D([], [], color="k", marker="D", ls="", label="dUptake_10y (sum of the bars)"),
                Line2D([], [], markerfacecolor="w", markeredgecolor="k", marker="o", ls="", label="dNEP, pulse year")]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 4), fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compute(design):
    controls = {tuple(r[k] for k in KEYS): r["run_id"] for r in design if r["pulse"] == "none"}
    rows, missing = [], []
    for r in design:
        if r["pulse"] == "none":
            continue
        crun = controls[tuple(r[k] for k in KEYS)]
        if not (finished(r["run_id"]) and finished(crun)):
            missing.append(r["run_id"])
            continue
        rows.append({**{k: r[k] for k in ["run_id", *KEYS, "pulse", "pulse_year"]}, "control": crun,
                     **metrics(r["run_id"], crun, int(r["pulse_year"]), months_of(r.get("pulse_months", "7-9")))})
    if not rows:
        raise SystemExit(f"no finished pulse/control pairs ({len(missing)} missing)")
    res = pd.DataFrame(rows)
    for k in ("elev_m", "pstemp_high", "phase", "pulse_year"):
        res[k] = res[k].astype(int)
    res.to_csv(OUT / "results.csv", index=False)
    print(f"{len(res)} pulse runs analysed, {len(missing)} missing (unfinished pulse or control run)")
    return res


def same_age_table(allh, m):
    """means by age bin x management with n in brackets"""
    piv = allh.pivot_table(index="age_bin", columns="mgmt", values=m, aggfunc="mean", observed=True)
    cnt = allh.pivot_table(index="age_bin", columns="mgmt", values=m, aggfunc="count", observed=True)
    out = pd.DataFrame(index=piv.index)
    for c in piv.columns:
        out[c] = [f"{v:.3f} ({int(n)})" if n > 0 else "-" for v, n in zip(piv[c], cnt[c].fillna(0))]
    return out


def report(res, design):
    pd.set_option("display.width", 220)
    bad = res[res.pre_maxdiff > 1e-6]
    print("pre-pulse check: " + (f"WARNING {len(bad)} pairs differ before the pulse year: {bad.run_id.head(3).tolist()}"
                                 if len(bad) else "all pairs identical before the pulse year"))
    nan_pct = res[res.dNPP_pct.isna()]
    if len(nan_pct):
        print(f"note: {len(nan_pct)} runs have dNPP_pct = NaN (control NPP rounds to 0), e.g. "
              f"{nan_pct.run_id.head(2).tolist()}; they are left out of the % table and the NPP panels")

    for m, label in (("dNPP_pct", "NPP change in the pulse year (%)"),
                     ("dUptake_10y", "uptake change, pulse year + 10 years (kg C m-2) = dEcoC_10y + dRemoval_10y"),
                     ("dUptake_10y_pctEcoC", "uptake change, pulse year + 10 years, % of the pre-pulse ecosystem C")):
        print(f"\n== {label}; mean over GCMs and phases ==")
        print(res.pivot_table(index=["pstemp_high", "elev_m", "mgmt"], columns=["ssp", "pulse"], values=m)
              .round(3).to_string())

    h = hot_only(res)
    acc = h.groupby("mgmt").agg(dEcoC_10y=("dEcoC_10y", "mean"), dRemoval_10y=("dRemoval_10y", "mean"),
                                dUptake_10y=("dUptake_10y", "mean"),
                                share_with_harvest_in_window=("harv_events_10y", lambda x: (x > 0).mean()))
    print("\n== where the 10-year loss ends up, hot pulses (kg C m-2): ecosystem vs lower removals ==")
    print(acc.round(4).to_string())
    pos = h[h.dEcoC_10y > 0]
    if len(pos):
        print(f"   dEcoC_10y > 0 in {len(pos)} runs, {int((pos.harv_events_10y > 0).sum())} of them with harvest "
              f"in the window")
    t = h.assign(next_clearcut=pd.cut(h.yrs_to_cc.fillna(999), [0, 3, 6, 10, 1000],
                                      labels=["1-3 yr", "4-6 yr", "7-10 yr", ">10 yr or none"]))
    print("\n== uptake loss by years from the pulse to the control's next clearcut, hot pulses ==")
    print(t.groupby("next_clearcut", observed=True).agg(
        n=("run_id", "size"), dNEP_pulse_year=("dNEP", "mean"), dUptake_10y=("dUptake_10y", "mean"),
        share_kept=("dUptake_10y", lambda x: float(np.nanmedian(x / t.loc[x.index, "dNEP"])))).round(3).to_string())
    print("   share_kept: median dUptake_10y / pulse-year dNEP; the rest is made up within the 10 years")

    ph, land = landscape(res)
    land.to_csv(OUT / "landscape.csv", index=False)
    ph.to_csv(OUT / "phase_means.csv", index=False)
    for m, unit in (("dUptake_10y", "kg C m-2"), ("dUptake_10y_pctEcoC", "% of the landscape's pre-pulse ecosystem C")):
        print(f"\n== landscape (phase means, W = {W}): uptake change, pulse year + 10 years ({unit}), "
              "mean over GCMs ==")
        print(land.pivot_table(index=["pstemp_high", "elev_m", "landscape"], columns=["ssp", "pulse"], values=m)
              .round(3).to_string())

    # flux decomposition: identities, then the pathways of the loss
    r_y = (res.dNEP - (res.dNPP - res.dRh - res.dOther.fillna(0))).abs().max()
    r_w = (res.dUptake_10y - (res.dNPP_10y - res.dRh_10y - res.dOther_10y.fillna(0))).abs().max()
    r_d = (res.dLitSoilC_10y - (res.dLitterIn_10y - res.dRh_10y)).abs().max()
    print("\n== budget identities, max |residual| over all runs (rounding of the outputs is 0.001 per year) ==")
    print(f"   dNEP = dNPP - dRh - dOther: {r_y:.4f};  dUptake_10y = dNPP_10y - dRh_10y - dOther_10y: {r_w:.4f};  "
          f"dLitSoilC_10y = dLitterIn_10y - dRh_10y: {r_d:.4f} (fire, if any, is in the residual)")
    if res.dRh.isna().all():
        print("   no Soil column in cflux.out: the Rh columns are NaN")
    if res.dGPP.isna().all():
        print("   no agpp.out in the runs: GPP and Ra columns are NaN (add file_agpp to the output section)")
    if res.dGPP_pm.isna().all():
        print("   no mgpp/mra/mrh.out in the runs: the pulse-month columns are NaN")
    pw = pathways(res)
    pw.to_csv(OUT / "pathways.csv")
    print("\n== pathways of the 10-year uptake loss, hot pulses, kg C m-2 (means; 'Y' = pulse year, 'post' = Y+1..Y+10;"
          " dUptake_10y = Y: dNEP + post: dNEP) ==")
    print(pw.round(3).to_string())
    print("   legacy = the loss after the pulse year, split into lower NPP (lost capacity after mortality/biomass\n"
          "   loss) and higher Rh (decomposition of the extra dead material: dLitterIn_10y, with dMort for the\n"
          "   mortality part); the 'via' shares add up to 1 with the other term and lie outside 0..1 when the two\n"
          "   paths have opposite signs (e.g. growth release of the survivors); dLitSoilC_10y is the extra dead\n"
          "   carbon still in litter + soil at Y+10, i.e. respiration still to come; in the pulse year, dRh_pm\n"
          "   (pulse months) is mostly the direct temperature response of decomposition, the rest of the year\n"
          "   includes the fresh litter")
    if res.dGPP_pm.notna().any():
        pmt = h.groupby("mgmt")[["dGPP", "dGPP_pm", "dRa", "dRa_pm", "dRh", "dRh_pm", "dNEP"]].mean()
        pmt = pmt.assign(dGPP_rest=pmt.dGPP - pmt.dGPP_pm, dRa_rest=pmt.dRa - pmt.dRa_pm, dRh_rest=pmt.dRh - pmt.dRh_pm)
        print("\n== pulse year split: pulse months (_pm) vs the rest of the year (_rest), hot pulses ==")
        print(pmt.round(3).to_string())

    # management vs stand age: the confounding, the same-age table, matched pairs, mediation
    allh, cs = with_age_bins(res)
    conf = h.groupby("mgmt").agg(n=("run_id", "size"), age=("age", "mean"), age_min=("age", "min"),
                                 age_max=("age", "max"), vegc_pre=("vegc_pre", "mean"), ecoc_pre=("ecoc_pre", "mean"),
                                 nep_pre=("nep_pre", "mean"),
                                 dNPP_pct=("dNPP_pct", "mean"), dNEP_pctNEP=("dNEP_pctNEP", "mean"),
                                 n_NEPratio=("dNEP_pctNEP", "count"), dUptake_10y=("dUptake_10y", "mean"),
                                 dUptake_10y_yrNEP=("dUptake_10y_yrNEP", "mean"),
                                 dUptake_10y_pctEcoC=("dUptake_10y_pctEcoC", "mean"),
                                 dEcoC_10y_pctEcoC=("dEcoC_10y_pctEcoC", "mean"),
                                 dVegC_5y_pctVegC=("dVegC_5y_pctVegC", "mean"))
    print("\n== the stands the pulse hits, and absolute vs relative losses, by management (hot pulses) ==")
    print(conf.round(3).to_string())
    print(f"   ratios to NEP only where the pre-pulse NEP >= {NEP_MIN} kg C m-2 yr-1 (n_NEPratio runs)")
    print("\n== same-age comparison: means by age bin x management, hot pulses (n in brackets) ==")
    for m in REL:
        print(f"   {m}")
        print("   " + same_age_table(allh, m).to_string().replace("\n", "\n   "))
    ma = matched_age(res)
    ma.to_csv(OUT / "matched_age.csv", index=False)
    print(f"\n== matched pairs: managed vs stop run of the same background and the same age (+- {AGE_TOL} yr), "
          "hot pulses ==")
    if len(ma) and ma.n_pairs.max() > 0:
        print(ma.round(3).to_string(index=False))
        print("   mean_diff: managed minus stop over the pairs; unmatched_diff: difference of the means over all\n"
              "   runs; share_managed_lower: share of pairs in which the managed stand loses more (more negative\n"
              "   metric)")
    else:
        print("   no pairs: no managed and stop stands of the same age under the same background")

    av, d = age_vs_mgmt(res)
    av.to_csv(OUT / "age_vs_mgmt.csv", index=False)
    print("\n== does management matter beyond stand age? hot pulses only (OLS with age bins, elevation, heat parameter,"
          " GCM, SSP, pulse; age bins shared by >= 2 managements) ==")
    if len(av):
        print(av.round(4).to_string(index=False))
        print("   runs per age bin and management in the test:")
        print("   " + pd.crosstab(d.age_bin, d.mgmt).to_string().replace("\n", "\n   "))
    else:
        print("   not enough overlap between managements")
    md = mediation(res)
    md.to_csv(OUT / "mediation.csv", index=False)
    print("\n== mediation: management -> age structure -> pulse effect? common age support, hot pulses, "
          "stop = reference ==")
    if len(md):
        print(md.round(3).to_string(index=False))
        print("   total: management effect with the background and pulse factors as controls; direct_*: the same with\n"
              "   age bins / log pre-pulse biomass / both held fixed; prop_mediated = (total - direct) / total.\n"
              "   A direct effect near 0 with a sizeable total effect means the management difference runs through\n"
              "   the age structure; compare |total| with the placebo noise floor below before reading a proportion")
    else:
        print("   not enough overlap between managements")

    band = {m: placebo_band(res, m) for m, _, _ in ROWS}
    pl = res[res.pulse.astype(str).str.startswith("plc")]
    if len(pl):
        print("\n== placebo noise floor (95th percentile of |x|) vs median |x| of the hot pulses ==")
        for m in ("dNEP", "dVegC_5y", "dEcoC_10y", "dUptake_10y", "dUptake_10y_pctEcoC"):
            print(f"   {m:20s} placebo {np.nanpercentile(pl[m].abs(), 95):.4f} (max {pl[m].abs().max():.4f})"
                  f"   hot {h[m].abs().median():.4f}")

    if len(h):
        hh = h.dropna(subset=["dNPP_pct"])
        print("\n== share of the variance of the pulse-year NPP change (%), hot pulses ==")
        stand = hh.assign(stand=hh.mgmt + "|" + hh.phase.astype(str) + "|" + hh.pulse)
        stand["bg"] = stand.gcm + "|" + stand.ssp + "|" + stand.elev_m.astype(str) + "|" + stand.pulse_year.astype(str)
        for m in ("dNPP_pct", "dNPP"):
            r_st, r_bg = r2_dummies(stand, m, ["stand"]), r2_dummies(stand, m, ["bg"])
            r_all = r2_dummies(stand, m, ["stand", "bg", "pstemp_high"])
            print(f"   {m:9s} stand (management x phase) alone {r_st:.2f}; background (GCM x SSP x band x year) "
                  f"alone {r_bg:.2f}; both + pstemp_high {r_all:.2f}")

    figs = []
    bt = background_temperature(h, design) if len(h) else None
    if bt is not None:
        hb = h.merge(bt, on=["gcm", "ssp", "elev_m", "pulse_year"], how="left")
        cond = hb.groupby(["gcm", "ssp", "elev_m", "pulse_year"]).agg(dNPP_pct=("dNPP_pct", "mean"),
                                                                      tas_bg=("tas_bg", "first"))
        if len(cond) > 2:
            b, _ = np.polyfit(cond.tas_bg, cond.dNPP_pct, 1)
            print(f"\n== pulse-year NPP change vs background temperature of the pulse months ({len(cond)} GCM x SSP x "
                  f"band x year means): {b:+.2f} % per K, r = {cond.tas_bg.corr(cond.dNPP_pct):.2f} ==")
        fig_temperature(hb, band, OUT / "fig_pulse_temp.png")
        figs.append("fig_pulse_temp.png")
    if len(h):
        fig_age(h, band, OUT / "fig_pulse_age.png")
        figs.append("fig_pulse_age.png")
        fig_pathway(pw, OUT / "fig_pulse_pathway.png")
        figs.append("fig_pulse_pathway.png")
    print(f"\nWritten: {OUT}/results.csv, landscape.csv, phase_means.csv, age_vs_mgmt.csv, matched_age.csv, "
          "mediation.csv, pathways.csv" + (", background_jas.csv" if bt is not None else "")
          + (", " + ", ".join(figs) if figs else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="paired heat-pulse metrics, tables and figures")
    ap.add_argument("--replot", action="store_true", help="tables and figures only, from an existing pulse/results.csv")
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    design = list(csv.DictReader(open(BASE / "design.csv")))
    if a.replot:
        res = pd.read_csv(OUT / "results.csv")
        lacking = [c for c in NEW_COLS + DECOMP_COLS if c not in res.columns]
        if lacking:
            raise SystemExit(f"pulse/results.csv lacks {lacking}: run once without --replot")
    else:
        res = compute(design)
    report(res, design)
