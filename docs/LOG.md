# Lishui LPJ-GUESS project log
Format: date - what - why - result/decision. Newest entries at the bottom.

## Setup and calibration (up to 2026-09-26)
- Model: LPJ-GUESS r13130 (forestry module), libraries from 4.1.1, one-line patch in cfinput.cpp.
- Historical climate: ISIMIP3a GSWP3-W5E5 obsclim, daily 0.5 deg, 1901-2019 (tas, pr, rsds), box 118.5-120.5E 27.0-29.5N.
- CO2: ISIMIP3a obsclim 1850-2021. N deposition: Lamarque (historic). Soil code 3 (soils_lpj.dat).
- Representative cell 119.75E 28.25N (cf gridlist index x=2 y=2).
- Calibration: respcoeff TeBE and TeNE 1.0 -> 0.3 (NPP 0.71, NPP/GPP 0.41).
- Drought: TeNE root_beta 0.982 -> 0.96 (range 0.95-0.97); default roots did not register the 2003 drought.
- Findings: heat (+2 C Jul-Sep: about -7.5 % NPP) >> drought (-30 % rain: about -2 %); density effect <= ~3 points;
  rotation controls carbon stock (16 yr 166, 26 yr 181, 41 yr 199, no harvest 294 t C/ha); harvest peaks at 26 yr.
- Open issues: model NPP high but biomass low vs field data; generic temperate heat parameters; no drought legacy.

## Decision: main experiment
- Elevation band (lapse-rate offset) x management from 2025 x heat parameter x GCM x SSP, 4 planting phases (design.csv).

## 2026-09-26 Literature and policy review
- Calibration targets added (docs/calib_targets.csv); model VegC ~1.4-1.9x and soil+litter ~1.5x above inventory means.
- Elevation bands changed to 200/500/800 m (fir mostly below 800 m); area weights still missing.
- Extended rotation changed from 41 to 36 yr (large-timber main felling age in LY/T 2908 implementation).
- BAU harvest intensity: observed ~1.5 %/yr vs 3.8 %/yr for 26-yr rotation -> commercial share w = 0.38 applied in analysis.

## 2026-09-26 Growth calibration and future-climate check
- TeNE respcoeff = 0.40 (RMSE 0.12 kg C m-2; middle 2.18 and near-mature 3.45 within target ranges; young 0.71, -19 %;
  NPP/GPP 0.43). Repeated run gave identical results. Soil+litter 13.1 kg C m-2 stays above field values (7.3-10.2):
  report changes, not absolute soil stocks. alphaa is not a PFT parameter in r13130.
- ISIMIP3b download complete (18 zips; 12 merged files 1901-2100, no gaps or missing values; CO2 1850-2100).
  Lishui cell 2071-2100 vs 1981-2010: +1.7 (GFDL 126), +3.7 (GFDL 585), +3.3 (UKESM 126), +7.7 C (UKESM 585);
  precipitation +13 to +25 %. UKESM 585 is outside the tested temperature range: treat as upper bound.

## 2026-09-26 Restart mechanism
- state_year counts simulation years from the first spin-up year (framework.cpp: date.year == state_year-1;
  default state_year = nyear_spinup). Calendar year Y -> state_year = nyear_spinup + (Y - 1901).
- First attempt with state_year 2000 saved nothing (0.state 8 bytes); restart failed "failed to find element to deserialize".

## 2026-09-26 Restart abandoned; N deposition patch
- Restarted BAU: same management events and stem numbers as the uninterrupted run, seed is serialized, but NPP differs
  from the first restarted year and is ~25 % higher in 2015-2019 -> some state is not restored in r13130.
  Decision: no restart; switch management in 2025 with management types (mt) inside uninterrupted 1901-2100 runs.
- cfinput.cpp line 634: Lamarque::RCP60 -> Lamarque::parse_timeseries(ndep_timeseries) (docs/patches/).
  Observed-climate base set to ndep_timeseries "rcp60" (identical to pre-patch); future: ssp126 -> rcp26, ssp585 -> rcp85.

## 2026-09-26 N deposition data and management switching
- ndep rcp26/rcp85 need GlobalNitrogenDepositionRCP26.bin / RCP85.bin next to the base file; the Zenodo package
  (LPJ-GUESS-CNP driver files) only has the base file (historic + RCP6.0). Decision: RCP6.0 N deposition for all future
  scenarios (common N forcing), "historic" as N sensitivity; request RCP files from the LPJ-GUESS team if needed.
- Management switching via mt blocks: management1 (26-yr) -> management2 in calendar year startman2 (rot_wait_for_cc 0),
  management2 with cutfirstyear 0 and cutfirstyear_nonsel 0. Test: scripts/make_mt_test.py (switch in 2000).

## 2026-09-26 Management switching validated; main runner
- mt switching test (switch in 2000): switch to an identical copy reproduces the reference exactly (max diff 0);
  ext36 skips the 2005 clearcut and cuts in 2015 (age 36); stop has no harvest after 2000.
  Restart runs had overestimated growth by 4-11 % (rs_ext 2007 5.00 vs 4.80; rs_stop 2019 7.69 vs 6.95).
- N deposition: historic vs rcp60 changes NPP 2006-2019 by -1.1 %. Main runs use rcp60.
- Main runner scripts/run_main.py: uninterrupted 1901-2100 runs from design.csv, management switch in 2025.
