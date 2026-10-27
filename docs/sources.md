# Sources and extracted values (retrieved 2026-09-26)

## Calibration targets (values in docs/calib_targets.csv)
- Diao J, Liu J, Zhu Z, Wei X, Li M (2022) Active forest management accelerates carbon storage in plantation
  forests in Lishui, southern China. Forest Ecosystems 9:100004. doi:10.1016/j.fecs.2022.100004
  Remote sensing 1989-2019: Chinese fir disturbed 4484 ha/yr; 14 % of area cut at least once.
  Historic Lishui fir rotation 26 yr; authors recommend ~50 yr and conifer-to-broadleaf conversion.
- Li Y et al. (2016) Carbon storage and its distribution of forest ecosystems in Zhejiang Province, China.
  Chinese Journal of Plant Ecology 40:354-363. doi:10.17521/cjpe.2015.0193
  Chinese fir volume-to-biomass: B = 0.4652 V + 19.141 (t/ha, V in m3/ha).

## Elevation of Chinese fir plantations (qualitative only)
- Longquan: fir plantations interspersed in the zone below 800 m; Huangshan pine dominates 800-1300 m.
- Jinyun: Huangshan pine above 700-800 m, Masson pine below.
- Decision: elevation bands 200 / 500 / 800 m. Area weights pending (inventory summary or map x DEM overlay).

## Policy constraints
- LY/T 2908-2017 (Fujian 2023 rules): fir main felling age >= 26 (medium), >= 31 (intensive large), >= 36 (general large timber).
- Zhejiang 2021 / 2025: main felling age set by operator; clearcut volume <= 50 % of main felling;
  provincial public-welfare forest: selective regeneration felling only; selective/thinning <= 40 % of volume.
- NFGA near-natural model (pine/fir): fir target diameter 40 cm at 40-45 yr; fir phased out by 65-70 yr.
- Lishui forest carbon action plan 2022-2025: forest carbon storage 80 Mt by 2025.

## Derived parameters
- Commercial (rotation-managed) share of fir area: w = (4484 / 304200) * 26 = 0.38
  (paper reports the fir area as "30.42 Mha", interpreted as 30.42e4 ha).
  Landscape BAU = w * bau26 + (1 - w) * stop; landscape EXT = w * ext36 + (1 - w) * stop; sensitivity w = 0.3-0.6.
- Model check vs targets (before calibration): VegC at age ~13 about 3.0 and at ~23 about 5.0 kg C m-2
  (targets ~2.1 and ~3.4); soil + litter ~15 kg C m-2 (target ~10). Calibrate before the main runs.
