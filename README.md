# SHyMA — Sounding-rocket Hypersonic Mission Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24568161.svg)](https://doi.org/10.5281/zenodo.24568161)

SHyMA is a framework of computationally light engineering models for the **pre-Phase A mission analysis of hypersonic (scramjet) flight experiments boosted by sounding rockets**. It chains these analyses:

- a stakeholder test window (Mach number, altitude, dynamic pressure, dwell time);
- a pseudo-five-degree-of-freedom ascent and point-mass flight over the WGS-84 Earth;
- an inverse design of the launch settings, deterministic or under uncertainty;
- a residual-dispersion Monte Carlo;
- hot-wall stagnation heating;
- telemetry coverage and the useful (continuous, link-closed) test time.

This repository accompanies the article *"Requirement-driven pre-Phase A mission analysis of sounding-rocket-boosted scramjet flight experiments"* (submitted to Acta Astronautica). It contains:

- the code;
- the public input data (VSB-30 flight table, HyShot 2 re-entry state, NASA NESC reference trajectories);
- the result files behind every number, table and figure of the article;
- test cases that reproduce the main results of the article (the long analyses are reproduced by `run_all.sh`);
- a single-file reviewer workbench that runs the same code in a web browser.

## Quick start for reviewers

**In the browser, no installation:** download `web/shyma_workbench.html` and open it. The first load downloads the Python runtime (Pyodide, about 10 MB), so an internet connection is needed. Then:

1. Pick a **paper test case** at the top of the page. The inputs of the corresponding tab are filled with the article's values, and the expected results are listed.
2. Press the tab's **Run** button and compare.
3. Change any input to explore: vehicle mass, launch settings, windows, dispersions, optimization goals.

**From the command line (Python 3.9 or later):**

```bash
pip install -r requirements.txt
python3 cases/run_case.py            # quick paper cases, about 3 min
python3 cases/run_case.py --all      # also the inverse design and the full 1000-run Monte Carlo sets (hours)
pytest -q tests/                     # the quick cases as tests, plus quick recomputations of Tables S3 (powered JSBSim case), S11 and S14
```

Each case file in `cases/` holds the inputs of one run and the values it must reproduce:

| Case | Article | Time |
|---|---|---|
| `nominal_A` … `nominal_D` | Table 4 (nominal concepts) | 5–15 s each |
| `montecarlo_X_first20` | Table 6: the first 20 runs of the paper set are reproduced to 1e-9 s | about 15–30 s each |
| `montecarlo_X_full` | Table 6, full 1000-run set (slow) | 20–60 min each |
| `optimizer_A` … `optimizer_D` | Table 5 (inverse design, slow) | 0.5–3 min each |
| `vsb30_validation` | Fig. 2, Table S19 | about 15 s |
| `nesc_checkcases` | Table S2, Fig. S1 | about 5 s |
| `formula_checks` | Equations of Secs. 3–4 and Supplementary S1 against independent references | under 1 s |

The same JSON files can be loaded in the workbench ("or a case file").

## Reproduce everything

```bash
./run_all.sh        # recomputes results/*.json, figures/*.pdf and results/manuscript_values.json
```

The result files of this release were computed in one run of `run_all.sh` on an eight-core Linux machine (`MC_PROCS=8`) in about 1 h 45 min, about 8 h on two cores (the Monte Carlo, selection and mass-ceiling scripts dominate; `MC_PROCS` sets the number of worker processes); the environment is recorded in `results/environment_run.txt` and the log in `results/chain_run.log`. The re-flight of Concept D under the ten-flight calibration in `results/robust_epistemic.json` was recomputed afterwards on two cores (`EPI_SETS=calibration_multiflight`); the runs reproduced on both machines agree bit for bit. The random seeds are fixed, so the result files are reproduced to round-off; the timing fields (`cpu_s`, `wall_s`) differ between runs.

## Where each result of the article comes from

| Article | Script | Result file |
|---|---|---|
| Fig. 1, Table 1 | diagram and literature review, no script | — |
| Table 2 (dispersions) | `shyma/dispersion.py` (`SIG`) | — |
| Table 3 (acceptance criteria), Table S1 (V&V coverage) | synthesis of the rows below | several, below |
| Sec. 5, Table S3 (code verification, JSBSim ballistic and powered) | `scripts/run_verification.py`, `scripts/run_weathercock_check.py`, `jsbsim_ref/run_jsbsim_ref.py`, `jsbsim_ref/run_jsbsim_thrust.py` | `results/verification.json`, `results/weathercock_check.json`, `results/jsbsim_verification.json`, `results/jsbsim_thrust.json` |
| Table S2, Fig. S1 (NESC) | `scripts/run_nesc_checkcases.py` | `results/nesc_checkcases.json` |
| Sec. 5, Fig. 2, Table S19 (VSB-30, residual tests, model-error bound) | `scripts/run_vsb30_validation.py`, `scripts/run_validation_stats.py` | `results/vsb30_validation.json`, `results/vsb30_zstats.json` |
| Sec. 5, Table S14 (launch-site atmosphere of the VSB-30 flights) | `scripts/run_site_bias.py` | `results/site_bias.json` |
| Sec. 5, Table S18 (ten-flight calibration) | `scripts/run_multiflight_calibration.py` | `results/multiflight_calibration.json` |
| Sec. 5, Fig. S2 (k_w calibration) | `scripts/run_wind_calibration.py` | `results/wind_calibration.json` |
| Sec. 5, Fig. S7, Table S4 (HyShot 2) | `scripts/run_hyshot_validation.py` | `results/hyshot_validation.json` |
| Sec. 5, Table S20 (validation domain) | `scripts/run_validation_domain.py` | `results/validation_domain.json` |
| Sec. 5, Tables S5–S6 (thrust shape, aerodynamic representation) | `scripts/run_thrust_shape.py`, `scripts/run_aero_sensitivity.py` | `results/thrust_shape.json`, `results/aero_sensitivity.json` |
| S2-H, S4-B (geometry-based stability derivatives and drag) | `scripts/run_geometric_aero.py`, `scripts/run_geometric_drag.py`, `scripts/run_stability_check.py` | `results/geometric_aero.json`, `results/geometric_drag.json`, `results/stability_check.json` |
| Sec. 5, Table S10 (nozzle exit and ambient-pressure thrust loss; held-out statistics of each variant) | `scripts/run_thrust_pressure.py` | `results/thrust_pressure.json` |
| Sec. 5.4, Table S13 (residual dispersion as an equivalent tip-off) | `scripts/run_tipoff_check.py` | `results/tipoff_check.json` |
| Sec. 4.2, S4-A (stage-mass sensitivity, recalibrated on V01) | `scripts/run_stage_mass_sensitivity.py` | `results/stage_mass_sensitivity.json` |
| Sec. 5 (effect of the HyShot 2 q bias on Concepts B and C) | `scripts/run_q_bias.py` | `results/q_bias.json` |
| Table 4, Figs. 3–4, Table S7, Figs. S3 and S8 (case study) | `scripts/run_casestudy.py` | `results/casestudy.json` |
| Sec. 4.3, S4-A (heating checks) | `scripts/run_thermal_model_check.py`, `scripts/run_aerothermal_check.py` | `results/thermal_model_check.json`, `results/aerothermal_check.json` |
| Sec. 6.3, S4-A, Table S11 (transient leading-edge temperature, heat loads) | `scripts/run_le_transient.py` | `results/le_transient.json` |
| Table 5 (inverse design, nominal mass ceilings) | `scripts/run_optimizer.py` | `results/optimizer.json` |
| Table 6, Figs. S4 and S9 (Monte Carlo; every set and concept uses the same draws, paired tests) | `scripts/run_montecarlo.py` | `results/montecarlo.json` |
| S4-C, Table S12 (Concept D: attitude error at S30 ignition) | `scripts/run_ignition_check.py` | `results/ignition_check.json` |
| Table S8 (rank correlations) | `scripts/run_mc_sensitivity.py` | `results/mc_sensitivity.json` |
| Sec. 6.4, Table S9, Fig. S5 (selection under uncertainty, useful test time) | `scripts/run_robust_selection.py` | `results/robust_selection.json` |
| Sec. 6.4, Table S17 (selected settings of D under the epistemic alternatives) | `scripts/run_robust_epistemic.py` | `results/robust_epistemic.json` |
| Sec. 6.4, Table S15 (local screen of Concept D, independent confirmation) | `scripts/run_robust_local.py` | `results/robust_local.json` |
| S3-E, S4-C (mass margins, fixed-β and fixed-geometry vehicles) | `scripts/run_mass_budget.py` | `results/mass_budget.json` |
| Table 5, Table S16, Fig. S6 (probabilistic mass ceilings, interpolated grid proportions with bootstrap band) | `scripts/run_mass_ceiling_mc.py` | `results/mass_ceiling_mc.json`, `results/mass_ceiling_mc_outcomes.json` |
| Table 7 (concept comparison) | synthesis of the results above | — |
| All figures | `scripts/make_figures.py` | `figures/` (generated) |

`results/manuscript_values.json` lists every computed number printed in the article and its Supplementary Material; its keys are the placeholder names used in the manuscript source. It is regenerated from the result files above by `scripts/fill_manuscript.py` (the last step of `run_all.sh`; without the manuscript templates, which are not part of this repository, the script writes only this file), so each printed value can be traced to its result file and script.

## Code

| Path | Content |
|---|---|
| `shyma/atmosphere.py` | US76 (0–86 km), upper-atmosphere density tables (> 100 km), optional NRLMSISE-00 (site and date) via `pymsis` |
| `shyma/earth.py` | WGS-84, spherical or round Earth; J2 gravity, frames, stations (CLA, CLBI, Esrange, Woomera), S-band link margin |
| `shyma/vehicle.py` | Multistage solid-motor model (VSB-30): constant or linear thrust law with ambient-pressure loss (exit diameter 0.45 m), test vehicle with ballistic coefficient (its excess drag over the generic nose is added to the attached stack), captive option |
| `shyma/trajectory.py`, `shyma/trajectory5.py`, `shyma/flight.py` | 3-DOF ECI propagation with wind (RK4, event handling); pseudo-5-DOF ascent; flight driver |
| `shyma/wind.py`, `shyma/dispersion.py` | Residual-dispersion model, wind weighting, Monte Carlo dispersions and statistics |
| `shyma/analysis.py` | Test-window metrics, hot-wall Sutton–Graves heating, coverage |
| `shyma/calibration.py` | Two-parameter calibration (Isp, S31 burn time) on one flight (`scripts/run_multiflight_calibration.py`: on all ten flights) |
| `shyma/geometry.py` | Parametric vehicle geometry (body profile, trapezoidal fin sets, stage configurations), mass properties, VSB-30 estimate; body profile from surface points (e.g. STL vertices) |
| `shyma/stability.py` | Static and dynamic pitch derivatives (C_N_alpha, centre of pressure, C_m_alpha, C_mq) from the geometry by component build-up (slender body, Diederich and supersonic linear theory, Pitts–Nielsen–Kaattari interference); damped-oscillation identification for the inverse check |
| `shyma/drag.py` | Zero-lift drag from the geometry: Van Driest II friction, shock-expansion and Taylor–Maccoll body wave drag, fin wave and leading-edge drag, base drag with jet reduction |
| `scripts/` | The analyses of the article (table above); `launch_optimizer.py` is the inverse design, with a three-state requirement report (met, violated, not evaluated); `run_casestudy.useful_window` gives the useful test time |
| `web/` | Reviewer workbench: `shyma_workbench.html` (built), `template.html`, `webdriver.py` (browser driver of the unmodified code), `build_html.py` |
| `cases/`, `tests/` | Paper test cases and their runner; pytest wrapper |

## Data and sources

| File | Content and source |
|---|---|
| `data/vsb30_flights_garcia2011.csv` | First ten VSB-30 flights: payload mass, launcher angles, operators' nominal and measured apogee, range and time above 100 km. Garcia et al., *J. Aerosp. Technol. Manag.* 3(3), 2011, doi:10.5028/jatm.2011.03032211. |
| `data/hyshot2_flight.json` | HyShot 2 re-entry state and test window. Hass, Smart and Paull, AIAA 2005-3354 (NASA NTRS 20050215040); Smart, Hass and Paull, *AIAA J.* 44(10), 2006, doi:10.2514/1.20661. |
| `validation_data/` | NASA Engineering and Safety Center atmospheric check-case reference trajectories, cases 5, 9, 10 and 17 (NASA/TM-2015-218675 and the NESC check-case data distribution), publicly released by NASA and redistributed unmodified. |
| `jsbsim_ref/aircraft/sphere/`, `jsbsim_ref/aircraft/thruster/` | Ballistic and powered point-mass definitions for the JSBSim comparisons. |

The model uses only public data. No data of the Brazilian 14-X flights and no proprietary motor data are included.

**Terms of the third-party data.** The VSB-30 and HyShot 2 values are numerical facts transcribed from the cited publications, which must be cited when they are reused. The NESC reference trajectories are a public NASA data release, redistributed unmodified with their citation. The workbench loads, at run time and from public CDNs, Pyodide 0.26.4 (MPL-2.0; cdn.jsdelivr.net) and Chart.js 4.4.1 (MIT; cdnjs.cloudflare.com); it therefore needs an internet connection and those two fixed URLs.

## Integrity

`HASHES.sha256` lists the SHA-256 of every file. Verify with `shasum -a 256 -c HASHES.sha256` (macOS) or `sha256sum -c HASHES.sha256` (Linux). The workbench shows the hashes of the code it embeds (About tab).

Results were produced with Python 3.13 on Linux (eight cores), and the paper test cases and the test suite were re-run on a second Linux machine; `requirements-lock.txt` lists the exact versions (NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, pymsis 0.13.0, JSBSim 1.3.3, pytest 9.1.1). The code runs on Python 3.9 or later by syntax, but only Python 3.13 was tested. The trajectory core uses only the standard library (`random.Random(seed)`), so the results are reproduced to round-off; agreement to 1e-9 can depend on the platform's math library.

## Archive

Each release of this repository is archived on Zenodo under the concept DOI [10.5281/zenodo.24568161](https://doi.org/10.5281/zenodo.24568161), which resolves to the latest version. Version 2.5.1 accompanies the Acta Astronautica submission; it differs from version 2.5 (DOI [10.5281/zenodo.24573924](https://doi.org/10.5281/zenodo.24573924)) only in the script that writes the numbers of the paper (conventional half-up rounding, elevations printed with their decimals, two added keys) and in figure labels, with the same results. Version 2.5: it adopts the stage masses of the IAE motor table (Palmerio, 2017; Kasemodel, 2010), sets the nominal elevation of Concept D to the optimum of the inverse design, and recomputes every result with the released code; the run was resumed once, from `scripts/run_robust_selection.py`, after the candidate grid of Concept D was centred on the new nominal setting (`results/chain_run.log`). Version 2.4, the earlier submission version, has the DOI [10.5281/zenodo.24568162](https://doi.org/10.5281/zenodo.24568162). The table, figure and section references of the paper test cases and of the workbench follow the submitted manuscript.

## Authors

- Pedro Kukulka de Albuquerque, Aeronautics Institute of Technology (ITA) and Institute for Advanced Studies (IEAv/DCTA), São José dos Campos, Brazil. ORCID [0009-0007-0275-9183](https://orcid.org/0009-0007-0275-9183)
- Lucas Galembeck, Aeronautics Institute of Technology (ITA) and Institute for Advanced Studies (IEAv/DCTA), São José dos Campos, Brazil. ORCID [0009-0002-1000-1069](https://orcid.org/0009-0002-1000-1069)

## License and citation

The code is released under the MIT License (`LICENSE`). If you use SHyMA, please cite the article; `CITATION.cff` has the details. Third-party data keep their original terms (see *Data and sources*).
