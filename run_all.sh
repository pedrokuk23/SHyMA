#!/usr/bin/env bash
# Recomputes every result file of results/ and every figure of the article (figures/).
# About 1.5 h on an Apple M5 Pro and several hours on an older 2-core laptop; the Monte Carlo, the selection under uncertainty, the wind calibration and the
# thrust-shape studies dominate.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m shyma.calibration                   # Isp, t_b1 on VSB-30 flight V01 (WGS-84, 5-DOF ascent)
python3 scripts/run_verification.py            # Kepler order, step convergence, window dwell convergence
python3 jsbsim_ref/run_jsbsim_ref.py           # code-to-code vs JSBSim (needs: pip install jsbsim)
python3 jsbsim_ref/run_jsbsim_thrust.py         # code-to-code vs JSBSim with propulsion (thrust, mass flow, pressure loss)
python3 scripts/run_nesc_checkcases.py         # NASA NESC check-cases 5, 9, 10, 17
python3 scripts/run_weathercock_check.py       # 5-DOF attitude dynamics vs linear short-period solution
python3 scripts/run_aero_sensitivity.py        # stability coefficients and drag curve, frozen vs recalibrated propulsion
python3 scripts/run_stability_check.py         # geometry-based derivatives: exact limits and inverse dynamic check
python3 scripts/run_geometric_aero.py          # VSB-30 geometry-based derivatives and their mission effect
python3 scripts/run_geometric_drag.py          # VSB-30 geometry-based drag (with and without the derivatives)
python3 scripts/run_wind_calibration.py        # residual-dispersion amplitude k_w (N = 600, seed 7)
python3 scripts/run_vsb30_validation.py        # held-out flights, leave-one-out, statistics, normalised residuals
python3 scripts/run_validation_stats.py        # residual test including the calibration uncertainty
python3 scripts/run_hyshot_validation.py       # HyShot 2 re-entry reconstruction
python3 scripts/run_casestudy.py               # concepts A-D: nominal, sensitivities, trade grids
python3 scripts/run_thermal_model_check.py     # thermally vs calorically perfect air
python3 scripts/run_validation_domain.py       # flight envelope of the validation flights vs the concepts
python3 scripts/run_aerothermal_check.py       # Sutton-Graves vs Fay-Riddell at the window states
python3 scripts/run_thrust_shape.py            # dwell of Concepts A and D vs thrust-curve shape (recalibrated)
python3 scripts/run_optimizer.py               # inverse design of the launch settings, mass ceilings
python3 scripts/run_multiflight_calibration.py # ten-flight calibration with site atmospheres, leave-one-out
MC_N=1000 MC_N_2X=300 python3 scripts/run_montecarlo.py   # Monte Carlo (seed 20261006)
python3 scripts/run_mc_sensitivity.py          # rank correlations of the Monte Carlo dwell
python3 scripts/run_robust_selection.py        # settings selected under uncertainty (useful test time)
python3 scripts/run_robust_epistemic.py       # settings of Concept D under the epistemic alternatives
python3 scripts/run_mass_budget.py             # mass margins, fixed-beta and fixed-geometry vehicles
python3 scripts/run_thrust_pressure.py         # ambient-pressure thrust loss (recalibrated booster)
python3 scripts/run_q_bias.py                  # effect of the HyShot 2 dynamic-pressure bias on Concepts B and C
python3 scripts/run_robust_local.py            # local screen around the setting of Concept D selected under uncertainty
python3 scripts/run_mass_ceiling_mc.py         # probabilistic mass ceilings of the nominal settings
python3 scripts/run_le_transient.py           # transient lumped leading-edge temperature vs equilibrium bound
python3 scripts/run_ignition_check.py         # Concept D: attitude error at S30 ignition, weaker static stability
python3 scripts/run_tipoff_check.py           # residual dispersion as an equivalent tip-off instead of k_w
python3 scripts/run_site_bias.py              # launch-site atmosphere (NRLMSISE-00) and the CLA/Esrange split
python3 scripts/make_figures.py                # figures/*.pdf
python3 scripts/fill_manuscript.py             # results/manuscript_values.json: every number printed in the article
python3 web/build_html.py                      # single-file reviewer workbench
python3 cases/run_case.py                      # paper test cases
./make_hashes.sh
