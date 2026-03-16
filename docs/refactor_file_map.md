# Refactor File Map

This file records where the pre-refactor files moved during the `src/` cleanup.

## Core Source Files

| Old path | New path | Notes |
|---|---|---|
| `src/Reference_Trajectory_Calculation.m` | `examples/reference_trajectory_calculation.m` | Example/reference script |
| `src/Sat_template.m` | `src/+vleo/+dynamics/sat_dynamics_nonlinear.m` | Main nonlinear dynamics model |
| `src/Sat_template2_linear.m` | `src/+vleo/+dynamics/sat_dynamics_linearized.m` | Linearized dynamics model |
| `src/Sat_template_control.m` | `src/+vleo/+dynamics/sat_dynamics_controlled.m` | Controlled pointing dynamics |
| `src/build_render_gif_binary.sh` | `tools/build_render_gif_binary.sh` | GIF build helper |
| `src/computeAeroForces.m` | `src/+vleo/+aero/compute_aero_forces.m` | Renamed and packaged |
| `src/control_test3.m` | `scripts/volcano_tracking_demo.m` | Runnable scenario script |
| `src/controllers/control_torques.m` | `src/+vleo/+control/attitude_pd_controller.m` | Primary controller implementation |
| `src/controllers/myController.m` | `src/+vleo/+control/zero_torque_controller.m` | Null/placeholder controller |
| `src/keplerian_to_ijk_safe.m` | `src/+vleo/+analysis/keplerian_to_eci_safe.m` | Same job, clearer name |
| `src/render_gif` | `tools/bin/render_gif` | Built helper binary |
| `src/render_gif.py` | `tools/render_gif.py` | Python GIF renderer |
| `src/sat_template_gui.m` | `src/+vleo/+dynamics/sat_dynamics_gui.m` | GUI dynamics entry point |
| `src/state_to_observation.m` | `src/+vleo/+analysis/state_to_observation.m` | Observation helper |

## GUI Files

| Old path | New path | Notes |
|---|---|---|
| `src/gui/applyButtonStyle.m` | `src/+vleo/+gui/apply_button_style.m` | Renamed to snake_case |
| `src/gui/displayResults.m` | `src/+vleo/+gui/display_results.m` | Results plotting entry point |
| `src/gui/goBack.m` | `src/+vleo/+gui/go_back.m` | Renamed to snake_case |
| `src/gui/openLoadSimulation.m` | `src/+vleo/+gui/open_load_simulation.m` | Renamed and packaged |
| `src/gui/openNewSimulation.m` | `src/+vleo/+gui/open_new_simulation.m` | Renamed and packaged |
| `src/gui/openResults.m` | `src/+vleo/+gui/open_results_window.m` | Renamed for clarity |
| `src/gui/openSimulationGUI.m` | `src/+vleo/+gui/open_simulation_gui.m` | Main GUI function |
| `src/gui/openSimulationGUI.m` | `scripts/open_simulation_gui.m` | New runnable wrapper script |
| `src/gui/runSimulation.m` | `src/+vleo/+gui/run_simulation.m` | Renamed and packaged |
| `src/gui/functionalButtons.m` | `src/+vleo/+gui/+dialogs/` | Split into dedicated dialog files |
| `src/gui/AttitudeAnimatorGUI.m` | `src/+vleo/+viz/animate_results.m` | GUI animation logic consolidated here |
| `src/gui/AttitudeAnimator.m` | `examples/legacy/attitude_animator_demo.m` | Example/demo wrapper around new animation flow |
| `src/gui/AttitudeAnimator.m` | `src/+vleo/+viz/animate_results.m` | Core animation implementation |
| `src/gui/Copy_of_AttitudeAnimator.m` | `src/+vleo/+viz/animate_results.m` | Duplicate removed; canonical replacement |

## Setup and Path Bootstrapping

| Old path | New path | Notes |
|---|---|---|
| `src/ensure_project_setup.m` | `setup_project.m` | Duplicate removed; use root setup |
| `src/controllers/ensure_project_setup.m` | `setup_project.m` | Duplicate removed; use root setup |
| `src/gui/ensure_project_setup.m` | `setup_project.m` | Duplicate removed; use root setup |

## Data and Assets

| Old path | New path | Notes |
|---|---|---|
| `src/dynamics/6U CubeSat.obj` | `data/geometry/6U CubeSat.obj` | Geometry asset |
| `src/dynamics/6U CubeSat.STL` | `data/geometry/6U CubeSat.STL` | Geometry asset |
| `src/dynamics/6U CubeSat.mat` | `data/geometry/6U CubeSat.mat` | Geometry asset |
| `src/gui/pitches.csv` | `data/legacy/gui_csv/pitches.csv` | Legacy GUI data |
| `src/gui/rolls.csv` | `data/legacy/gui_csv/rolls.csv` | Legacy GUI data |
| `src/gui/yaws.csv` | `data/legacy/gui_csv/yaws.csv` | Legacy GUI data |

## Tests and Examples Renamed Outside `src/`

| Old path | New path | Notes |
|---|---|---|
| `examples/example_computeAeroForces.m` | `examples/example_compute_aero_forces.m` | Renamed to match packaged function |
| `tests/test_computeAeroForces.m` | `tests/test_compute_aero_forces.m` | Renamed to match packaged function |
| `tests/test_template.m` | `tests/templates/test_template.m` | Moved into templates folder |
| `tests/Sputnik_test.m` | `examples/legacy/sputnik_test.m` | Better fit as legacy example |
| `tests/test.m` | `tests/legacy/incomplete_test.m` | Archived incomplete test |
| `src/test_Sat_template.m` | `tests/test_sat_dynamics_nonlinear.m` | Replaced by a clearer smoke/regression test |

## Removed or Replaced Without a 1:1 Move

| Old path | Replacement | Notes |
|---|---|---|
| `src/Sputnik_test2.m` | none | Removed as stale exploratory script |
| `src/.DS_Store` | none | Mac metadata file removed |

## Package Call Names

The new reusable code is intended to be called through MATLAB package namespaces, for example:

```matlab
vleo.aero.compute_aero_forces(...)
vleo.analysis.state_to_observation(...)
vleo.control.attitude_pd_controller(...)
vleo.dynamics.sat_dynamics_nonlinear(...)
vleo.gui.open_simulation_gui()
vleo.viz.animate_results(...)
```
