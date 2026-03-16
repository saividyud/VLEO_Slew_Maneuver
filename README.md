# VLEO Slew Maneuver

Computationally determining the feasibility of using ion thrusters for slew maneuvers in Very Low Earth Orbit (VLEO).

## Project Structure

```text
.
├── src/
│   └── +vleo/          # Namespaced MATLAB package code
│       ├── +aero/      # Aerodynamic and surface interaction models
│       ├── +analysis/  # Observation and orbit analysis helpers
│       ├── +control/   # Controller implementations and resolution helpers
│       ├── +dynamics/  # Spacecraft dynamics models
│       ├── +gui/       # GUI workflows and parameter dialogs
│       ├── +util/      # Shared path, quaternion, and config utilities
│       └── +viz/       # Animation and plotting utilities
├── scripts/            # Runnable top-level MATLAB scripts
├── examples/           # Reference examples and demos
│   └── legacy/         # Older exploratory examples kept for reference
├── tests/              # Verification scripts and regression tests
│   ├── legacy/         # Incomplete or archived tests
│   └── templates/      # Starter template for new MATLAB tests
├── data/
│   ├── geometry/       # Meshes and geometry assets used by models
│   ├── generated/      # Generated intermediate data
│   └── legacy/         # Legacy CSV/data files kept for compatibility
├── tools/              # Non-MATLAB helper tooling (GIF rendering, build helpers)
├── simulations/        # Saved simulation outputs and generated artifacts
├── assets/             # Exported media and legacy result files
├── lib/                # External libraries (HPOP, SGP4, etc.) - DO NOT MODIFY
├── docs/               # Documentation and theory
└── workspaces/         # Personal development sandboxes
```

## Setup

1. Open MATLAB in the repository root
2. Run the setup script:
   ```matlab
    setup_project
    ```
3. Call package functions with their namespace, for example:
   ```matlab
   vleo.gui.open_simulation_gui()
   odefun = @vleo.dynamics.sat_dynamics_nonlinear;
   ```

`setup_project` is silent by default. Use `setup_project(true)` to print added paths for debugging.

## Main Entry Points

- Open the GUI:
  ```matlab
  run('scripts/open_simulation_gui.m')
  ```
- Run the volcano observation scenario:
  ```matlab
  run('scripts/volcano_tracking_demo.m')
  ```
- Run the aerodynamic regression test:
  ```matlab
  results = test_compute_aero_forces();
  ```
- Run the nonlinear dynamics smoke test:
  ```matlab
  test_sat_dynamics_nonlinear
  ```
- Render the latest `control_test3` animation data to GIF:
  ```bash
  python tools/render_gif.py
  ```

## Project Standards

See [PROJECT_STANDARDS.md](PROJECT_STANDARDS.md) for complete guidelines.

### Quick Reference

| Item | Standard |
|------|----------|
| ECI Frame | EME2000 (J2000) |
| Units | SI (meters, seconds, radians) |
| Quaternion | Scalar-first `[q0; q1; q2; q3]` |
| Package Naming | `vleo.<domain>.<function>` |
| Variable Naming | `r_sat_eci`, `v_sat_eci`, `q_eci_to_body` |

### Compliance Checklist

Before merging reusable code into `src/+vleo/`:
- [ ] Test file in `tests/`
- [ ] Example file in `examples/` for each new custom repo function
- [ ] Theory documented in `docs/theory.tex`
- [ ] Runnable scripts placed in `scripts/` instead of `src/`
- [ ] No magic numbers
- [ ] SI units throughout

## Contributors

- Thanadis Charoenrujijin
- Sai Vidyud Senthil Nathan
- Connor Lashley

## License

This repository is licensed under the GNU General Public License, version 3 or later (`GPL-3.0-or-later`). See `LICENSE`.

Parts of the aerodynamic modeling workflow are adapted from ADBSat, which is also distributed under GPL-3.0-compatible terms. See `THIRD_PARTY_LICENSES.md` for attribution and compliance notes.
