# VLEO Slew Maneuver

Exploring the feasibility of using ion thrusters for slew maneuvers in Very Low Earth Orbit (VLEO).

## Project Structure

The repository is organized as follows:

```text
.
├── src/                # Core source code and utilities
│   ├── utils/          # Mathematical conversions (e.g., DCMfromQ)
│   └── dynamics/       # Physics models (e.g., J2, control_torques)
├── lib/                # External libraries (ADBSat, HPOP, SGP4)
├── tests/              # Unit tests and benchmarks
├── workspaces/         # User-specific working directories
│   ├── Nill/
│   └── Sai/
└── docs/               # Documentation and reports
```

## Setup and Usage

To ensure all functions and libraries are accessible, you **must** configure the MATLAB path before running any scripts.

1.  Open MATLAB in the root directory of this repository.
2.  Run the setup script in the Command Window:
    ```matlab
    setup_project
    ```
    This script will add all necessary folders (`src`, `lib`, `tests`, etc.) to your MATLAB path.

3.  You can now run any script or function (e.g., from `tests/` or your workspace) without worrying about file locations.