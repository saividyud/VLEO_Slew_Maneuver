# VLEO Slew Maneuver

Computationally determining the feasibility of using ion thrusters for slew maneuvers in Very Low Earth Orbit (VLEO).

## Project Structure

```text
.
├── src/                # Core source code (must meet compliance standards)
│   ├── utils/          # Mathematical conversions (DCM, quaternions, orbital elements)
│   └── dynamics/       # Physics models (J2, control torques)
├── lib/                # External libraries (HPOP, SGP4, ADBSat) - DO NOT MODIFY
├── tests/              # Unit tests and verification scripts
├── examples/           # Example use cases for each function
├── docs/               # Documentation and theory
│   └── theory.md       # Mathematical foundations and citations
├── data/               # Shared data files (TLE, ephemeris)
└── workspaces/         # Personal development sandboxes
    ├── Nill/
    └── Sai/
```

## Setup

1. Open MATLAB in the repository root
2. Run the setup script:
   ```matlab
   setup_project
   ```
3. All paths are now configured

## Project Standards

See [PROJECT_STANDARDS.md](PROJECT_STANDARDS.md) for complete guidelines.

### Quick Reference

| Item | Standard |
|------|----------|
| ECI Frame | EME2000 (J2000) |
| Units | SI (meters, seconds, radians) |
| Quaternion | Scalar-first `[q0; q1; q2; q3]` |
| Naming | `r_sat_ECI`, `v_sat_ECI`, `q_ECI_to_body` |

### Compliance Checklist

Before merging code to `src/`:
- [ ] Test file in `tests/`
- [ ] Example file in `examples/`
- [ ] Theory documented in `docs/theory.md`
- [ ] No magic numbers
- [ ] SI units throughout

## Contributors

- Nill
- Sai Vidyud
