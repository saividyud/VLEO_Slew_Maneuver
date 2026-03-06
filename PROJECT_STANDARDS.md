# VLEO Slew Maneuver Project Standards

**Project Goal**: Computationally determine if an ion thruster is a viable option for VLEO slew maneuvers.

**Team Size**: 3 contributors

**Last Updated**: 2025-12-30

---

## Table of Contents
1. [Architecture](#1-architecture)
2. [Compliance](#2-compliance)
3. [Standards](#3-standards)
4. [Code Style](#4-code-style)
5. [Git Workflow](#5-git-workflow)
6. [Review Checklist](#6-review-checklist)

---

## 1. Architecture

### 1.1 Modular Structure
The codebase **SHALL** maintain a modular structure enabling simultaneous independent development.

**Directory Structure**:
```
VLEO_Slew_Maneuver/
├── src/                    # Core reusable source code
│   ├── controllers/        # Control laws and torque models
│   ├── dynamics/           # Physics models and assets
│   ├── gui/                # GUI and visualization code
│   ├── propagators/        # Orbit propagation methods
│   └── analysis/           # Analysis and post-processing functions
├── lib/                    # External libraries (HPOP, SGP4, etc.) - DO NOT MODIFY
├── tests/                  # Verification and unit tests
├── examples/               # Examples for custom repo functions and demos
├── docs/                   # Documentation and theory
├── workspaces/             # Personal development sandboxes
│   ├── Nill/
│   ├── Sai/
│   └── [Name]/
└── data/                   # Shared data files (TLE, ephemeris, etc.)
```

### 1.2 Module Independence Criteria
A module is considered **properly decoupled** when:
- [ ] It can be tested in isolation (no dependencies on workspace-specific code)
- [ ] It has a clearly defined interface (inputs/outputs documented)
- [ ] It does not rely on global variables or workspace state
- [ ] Changes to the module do not break unrelated modules
- [ ] Another team member can use it without reading the implementation

### 1.3 Code Promotion Path
```
workspaces/[Name]/  →  (review)  →  src/  →  (merge to main)
     ↓
   Personal experimentation allowed
   No standards enforcement required
```

**Rule**: Code in `src/` MUST meet all compliance requirements. Code in `workspaces/` is for experimentation and does not require compliance.

### 1.4 Dependency Management
- External libraries **SHALL** reside in `lib/` and remain unmodified
- If a library function needs modification, create a wrapper in `src/` with clear documentation of changes
- All dependencies **SHALL** be documented in `setup_project.m`

---

## 2. Compliance

### 2.1 Testing Requirements
Every function in `src/` **MUST** have:

1. **Verification Test** (`tests/test_<function_name>.m`)
   - Tests against known analytical solutions or published data
   - Tests edge cases (zero inputs, singularities, limits)
   - Returns PASS/FAIL status

Every **custom repo function** in `src/` **MUST** also have:

2. **Example Use Case** (`examples/example_<function_name>.m`)
   - Demonstrates typical usage
   - Contains explanatory comments
   - Produces interpretable output (plots or printed values)

**Exception**: Direct use of MATLAB or Aerospace Toolbox functions does not create a new example-file requirement.

### 2.2 No Redundant Functions
**DO NOT** create custom functions if MATLAB built-in or well-established library functions exist.

| Instead of...           | Use...                          |
|-------------------------|---------------------------------|
| Custom cross product    | `cross(a, b)`                   |
| Custom matrix inverse   | `inv(A)` or `A \ b`             |
| Custom norm             | `norm(v)`                       |
| Custom quaternion mult  | Aerospace Toolbox or document why custom needed |

**If a custom implementation is necessary**: Document the reason in the function header (performance, edge case handling, licensing).

### 2.3 Variable Naming
All variables **MUST** have self-explanatory names understandable by anyone with aerospace background.

**Good**:
```matlab
orbital_period_sec = 2 * pi * sqrt(semi_major_axis_m^3 / mu_earth);
quaternion_inertial_to_body = [q0; q1; q2; q3];
angular_velocity_body_rad_s = [omega_x; omega_y; omega_z];
```

**Bad**:
```matlab
T = 2 * pi * sqrt(a^3 / mu);  % What is T? What is a?
q = [q0; q1; q2; q3];          % q of what?
w = [wx; wy; wz];              % w in what frame?
```

**Exception**: Loop indices (`i`, `j`, `k`) and temporary variables in limited scope.

### 2.4 Theory Documentation
All physics, engineering, and mathematical theory used in `src/` functions **MUST** be documented in `docs/theory.md`.

Each theory entry **SHALL** include:
- Name/title of the concept
- Mathematical formulation
- Source citation (textbook, paper, or URL)
- List of functions that implement this theory

**Example entry**:
```markdown
### J2 Perturbation
The J2 perturbation accounts for Earth's oblateness...

**Equation**:
$$a_{J2} = \frac{3}{2} J_2 \frac{\mu R_E^2}{r^4} \begin{bmatrix} ... \end{bmatrix}$$

**Source**: Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications, 4th ed., p. 596

**Implemented in**: Aerospace Toolbox `gravityzonal` as used in `src/Sat_template.m` and `src/sat_template_gui.m`
```

### 2.5 No Magic Numbers
Hard-coded numerical values **SHALL NOT** appear in code unless:

1. **Well-known physical constants** (document source):
   ```matlab
   mu_earth = 3.986004418e14;  % [m^3/s^2] Earth gravitational parameter (WGS84)
   R_earth = 6378137;          % [m] Earth equatorial radius (WGS84)
   J2 = 1.08263e-3;            % [-] Earth J2 coefficient (EGM96)
   ```

2. **Mathematical constants**:
   ```matlab
   pi, exp(1), sqrt(2)         % Use MATLAB built-ins
   ```

3. **Algorithm-specific values** with documented justification:
   ```matlab
   tolerance = 1e-12;  % Convergence tolerance for Newton-Raphson iteration
                       % Chosen based on machine precision (eps ~ 2e-16)
   ```

**For configurable parameters**: Pass as function arguments or define in a configuration structure.

```matlab
% BAD
altitude = 300e3;  % Why 300 km? Is this always valid?

% GOOD
function result = analyze_orbit(altitude_m)
    % altitude_m: Orbital altitude [m], valid range [200e3, 2000e3] for VLEO/LEO
```

### 2.6 Error Handling
Functions **SHALL**:
- Validate input dimensions and types
- Return meaningful error messages
- Handle edge cases gracefully (circular orbits, singularities, etc.)

```matlab
function oe = state_to_elements(r, v, mu)
    % Input validation
    if length(r) ~= 3 || length(v) ~= 3
        error('state_to_elements:InvalidInput', 'r and v must be 3-element vectors');
    end
    if mu <= 0
        error('state_to_elements:InvalidInput', 'mu must be positive');
    end
    % ...
end
```

---

## 3. Standards

### 3.1 Reference Frames

#### ECI Frame: EME2000 (Earth Mean Equator and Equinox of J2000)
- **Origin**: Earth center of mass
- **X-axis**: Points to vernal equinox at J2000.0 epoch (2000-01-01 12:00:00 TT)
- **Z-axis**: Points to celestial north pole (perpendicular to mean equator at J2000.0)
- **Y-axis**: Completes right-handed system
- **Property**: Non-rotating (inertial for most LEO applications)

**Rationale**: EME2000 is the de facto standard for astrodynamics, compatible with HPOP, SGP4, and most ephemeris data.

#### Satellite Body Frame
- **Origin**: Satellite center of mass
- **X-axis (x_body)**: Perpendicular to the 6-face side (primary thrust axis)
- **Y-axis (y_body)**: Perpendicular to the 3-face side
- **Z-axis (z_body)**: Perpendicular to the 2-face side
- **Property**: Right-handed coordinate system

```
        z_body
           ↑
           |    2-face
           |   ┌─────┐
           |  /      /|
           | /  6   / | 3-face
           |/______/  |
           └──────┴──→ y_body
          /
         ↙ x_body
```

#### LVLH Frame (Local Vertical Local Horizontal)
- **Origin**: Satellite center of mass
- **Z-axis (nadir)**: Points toward Earth center (`-r/|r|`)
- **Y-axis (cross-track)**: Opposite to orbital angular momentum (`-h/|h|`)
- **X-axis (velocity)**: Completes right-handed system (approximately velocity direction for circular orbits)

### 3.2 Naming Conventions

#### Vector Notation
Format: `<quantity>_<of_what>_<in_frame>`

| Symbol | Meaning |
|--------|---------|
| `r`    | Position vector |
| `v`    | Velocity vector |
| `a`    | Acceleration vector |
| `omega` or `w` | Angular velocity vector |
| `q`    | Quaternion |
| `T`    | Torque vector |
| `F`    | Force vector |
| `h`    | Angular momentum vector |

**Examples**:
```matlab
r_sat_ECI       % Position of satellite in ECI frame [m]
v_sat_ECI       % Velocity of satellite in ECI frame [m/s]
omega_body      % Angular velocity of body frame w.r.t. inertial, expressed in body [rad/s]
q_ECI_to_body   % Quaternion rotating from ECI to body frame
T_control_body  % Control torque expressed in body frame [N·m]
r_sun_ECI       % Position of Sun in ECI frame [m]
```

#### Transformation Notation
Format: `<type>_<from>_to_<to>` or `<type>_<from>_<to>`

```matlab
DCM_ECI_to_body    % Direction Cosine Matrix from ECI to body
q_ECI_to_LVLH      % Quaternion from ECI to LVLH
R_1                % Elementary rotation about axis 1
```

#### Scalar Quantities
Include units in name or comment:
```matlab
altitude_m = 400e3;           % Altitude in meters
period_sec = 5400;            % Orbital period in seconds
inclination_rad = deg2rad(51.6);  % Inclination in radians
mass_kg = 100;                % Satellite mass in kilograms
```

### 3.3 Units
All computations **SHALL** use SI units unless explicitly documented.

| Quantity | SI Unit | Symbol |
|----------|---------|--------|
| Length | meters | m |
| Mass | kilograms | kg |
| Time | seconds | s |
| Angle | radians | rad |
| Force | Newtons | N |
| Torque | Newton-meters | N·m |
| Angular velocity | rad/s | rad/s |
| Specific impulse | seconds | s |
| Thrust | Newtons | N |

**Conversion functions**: Use MATLAB's `deg2rad()`, `rad2deg()` for angles. Document conversions explicitly.

### 3.4 Quaternion Convention
- **Format**: `q = [q0; q1; q2; q3]` where `q0` is the scalar part
- **Normalization**: `|q| = 1`
- **Convention**: Hamilton multiplication, passive rotation (frame transformation)
- **Identity**: `q = [1; 0; 0; 0]` represents no rotation

### 3.5 Time Systems
- **Primary**: UTC for input/output, TT (Terrestrial Time) for precise calculations
- **Epoch format**: Julian Date (JD) or Modified Julian Date (MJD) for internal calculations
- **Display format**: ISO 8601 (`YYYY-MM-DDTHH:MM:SS.sssZ`)

### 3.6 Orbital Elements
Classical Keplerian elements in order:
```matlab
oe = [a; e; i; RAAN; omega; f];
%      │  │  │   │      │    └── True anomaly [rad]
%      │  │  │   │      └─────── Argument of periapsis [rad]
%      │  │  │   └────────────── Right ascension of ascending node [rad]
%      │  │  └─────────────────── Inclination [rad]
%      │  └────────────────────── Eccentricity [-]
%      └───────────────────────── Semi-major axis [m]
```

---

## 4. Code Style

### 4.1 Function Header Template
```matlab
function [output1, output2] = function_name(input1, input2, options)
%FUNCTION_NAME Brief one-line description
%
%   [OUT1, OUT2] = FUNCTION_NAME(IN1, IN2) detailed description of what
%   the function does and how it works.
%
%   [OUT1, OUT2] = FUNCTION_NAME(IN1, IN2, OPTIONS) description of optional
%   parameters.
%
%   Inputs:
%       input1  - Description, units, dimensions [m], [3x1]
%       input2  - Description, units, dimensions [rad], scalar
%       options - (optional) Description of optional parameters
%
%   Outputs:
%       output1 - Description, units, dimensions [m/s], [3x1]
%       output2 - Description, units, dimensions [-], scalar
%
%   Theory:
%       Brief explanation or reference to docs/theory.md section
%
%   Example:
%       [v, a] = function_name([1e6; 0; 0], 0.5);
%
%   See also: RELATED_FUNCTION1, RELATED_FUNCTION2
%
%   References:
%       [1] Author, "Title", Year, Page
%
%   Author: Name
%   Date: YYYY-MM-DD
%   Version: 1.0

% Implementation...
end
```

### 4.2 File Organization
```matlab
% 1. Function signature
function output = my_function(input)

% 2. Header documentation (as above)

% 3. Input validation
arguments  % (R2019b+) or manual validation
    input (3,1) double
end

% 4. Constants (if needed)
MU_EARTH = 3.986004418e14;  % [m^3/s^2]

% 5. Main computation
% ... implementation ...

% 6. Output formatting
output = result;

end

% 7. Local helper functions (if needed)
function y = helper(x)
    % ...
end
```

### 4.3 Commenting
- **Code comments**: Explain *why*, not *what*
- **Section dividers**: Use `%% Section Name` for logical sections
- **TODO/FIXME**: Use format `% TODO(name): description`

---

## 5. Git Workflow

### 5.1 Branch Naming
```
main                    # Protected, always working
feature/<description>   # New features
fix/<description>       # Bug fixes
experiment/<name>/<desc> # Personal experiments (can be messy)
```

### 5.2 Commit Messages
Format: `<type>: <description>`

Types:
- `Add`: New feature or file
- `Fix`: Bug fix
- `Update`: Enhancement to existing feature
- `Refactor`: Code restructure without behavior change
- `Docs`: Documentation only
- `Test`: Test additions or fixes

**Examples**:
```
Add: J3 perturbation model in dynamics
Fix: Handle circular-orbit singularity in state-to-elements conversion
Docs: Update theory.md with quaternion kinematics derivation
Test: Add quaternion conversion verification against Vallado examples
```

### 5.3 Pull Request Requirements
Before merging to `main`:
- [ ] All tests pass
- [ ] Code follows naming conventions
- [ ] Functions have proper headers
- [ ] Theory documented (if applicable)
- [ ] At least one team member reviewed

---

## 6. Review Checklist

### For Code Authors
Before requesting review:
- [ ] Function has complete header documentation
- [ ] Variable names follow conventions (Section 3.2)
- [ ] No magic numbers (Section 2.5)
- [ ] Units are SI (Section 3.3)
- [ ] Test file exists in `tests/`
- [ ] Example file exists in `examples/` for each new custom repo function
- [ ] Theory added to `docs/theory.md` (if new physics)
- [ ] Code runs without errors
- [ ] Edge cases handled

### For Code Reviewers
When reviewing:
- [ ] Can I understand this code without asking questions?
- [ ] Are the variable names clear?
- [ ] Is the theory correct and properly cited?
- [ ] Do the tests cover edge cases?
- [ ] Does this duplicate existing functionality?
- [ ] Would I be comfortable maintaining this code?

---

## Appendix A: Physical Constants

Standard values to use throughout the project (WGS84/EGM96 unless noted):

```matlab
% Gravitational parameters
mu_earth = 3.986004418e14;     % [m^3/s^2] Earth
mu_sun = 1.32712440018e20;     % [m^3/s^2] Sun
mu_moon = 4.9028695e12;        % [m^3/s^2] Moon

% Earth parameters
R_earth = 6378137;             % [m] Equatorial radius
f_earth = 1/298.257223563;     % [-] Flattening
J2 = 1.08263e-3;               % [-] J2 coefficient
J3 = -2.5327e-6;               % [-] J3 coefficient
omega_earth = 7.2921150e-5;    % [rad/s] Earth rotation rate

% Other
c_light = 299792458;           % [m/s] Speed of light
AU = 1.495978707e11;           % [m] Astronomical unit
```

---

## Appendix B: Ion Thruster Reference Parameters

Typical ranges for VLEO ion thruster analysis:

| Parameter | Typical Range | Notes |
|-----------|---------------|-------|
| Thrust | 1 - 100 mN | Hall thruster / ion engine range |
| Specific impulse | 1000 - 4000 s | Higher than chemical propulsion |
| Power | 100 W - 10 kW | Depends on thrust level |
| Propellant | Xenon, Krypton | Xenon most common |
| Slew rate requirement | 1 - 10 deg/s | Mission dependent |
| VLEO altitude | 200 - 450 km | High drag environment |

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2025-12-30 | Team | Initial draft |
