# Theoretical Foundations

This document contains the mathematical and physical theory underlying the VLEO Slew Maneuver codebase. Each section corresponds to implementations in `src/`.

---

## Table of Contents
1. [Coordinate Transformations](#1-coordinate-transformations)
2. [Orbital Mechanics](#2-orbital-mechanics)
3. [Attitude Dynamics](#3-attitude-dynamics)
4. [Perturbations](#4-perturbations)
5. [Propulsion](#5-propulsion)
6. [References](#6-references)

---

## 1. Coordinate Transformations

### 1.1 Quaternion to Direction Cosine Matrix

A quaternion `q = [q0; q1; q2; q3]` (scalar-first convention) can be converted to a Direction Cosine Matrix (DCM) using:

```
       ┌                                                              ┐
       │ q0² + q1² - q2² - q3²    2(q1q2 + q0q3)      2(q1q3 - q0q2)  │
DCM =  │ 2(q1q2 - q0q3)          q0² - q1² + q2² - q3²  2(q2q3 + q0q1) │
       │ 2(q1q3 + q0q2)          2(q2q3 - q0q1)      q0² - q1² - q2² + q3² │
       └                                                              ┘
```

**Source**: Wertz, J. R., "Spacecraft Attitude Determination and Control", 1978, Eq. 12-12

**Implemented in**: `src/utils/DCMfromQ.m`

---

### 1.2 Direction Cosine Matrix to Quaternion

The inverse transformation uses Shepperd's method for numerical stability:

1. Compute trace and diagonal elements
2. Find maximum of {trace, C11, C22, C33}
3. Use corresponding extraction formula to avoid division by small numbers

**Source**: Shepperd, S.W., "Quaternion from Rotation Matrix", Journal of Guidance and Control, Vol. 1, No. 3, 1978

**Implemented in**: `src/utils/QfromDCM.m`

---

### 1.3 Quaternion and Principal Axis/Angle

**Quaternion from axis-angle**:
```
q0 = cos(θ/2)
q1 = e1 * sin(θ/2)
q2 = e2 * sin(θ/2)
q3 = e3 * sin(θ/2)
```

Where `e = [e1; e2; e3]` is the unit rotation axis and `θ` is the rotation angle.

**Axis-angle from quaternion**:
```
θ = 2 * acos(q0)
e = [q1; q2; q3] / sin(θ/2)    if θ ≠ 0
e = [1; 0; 0]                   if θ = 0 (arbitrary)
```

**Source**: Euler's rotation theorem, Kuipers, J., "Quaternions and Rotation Sequences", 1999

**Implemented in**: `src/utils/QfromPAT.m`, `src/utils/PATfromQ.m`

---

### 1.4 Fundamental Rotation Matrices

Elementary rotations about principal axes:

**Rotation about axis 1 (x-axis)**:
```
         ┌ 1    0        0     ┐
R1(θ) =  │ 0   cos(θ)   sin(θ) │
         └ 0  -sin(θ)   cos(θ) ┘
```

**Rotation about axis 2 (y-axis)**:
```
         ┌ cos(θ)   0  -sin(θ) ┐
R2(θ) =  │   0      1     0    │
         └ sin(θ)   0   cos(θ) ┘
```

**Rotation about axis 3 (z-axis)**:
```
         ┌ cos(θ)   sin(θ)   0 ┐
R3(θ) =  │-sin(θ)   cos(θ)   0 │
         └   0        0      1 ┘
```

**Source**: Schaub, H. and Junkins, J., "Analytical Mechanics of Space Systems", 4th ed., 2018, Ch. 3

**Implemented in**: `src/utils/FRE.m`

---

## 2. Orbital Mechanics

### 2.1 Orbital Elements from Position/Velocity

Given position `r` and velocity `v` in ECI frame, classical orbital elements are computed:

1. **Specific angular momentum**: `h = r × v`
2. **Node vector**: `n = [0; 0; 1] × h`
3. **Eccentricity vector**: `e_vec = ((v² - μ/r)r - (r·v)v) / μ`
4. **Semi-major axis**: `a = -μ / (2ε)` where `ε = v²/2 - μ/r`
5. **Inclination**: `i = acos(h_z / |h|)`
6. **RAAN**: `Ω = acos(n_x / |n|)`, adjusted for n_y sign
7. **Argument of periapsis**: `ω = acos(n·e / (|n||e|))`, adjusted for e_z sign
8. **True anomaly**: `f = acos(e·r / (|e||r|))`, adjusted for r·v sign

**Special cases**:
- Circular orbits (e ≈ 0): Use argument of latitude `u = ω + f`
- Equatorial orbits (i ≈ 0): Use longitude of periapsis `ϖ = Ω + ω`

**Source**: Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th ed., 2013, Algorithm 9

**Implemented in**: `src/utils/OEfromRV.m`

---

### 2.2 Position/Velocity from Orbital Elements

1. Compute position in perifocal frame (PQW):
   ```
   r_pqw = (a(1-e²)/(1+e*cos(f))) * [cos(f); sin(f); 0]
   v_pqw = sqrt(μ/(a(1-e²))) * [-sin(f); e+cos(f); 0]
   ```

2. Transform to ECI using rotation sequence:
   ```
   R = R3(-Ω) * R1(-i) * R3(-ω)
   r_ECI = R * r_pqw
   v_ECI = R * v_pqw
   ```

**Source**: Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th ed., 2013, Algorithm 10

**Implemented in**: `src/utils/RVfromOE.m`

---

## 3. Attitude Dynamics

### 3.1 Quaternion Kinematics

The time derivative of the quaternion representing body orientation:

```
q̇ = (1/2) * Ω(ω) * q
```

Where `Ω(ω)` is the quaternion rate matrix:
```
        ┌  0   -ω1  -ω2  -ω3 ┐
Ω(ω) =  │ ω1    0    ω3  -ω2 │
        │ ω2  -ω3    0    ω1 │
        └ ω3   ω2  -ω1    0  ┘
```

And `ω = [ω1; ω2; ω3]` is the angular velocity of the body frame with respect to the inertial frame, expressed in body coordinates.

**Source**: Wertz, J. R., "Spacecraft Attitude Determination and Control", 1978, Section 12.1

**Implemented in**: `src/Sat_template.m`

---

### 3.2 Euler's Rotational Equations

For a rigid body with inertia tensor `J`:

```
J * ω̇ = -ω × (J * ω) + T_ext
```

Expanded for diagonal inertia tensor `J = diag([Jx, Jy, Jz])`:

```
Jx * ω̇x = (Jy - Jz) * ωy * ωz + Tx
Jy * ω̇y = (Jz - Jx) * ωz * ωx + Ty
Jz * ω̇z = (Jx - Jy) * ωx * ωy + Tz
```

**Source**: Schaub, H. and Junkins, J., "Analytical Mechanics of Space Systems", 4th ed., 2018, Ch. 4

**Implemented in**: `src/Sat_template.m`

---

### 3.3 PD Attitude Control Law

Quaternion-based feedback control:

```
T_control = -Kp * q_error(1:3) - Kd * ω_error
```

Where:
- `q_error`: Quaternion error between current and desired attitude
- `ω_error`: Angular velocity error
- `Kp`: Proportional gain matrix
- `Kd`: Derivative gain matrix

**Source**: Wie, B., "Space Vehicle Dynamics and Control", 2nd ed., 2008, Ch. 7

**Implemented in**: `src/dynamics/control_torques.m`

---

## 4. Perturbations

### 4.1 J2 Perturbation (Earth Oblateness)

The acceleration due to Earth's J2 zonal harmonic:

```
a_J2 = (3/2) * J2 * (μ/r²) * (R_E/r)² * [ (x/r)(5(z/r)² - 1) ]
                                        [ (y/r)(5(z/r)² - 1) ]
                                        [ (z/r)(5(z/r)² - 3) ]
```

Where:
- `J2 = 1.08263 × 10⁻³` (dimensionless)
- `R_E = 6378137 m` (Earth equatorial radius)
- `μ = 3.986004418 × 10¹⁴ m³/s²` (Earth gravitational parameter)
- `r = |r|` (distance from Earth center)
- `(x, y, z)` are ECI coordinates

**Source**: Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th ed., 2013, Eq. 8-25

**Implemented in**: `src/dynamics/J2.m`

---

### 4.2 J3 Perturbation

The acceleration due to Earth's J3 zonal harmonic:

```
a_J3 = (1/2) * J3 * (μ/r²) * (R_E/r)³ * [ (x/r)(10(z/r)³ - (15z/r)/r) ]
                                         [ (y/r)(10(z/r)³ - (15z/r)/r) ]
                                         [ (4(z/r)³ - (3z/r)/r)         ]
```

Where `J3 = -2.5327 × 10⁻⁶`

**Source**: Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th ed., 2013, Eq. 8-26

**Implemented in**: `src/dynamics/J2.m`

---

### 4.3 Atmospheric Drag (VLEO)

```
a_drag = -(1/2) * ρ * Cd * (A/m) * |v_rel| * v_rel
```

Where:
- `ρ`: Atmospheric density (varies with altitude, solar activity)
- `Cd`: Drag coefficient (~2.2 for typical spacecraft)
- `A`: Cross-sectional area
- `m`: Spacecraft mass
- `v_rel`: Velocity relative to atmosphere

**Atmospheric models**:
- Exponential (simple)
- NRLMSISE-00 (high fidelity)
- JB2008 (operational)

**Source**: Vallado, D. A., "Fundamentals of Astrodynamics and Applications", 4th ed., 2013, Ch. 8.6

**Implemented in**: `workspaces/Nill/testing_ground/atmospheric_drag.m` (under development)

---

### 4.4 Solar Radiation Pressure

```
a_SRP = -P_sr * Cr * (A/m) * (r_sun / |r_sun|)
```

Where:
- `P_sr ≈ 4.56 × 10⁻⁶ N/m²` at 1 AU
- `Cr`: Reflectivity coefficient (1.0 for absorption, 2.0 for specular reflection)
- `A`: Area normal to sun
- `r_sun`: Position vector from spacecraft to Sun

**Source**: Montenbruck, O. and Gill, E., "Satellite Orbits", 2000, Ch. 3.4

**Implemented in**: `workspaces/Nill/testing_ground/solar_radiation_pressure.m` (under development)

---

## 5. Propulsion

### 5.1 Ion Thruster Fundamentals

**Thrust equation**:
```
F = ṁ * v_e = ṁ * g0 * Isp
```

Where:
- `F`: Thrust [N]
- `ṁ`: Mass flow rate [kg/s]
- `v_e`: Exhaust velocity [m/s]
- `g0 = 9.80665 m/s²`: Standard gravity
- `Isp`: Specific impulse [s]

**Power-thrust relationship**:
```
P = (1/2) * ṁ * v_e² * (1/η)
```

Where `η` is the thruster efficiency.

**Source**: Goebel, D. and Katz, I., "Fundamentals of Electric Propulsion: Ion and Hall Thrusters", JPL, 2008

**Implemented in**: TBD

---

### 5.2 Slew Maneuver Analysis

For a rest-to-rest slew maneuver of angle `θ`:

**Minimum time (bang-bang control)**:
```
t_min = 2 * sqrt(θ * J / T_max)
```

**Minimum energy (smooth trajectory)**:
```
E_min = (4/3) * J * θ² / t³
```

Where:
- `J`: Moment of inertia about rotation axis
- `T_max`: Maximum available torque

**Source**: Wie, B., "Space Vehicle Dynamics and Control", 2nd ed., 2008, Ch. 7

**Implemented in**: TBD

---

## 6. References

### Textbooks

1. Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*, 4th ed. Microcosm Press.

2. Wertz, J. R. (1978). *Spacecraft Attitude Determination and Control*. Kluwer Academic Publishers.

3. Schaub, H. and Junkins, J. (2018). *Analytical Mechanics of Space Systems*, 4th ed. AIAA.

4. Wie, B. (2008). *Space Vehicle Dynamics and Control*, 2nd ed. AIAA.

5. Montenbruck, O. and Gill, E. (2000). *Satellite Orbits: Models, Methods, and Applications*. Springer.

6. Kuipers, J. (1999). *Quaternions and Rotation Sequences*. Princeton University Press.

7. Goebel, D. and Katz, I. (2008). *Fundamentals of Electric Propulsion: Ion and Hall Thrusters*. JPL/Wiley.

### Papers

8. Shepperd, S.W. (1978). "Quaternion from Rotation Matrix". *Journal of Guidance and Control*, Vol. 1, No. 3.

### Standards

9. IERS Conventions (2010). IERS Technical Note No. 36.

10. WGS84 (2014). NGA.STND.0036_1.0.0_WGS84.

---

## Document History

| Date | Author | Changes |
|------|--------|---------|
| 2025-12-30 | Team | Initial template with existing implementations |

