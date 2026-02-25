# Ion Thruster Slew Feasibility Report

## Scenario

- Simulation case: `simulations/My_Simulation_torque_breakdown.csv`
- Model branch: nonlinear attitude dynamics with aerodynamic and control torques enabled in GUI flow
- Purpose: check whether ion-thruster-class thrust levels can provide the control torque seen in this scenario

## 1) What the simulation shows

From `simulations/My_Simulation_torque_breakdown.csv`:

- Peak control torque magnitude: `|tau_control|max = 2.1650635e-3 N m` at `t = 0 s`
- Peak aerodynamic torque magnitude: `|tau_aero|max = 3.5600937e-7 N m` at `t = 49 s`
- Aero/control peak ratio: `|tau_aero|max / |tau_control|max = 1.644e-4` (about `0.0164%`)
- Control torque magnitude at end of run: `3.6291109e-7 N m` (about `5966x` below initial peak)

Interpretation: in this run, control torque demand is dominated by the initial slew transient, while aerodynamic torque is much smaller.

## 2) Required thrust equivalent for the simulated peak torque

Using `tau = F * r`, the equivalent thrust required is:

`F_required = tau_required / r`

Assuming a representative moment arm of `r = 0.29 m` (half of the `0.58 m` characteristic size used in inertia approximation in `src/controllers/control_torques.m`):

- `F_required = 2.1650635e-3 / 0.29 = 7.466e-3 N = 7.466 mN`

Sensitivity to moment arm:

| Assumed moment arm r (m) | Required thrust (mN) |
|---:|---:|
| 0.29 | 7.466 |
| 0.20 | 10.825 |
| 0.15 | 14.434 |
| 0.10 | 21.651 |

## 3) Public ion-thruster thrust references (full links)

### Gridded ion thrusters

| Thruster | Type | Power range | Reported thrust | Notes | Source |
|---|---|---|---|---|---|
| NEXT-C | Gridded ion | 0.5-6.9 kW | 25-235 mN | NASA guidebook key requirements | [1] |
| QinetiQ T7 | Gridded ion | around 2.7-6.4 kW operating points shown | up to 250 mN in Table 1 operating point set | Performance target/expected map in IEPC paper | [2] |
| ArianeGroup RIT 2X | RF gridded ion | 2-8 kW | 70-260 mN | 2024 IEPC paper | [3] |
| QinetiQ T5 | Gridded ion | 700 W class | about 25 mN max per thruster | Flight heritage on GOCE/Artemis family | [4] |

### Hall thruster context (often grouped as "ion" in non-specialist discussions)

| Thruster | Type | Power range | Reported thrust | Source |
|---|---|---|---|---|
| HERMeS/AEPS development line | Hall effect | 6.3-12.5 kW | 396-613 mN | [5] |

## 4) Feasibility statement for this simulation scenario

Using the simulated peak required torque (`2.165e-3 N m`):

- At `r = 0.29 m`, required thrust is `7.466 mN`.
- Even a small gridded ion class at `25 mN` provides about `3.35x` thrust margin versus this requirement.
- Higher-power ion thrusters (`~235-260 mN`) provide roughly `31x-35x` margin for the same assumed arm.

Therefore, for this specific simulated torque profile, an ion-thruster-class actuator can provide sufficient moment magnitude for slew authority, assuming suitable actuator geometry and control architecture.

## 5) Important engineering caveats

- This check is magnitude-based (`tau = F * r`) and does not include actuator dynamics (throttle bandwidth, gimbal dynamics, valve delays, minimum stable thrust, duty cycling, plume constraints).
- Slew control requires signed torque about all commanded axes; practical implementation needs multiple thrusters and/or gimbal strategy.
- Traditional attitude control for agile slews is often reaction wheels/CMGs, with electric propulsion primarily for translation; this report only verifies torque magnitude feasibility for the simulated case.

## References

[1] NASA New Frontiers NEXT-C Guidebook (PDF): https://newfrontiers.larc.nasa.gov/NF4/PDF_FILES/NEXT-C_New_Frontiers_Guidebook_20161230_REV5.pdf

[2] J.-P. Luna et al., "T7 Thruster Design and Performance," IEPC 2019 (PDF): https://electricrocket.org/2019/356.pdf

[3] J.-P. Porst et al., "RIT 2X and RIT Family Radio Frequency Ion Thruster Development and Product Evolution," IEPC 2024 (PDF): https://electricrocket.org/IEPC_2024/602.pdf

[4] P. N. Randall et al., "QinetiQ T5 based Electric Propulsion System and Architectural Options for Future Applications," IEPC 2017 (PDF): https://electricrocket.org/IEPC/IEPC_2017_170.pdf

[5] R. Hofer et al., "Completing the Development of the 12.5 kW Hall Effect Rocket with Magnetic Shielding," IEPC 2019 (PDF): https://electricrocket.org/2019/193.pdf
