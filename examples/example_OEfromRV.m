%% Example: Converting Position/Velocity to Orbital Elements
% This example demonstrates how to use OEfromRV to convert ECI state vectors
% to classical Keplerian orbital elements.
%
% Author: Team
% Date: 2025-12-30

%% Setup
clear; clc;

% Add project paths
run('../setup_project.m');

%% Define constants
mu_earth = 3.986004418e14;  % [m^3/s^2] Earth gravitational parameter (WGS84)

%% Example 1: ISS-like orbit
% ISS orbits at approximately 420 km altitude, 51.6 deg inclination

% Position and velocity in ECI frame
r_sat_ECI = [6787746.891;   % [m] x-component
             1234567.890;   % [m] y-component
             2345678.901];  % [m] z-component

v_sat_ECI = [-1234.567;     % [m/s] vx-component
              6789.012;     % [m/s] vy-component
              3456.789];    % [m/s] vz-component

% Convert to orbital elements
oe = OEfromRV(r_sat_ECI, v_sat_ECI, mu_earth);

% Display results
fprintf('=== Example 1: Position/Velocity to Orbital Elements ===\n\n');
fprintf('Input State Vector (ECI):\n');
fprintf('  Position: [%.3f, %.3f, %.3f] km\n', r_sat_ECI'/1e3);
fprintf('  Velocity: [%.3f, %.3f, %.3f] km/s\n\n', v_sat_ECI'/1e3);

fprintf('Output Orbital Elements:\n');
fprintf('  Semi-major axis (a):    %.3f km\n', oe(1)/1e3);
fprintf('  Eccentricity (e):       %.6f\n', oe(2));
fprintf('  Inclination (i):        %.3f deg\n', rad2deg(oe(3)));
fprintf('  RAAN (Omega):           %.3f deg\n', rad2deg(oe(4)));
fprintf('  Arg. of periapsis (w):  %.3f deg\n', rad2deg(oe(5)));
fprintf('  True anomaly (f):       %.3f deg\n\n', rad2deg(oe(6)));

%% Example 2: Circular equatorial orbit (edge case)
% Tests handling of special cases: i ≈ 0, e ≈ 0

altitude_m = 400e3;  % [m] 400 km altitude
r_circular = 6378137 + altitude_m;  % [m] orbital radius
v_circular = sqrt(mu_earth / r_circular);  % [m/s] circular velocity

r_equatorial_ECI = [r_circular; 0; 0];      % [m] on x-axis
v_equatorial_ECI = [0; v_circular; 0];      % [m/s] in y-direction

oe_circular = OEfromRV(r_equatorial_ECI, v_equatorial_ECI, mu_earth);

fprintf('=== Example 2: Circular Equatorial Orbit (Edge Case) ===\n\n');
fprintf('Input: Circular orbit at %.0f km altitude\n\n', altitude_m/1e3);
fprintf('Output Orbital Elements:\n');
fprintf('  Semi-major axis (a):    %.3f km\n', oe_circular(1)/1e3);
fprintf('  Eccentricity (e):       %.6e (near-circular)\n', oe_circular(2));
fprintf('  Inclination (i):        %.3f deg (equatorial)\n', rad2deg(oe_circular(3)));

%% Example 3: Round-trip verification
% Convert OE -> RV -> OE to verify consistency

fprintf('\n=== Example 3: Round-Trip Verification ===\n\n');

% Define orbital elements directly
a_test = 7000e3;           % [m] semi-major axis
e_test = 0.001;            % [-] eccentricity (nearly circular)
i_test = deg2rad(45);      % [rad] inclination
RAAN_test = deg2rad(30);   % [rad] RAAN
omega_test = deg2rad(60);  % [rad] argument of periapsis
f_test = deg2rad(90);      % [rad] true anomaly

oe_original = [a_test; e_test; i_test; RAAN_test; omega_test; f_test];

% Convert to RV
[r_test, v_test] = RVfromOE(oe_original, mu_earth);

% Convert back to OE
oe_recovered = OEfromRV(r_test, v_test, mu_earth);

% Compute errors
error_a = abs(oe_recovered(1) - oe_original(1));
error_e = abs(oe_recovered(2) - oe_original(2));
error_i = abs(oe_recovered(3) - oe_original(3));

fprintf('Round-trip errors:\n');
fprintf('  Semi-major axis error: %.3e m\n', error_a);
fprintf('  Eccentricity error:    %.3e\n', error_e);
fprintf('  Inclination error:     %.3e rad\n', error_i);

if error_a < 1e-6 && error_e < 1e-10 && error_i < 1e-10
    fprintf('\n  [PASS] Round-trip conversion successful!\n');
else
    fprintf('\n  [WARN] Errors exceed expected tolerance\n');
end

%% Notes
% - This function follows the EME2000 (J2000) reference frame convention
% - Units: SI (meters, seconds, radians)
% - See docs/theory.md Section 2.1 for mathematical derivation
% - See also: RVfromOE (inverse function)
