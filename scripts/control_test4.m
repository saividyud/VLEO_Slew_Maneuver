% Control Test 4 - Nadir Hold, Minimax Slew, and Volcano Tracking
% Holds nadir pointing before eruption, slews with a minimax torque profile,
% then matches the volcano-tracking attitude/rate so tracking can continue.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('includeAerodynamics', 'var')
    includeAerodynamics = true;
end

if ~exist('scenarioSeed', 'var')
    scenarioSeed = [];
end

if ~exist('targetMatchTime', 'var')
    targetMatchTime = [];
end

if ~exist('generateVisuals', 'var')
    generateVisuals = usejava('desktop');
end

clearvars -except includeAerodynamics scenarioSeed targetMatchTime generateVisuals
clc;

if ~(isscalar(includeAerodynamics) && (islogical(includeAerodynamics) || isnumeric(includeAerodynamics)))
    error('control_test4:InvalidAeroToggle', ...
        'includeAerodynamics must be a numeric or logical scalar.');
end
includeAerodynamics = logical(includeAerodynamics);

if isempty(scenarioSeed)
    rng('shuffle');
    scenarioSeedText = 'shuffle';
elseif isscalar(scenarioSeed) && isnumeric(scenarioSeed) && isfinite(scenarioSeed) && ...
        scenarioSeed >= 0 && scenarioSeed == floor(scenarioSeed)
    rng(double(scenarioSeed), 'twister');
    scenarioSeedText = num2str(scenarioSeed);
else
    error('control_test4:InvalidScenarioSeed', ...
        'scenarioSeed must be empty or a nonnegative integer scalar.');
end

if ~(isempty(targetMatchTime) || (isscalar(targetMatchTime) && isnumeric(targetMatchTime) && isfinite(targetMatchTime)))
    error('control_test4:InvalidTargetMatchTime', ...
        'targetMatchTime must be empty or a finite numeric scalar.');
end

if ~(isscalar(generateVisuals) && (islogical(generateVisuals) || isnumeric(generateVisuals)))
    error('control_test4:InvalidGenerateVisuals', ...
        'generateVisuals must be a numeric or logical scalar.');
end
generateVisuals = logical(generateVisuals);

% Parameters
params.mu = 3.986004418e14; % Earth gravitational parameter [m^3/s^2]
params.R_e = 6378137; % Earth radius [m]
params.bodyDims = [0.10; 0.20; 0.30]; % [x; y; z] body dimensions [m]
params.mass = 12.0; % Uniform 6U CubeSat mass [kg]
params.radius = 0.5 * norm(params.bodyDims); % Bounding radius [m]
params.J2 = 1.08263e-3; % J2 coefficient
params.I_CB = cuboid_inertia_matrix(params.mass, params.bodyDims); % [kg*m^2]
params.Kp_att = 10; % Proportional control gain [N*m]
params.omega_earth = 7.2921159e-5; % Earth rotation rate [rad/s]
params.includeAerodynamics = includeAerodynamics;
params.aeroSampleStep = 10; % Aerodynamic torque evaluation cadence [s]
params.aero.objFilePath = vleo.util.geometry_asset_path('6U CubeSat.obj');
params.aero.year = 2002;
params.aero.dayOfYear = 106;
params.aero.secondsOfDay = 43200;
params.aero.options.f107Average = 150;
params.aero.options.f107Daily = 150;
params.aero.options.magneticIndex = ones(1, 7) * 3;
params.aero.options.anomalousOxygen = 0;
params.aero.options.gsi_model = 'cook';
params.aero.options.alpha = 1;
params.aero.options.Tw = 300;
params.aero.options.enableShadow = 1;
params.aero.options.enableSolar = 0;
params.aero.options.sol_cR = 0.15;
params.aero.options.sol_cD = 0.25;

fprintf('=== VOLCANO OBSERVATION WITH MINIMAX ATTITUDE SLEW ===\n');
fprintf('Spacecraft: 6U CubeSat (10 x 20 x 30 cm, 12 kg, uniform density)\n');
fprintf('Aerodynamic disturbance torque estimate: %s\n', on_off_text(params.includeAerodynamics));
fprintf('Scenario seed: %s\n', scenarioSeedText);
fprintf('Principal inertias [kg m^2]: [%.4f, %.4f, %.4f]\n\n', diag(params.I_CB));

% Initial Conditions - Random position and velocity
altitude = 200e3; % 200 km
r_orbit = params.R_e + altitude;
v_circ = sqrt(params.mu / r_orbit);

% Random orbital plane orientation
inclination = rand() * pi;
RAAN = rand() * 2 * pi;
true_anomaly = rand() * 2 * pi;

R_RAAN = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
R_inc = [1, 0, 0; 0, cos(inclination), -sin(inclination); 0, sin(inclination), cos(inclination)];
R_ta = [cos(true_anomaly), -sin(true_anomaly), 0; sin(true_anomaly), cos(true_anomaly), 0; 0, 0, 1];
R_perifocal_to_ECI = R_RAAN * R_inc * R_ta;

r_perifocal = [r_orbit; 0; 0];
v_perifocal = [0; v_circ; 0];

r_eci = R_perifocal_to_ECI * r_perifocal;
v_eci = R_perifocal_to_ECI * v_perifocal;

% Initial attitude - nadir pointing
orbit_normal_eci = cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
z_body_eci = -r_eci / norm(r_eci);
x_body_eci = v_eci / norm(v_eci);
y_body_eci = cross(z_body_eci, x_body_eci);
R_Body_to_ECI = [x_body_eci, y_body_eci, z_body_eci];

q0 = dcm2quat(R_Body_to_ECI')';

omega_orbit_mag = v_circ / r_orbit;
omega_eci = omega_orbit_mag * orbit_normal_eci;

X0 = [r_eci; v_eci; q0; omega_eci];

% Calculate orbital period
orbital_period = 2 * pi * sqrt(r_orbit^3 / params.mu); % seconds

% Place volcano so it is visible at t = 0
r_sat_0 = X0(1:3);
r_projected = r_sat_0 / norm(r_sat_0) * params.R_e;

theta_max = acos(params.R_e / norm(r_sat_0));
theta = rand() * theta_max;
phi = rand() * 2 * pi;

axis1 = v_eci / norm(v_eci);
r_temp = r_projected * cos(theta) + ...
         cross(axis1, r_projected) * sin(theta) + ...
         axis1 * dot(axis1, r_projected) * (1 - cos(theta));

axis2 = r_projected / norm(r_projected);
r_volcano_0 = r_temp * cos(phi) + ...
              cross(axis2, r_temp) * sin(phi) + ...
              axis2 * dot(axis2, r_temp) * (1 - cos(phi));
r_volcano_0 = r_volcano_0 / norm(r_volcano_0) * params.R_e;

volcano_direction = r_volcano_0 / norm(r_volcano_0);
volcano_lat = rad2deg(asin(volcano_direction(3)));
volcano_lon = rad2deg(atan2(volcano_direction(2), volcano_direction(1)));

fprintf('Volcano erupts at t = 0\n');
fprintf('Volcano Location (at t=0): Lat=%.2f deg, Lon=%.2f deg\n', volcano_lat, volcano_lon);
fprintf('Orbital Period: %.2f minutes\n', orbital_period / 60);
fprintf('Orbital Altitude: %.2f km\n', altitude / 1000);
fprintf('Max horizon angle (theta_max): %.2f degrees\n', rad2deg(theta_max));
fprintf('Selected theta: %.2f degrees, phi: %.2f degrees\n\n', rad2deg(theta), rad2deg(phi));

visibility_at_eruption = dot(r_volcano_0, r_sat_0 - r_volcano_0);
fprintf('Visibility check at t=0: %.2e (should be > 0)\n', visibility_at_eruption);

zenith_at_volcano_0 = r_volcano_0 / norm(r_volcano_0);
volcano_to_sat_0 = r_sat_0 - r_volcano_0;
volcano_to_sat_0_unit = volcano_to_sat_0 / norm(volcano_to_sat_0);
cos_zenith_0 = dot(zenith_at_volcano_0, volcano_to_sat_0_unit);
zenith_distance_0 = rad2deg(acos(clamp_scalar(cos_zenith_0, -1, 1)));
fprintf('Zenith distance at t=0: %.2f degrees (should be < 90)\n\n', zenith_distance_0);

% Simulation window: half period before to half period after eruption
t_start = -orbital_period / 2;
t_end = orbital_period / 2;
t_eruption = 0;
dt = 1;
tspan = t_start:dt:t_end;

fprintf('Simulation from t = %.2f to %.2f seconds\n', t_start, t_end);
fprintf('Total simulation time: %.2f minutes\n\n', (t_end - t_start) / 60);

% Propagate the orbital state and free nadir attitude across the entire time window
fprintf('Propagating satellite trajectory with nadir hold before eruption...\n');
ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
ode_fun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, z_body_eci, params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false);

tspan_fwd = 0:dt:t_end;
[~, X_fwd] = ode45(ode_fun, tspan_fwd, X0, ode_opts);

tspan_bwd = 0:-dt:t_start;
[~, X_bwd] = ode45(ode_fun, tspan_bwd, X0, ode_opts);

X_orbit_history = [flipud(X_bwd(2:end, :)); X_fwd];
tspan = [fliplr(tspan_bwd(2:end)), tspan_fwd];
n_steps = size(X_orbit_history, 1);
idx_eruption = find(abs(tspan - t_eruption) < 1e-12, 1, 'first');

% Pre-allocate output arrays
ra_volcano_history = zeros(n_steps, 1);
dec_volcano_history = zeros(n_steps, 1);
visibility_history = false(n_steps, 1);
zenith_distance_history = zeros(n_steps, 1);
R_track_hist = zeros(3, 3, n_steps);
q_track_history = zeros(n_steps, 4);
omega_track_body_history = zeros(n_steps, 3);
omega_track_eci_history = zeros(n_steps, 3);
euler_track_history_deg = zeros(n_steps, 3);

% Compute the volcano-tracking reference history driven by orbital geometry
for i = 1:n_steps
    t = tspan(i);
    r_sat = X_orbit_history(i, 1:3)';
    v_sat = X_orbit_history(i, 4:6)';

    [r_volcano_t, v_volcano_t] = volcano_state_at_time(t, r_volcano_0, params.omega_earth);
    apparent_volcano_vec = r_volcano_t - r_sat;
    apparent_volcano_vec = apparent_volcano_vec / norm(apparent_volcano_vec);

    dec_rad = asin(apparent_volcano_vec(3));
    ra_rad = atan2(apparent_volcano_vec(2), apparent_volcano_vec(1));
    if ra_rad < 0
        ra_rad = ra_rad + 2 * pi;
    end
    ra_volcano_history(i) = rad2deg(ra_rad);
    dec_volcano_history(i) = rad2deg(dec_rad);

    visibility_check = dot(r_volcano_t, r_sat - r_volcano_t);
    visibility_history(i) = visibility_check > 0;

    zenith_at_volcano = r_volcano_t / norm(r_volcano_t);
    volcano_to_sat = r_sat - r_volcano_t;
    volcano_to_sat_unit = volcano_to_sat / norm(volcano_to_sat);
    cos_zenith = dot(zenith_at_volcano, volcano_to_sat_unit);
    zenith_distance_history(i) = rad2deg(acos(clamp_scalar(cos_zenith, -1, 1)));

    R_track_hist(:, :, i) = build_tracking_frame(r_sat, v_sat, r_volcano_t, v_volcano_t);
    q_track_history(i, :) = dcm2quat(R_track_hist(:, :, i)')';
    euler_track_history_deg(i, :) = rotation_matrix_to_euler_deg(R_track_hist(:, :, i));
end

q_track_history = normalize_quaternion_rows(q_track_history);
q_track_history = align_quaternion_signs(q_track_history);
omega_track_body_history = compute_body_rate_history_from_dcm_history(R_track_hist, dt);
alpha_track_body_history = compute_angular_acceleration_history(omega_track_body_history, dt);
tau_track_total_history = zeros(n_steps, 3);

for i = 1:n_steps
    omega_body = omega_track_body_history(i, :)';
    alpha_body = alpha_track_body_history(i, :)';
    tau_track_total_history(i, :) = rigid_body_total_torque(omega_body, alpha_body, params.I_CB)';
    omega_track_eci_history(i, :) = (R_track_hist(:, :, i) * omega_body)';
end

% Find visibility events
visibility_changes = diff([0; visibility_history]);
idx_start_visible = find(visibility_changes == 1, 1, 'first');
idx_end_visible = find(visibility_changes == -1, 1, 'first');

if ~isempty(idx_start_visible)
    t_start_visible = tspan(idx_start_visible);
    fprintf('Satellite starts seeing volcano at t = %.2f s (%.2f min)\n', ...
        t_start_visible, t_start_visible / 60);
else
    t_start_visible = NaN;
    if visibility_history(1)
        fprintf('Volcano visible at simulation start\n');
    else
        fprintf('Volcano not visible at simulation start\n');
    end
end

fprintf('Volcano erupts at t = %.2f s (%.2f min)\n', t_eruption, t_eruption / 60);

if ~isempty(idx_end_visible)
    t_end_visible = tspan(idx_end_visible);
    fprintf('Volcano goes out of sight at t = %.2f s (%.2f min)\n', ...
        t_end_visible, t_end_visible / 60);
else
    if visibility_history(end)
        error('control_test4:NoVisibilityEnd', ...
            ['Volcano remains visible at simulation end, so the default match-time ', ...
             'rule is undefined in this window. Increase the simulation window.']);
    end
    error('control_test4:NoVisibilityEnd', ...
        'Volcano is not visible after eruption in the current window.');
end

if t_end_visible <= t_eruption
    error('control_test4:InvalidVisibilityWindow', ...
        'The volcano is not visible long enough after eruption to define a post-eruption match time.');
end

t_match_default = t_eruption + 0.75 * (t_end_visible - t_eruption);
if isempty(targetMatchTime)
    t_match_request = t_match_default;
    matchSourceText = 'default 3/4 visibility rule';
else
    t_match_request = double(targetMatchTime);
    matchSourceText = 'user override';
end

if t_match_request <= t_eruption || t_match_request >= t_end_visible
    error('control_test4:MatchTimeOutOfRange', ...
        'targetMatchTime must be strictly between eruption time and end-visible time.');
end

[~, idx_match] = min(abs(tspan - t_match_request));
t_match = tspan(idx_match);
if idx_match <= idx_eruption || idx_match >= idx_end_visible
    error('control_test4:SnappedMatchTimeOutOfRange', ...
        'The requested match time snaps outside the valid post-eruption visible interval.');
end

fprintf('Requested match time (%s): %.2f s\n', matchSourceText, t_match_request);
fprintf('Using sampled match time: %.2f s (%.2f min)\n', t_match, t_match / 60);
fprintf('Maneuver duration after eruption: %.2f s\n\n', t_match - t_eruption);

% Initial nadir state at eruption and target volcano-tracking state at match time
X_eruption = X_orbit_history(idx_eruption, :)';
q_initial = X_eruption(7:10);
q_initial = q_initial / norm(q_initial);
R_ECI_to_Body_initial = quat2dcm(q_initial');
omega_initial_body = R_ECI_to_Body_initial * X_eruption(11:13);

q_target = q_track_history(idx_match, :)';
q_target = q_target / norm(q_target);
if dot(q_initial, q_target) < 0
    q_target = -q_target;
end
omega_target_body = omega_track_body_history(idx_match, :)';

fprintf('Target tracking state at t_match:\n');
fprintf('  Target camera RA/Dec: [%.3f, %.3f] deg\n', ...
    ra_volcano_history(idx_match), dec_volcano_history(idx_match));
fprintf('  Target body-rate magnitude: %.3e rad/s\n\n', norm(omega_target_body));

% Solve the rotational minimax slew problem
n_opt_intervals = 20;
fprintf('Solving fast bang-bang rotational slew with %d output intervals...\n', n_opt_intervals);
fprintf('Using an analytic bang-bang guess with optional low-dimensional refinement.\n');

opt_result = solve_minimax_attitude_slew(q_initial, omega_initial_body, q_target, ...
    omega_target_body, params.I_CB, t_match - t_eruption, n_opt_intervals);

fprintf('Slew solver: %s\n', opt_result.solverUsed);
fprintf('Optimization exit flag: %d\n', opt_result.exitflag);
fprintf('Optimization message: %s\n', opt_result.message);
fprintf('Peak total torque gamma: %.3e N m\n', opt_result.gamma);
fprintf('Grid max |tau|: %.3e N m\n', opt_result.maxTorqueNodeNorm);
fprintf('Max equality residual: %.3e\n', opt_result.maxEqualityResidual);
fprintf('Max inequality violation: %.3e\n', opt_result.maxInequalityViolation);
fprintf('Continuous terminal attitude error: %.3e deg\n', opt_result.terminalAngleErrorDeg);
fprintf('Continuous terminal body-rate error: %.3e deg/s\n\n', opt_result.terminalRateErrorDegPerSec);

% Build the full-timeline required total torque profile
tau_total_history = zeros(n_steps, 3);
if isfield(opt_result, 'interpMethod')
    maneuverInterpMethod = opt_result.interpMethod;
else
    maneuverInterpMethod = 'pchip';
end
tau_maneuver_history = interp1(opt_result.timeNodes, opt_result.tauNodes, ...
    tspan(idx_eruption:idx_match) - t_eruption, maneuverInterpMethod);
tau_total_history(idx_eruption:idx_match, :) = tau_maneuver_history;

if idx_match < idx_end_visible
    tau_total_history(idx_match+1:idx_end_visible, :) = ...
        tau_track_total_history(idx_match+1:idx_end_visible, :);
end

% First replay the total torque history with zero aerodynamic disturbance to obtain
% the intended full-state trajectory.
t_control = tspan(idx_eruption:end);
tau_total_interp = @(t) interp1(tspan, tau_total_history, t, 'pchip', 'extrap')';
zero_torque_interp = @(t) zeros(3, 1); %#ok<NASGU>

ode_fun_total = @(t, X) sat_dynamics_openloop(t, X, params, tau_total_interp, zero_torque_interp);
[~, X_post_total] = ode45(ode_fun_total, t_control, X_eruption, ode_opts);

X_total_history = X_orbit_history;
X_total_history(idx_eruption:end, :) = X_post_total;
X_total_history(:, 7:10) = normalize_quaternion_rows(X_total_history(:, 7:10));
X_total_history(:, 7:10) = align_quaternion_signs(X_total_history(:, 7:10));

% Desired state: nadir before eruption, optimal maneuver during the slew, then
% exact volcano-tracking attitude/rate after the match time while visible.
X_desired_history = X_total_history;
X_desired_history(1:idx_eruption-1, :) = X_orbit_history(1:idx_eruption-1, :);
X_desired_history(idx_match:idx_end_visible, 7:10) = q_track_history(idx_match:idx_end_visible, :);
X_desired_history(idx_match:idx_end_visible, 11:13) = omega_track_eci_history(idx_match:idx_end_visible, :);
X_desired_history(:, 7:10) = normalize_quaternion_rows(X_desired_history(:, 7:10));
X_desired_history(:, 7:10) = align_quaternion_signs(X_desired_history(:, 7:10));

R_desired_history = body_to_eci_history_from_quaternion_rows(X_desired_history(:, 7:10));

if params.includeAerodynamics
    fprintf('Evaluating aerodynamic torque history every %.0f seconds...\n', params.aeroSampleStep);
    tau_aero_history = estimate_aero_torque_history(tspan, X_desired_history, R_desired_history, params);
else
    fprintf('Aerodynamic torque estimate disabled; using zero disturbance torque.\n');
    tau_aero_history = zeros(n_steps, 3);
end

tau_control_history = tau_total_history - tau_aero_history;

% Replay with actuator plus aerodynamic torques so the total applied torque still
% matches the optimized/desired total rigid-body torque.
tau_control_interp = @(t) interp1(tspan, tau_control_history, t, 'pchip', 'extrap')';
tau_aero_interp = @(t) interp1(tspan, tau_aero_history, t, 'pchip', 'extrap')';
ode_fun_actual = @(t, X) sat_dynamics_openloop(t, X, params, tau_control_interp, tau_aero_interp);
[~, X_post_actual] = ode45(ode_fun_actual, t_control, X_eruption, ode_opts);

X_actual_history = X_orbit_history;
X_actual_history(idx_eruption:end, :) = X_post_actual;
X_actual_history(:, 7:10) = normalize_quaternion_rows(X_actual_history(:, 7:10));
X_actual_history(:, 7:10) = align_quaternion_signs(X_actual_history(:, 7:10));

% Post-process actual and desired attitude histories
plotTime = tspan(:);
plotDesiredEulerDeg = vleo.util.quat_history_to_euler_deg(X_desired_history(:, 7:10));
plotActualEulerDeg = vleo.util.quat_history_to_euler_deg(X_actual_history(:, 7:10));
plotDesiredEulerDeg = vleo.util.unwrap_angle_history_deg(plotDesiredEulerDeg);
plotActualEulerDeg = vleo.util.unwrap_angle_history_deg(plotActualEulerDeg);

omega_body_desired_history = eci_rates_to_body_history(X_desired_history(:, 11:13), X_desired_history(:, 7:10));
omega_body_actual_history = eci_rates_to_body_history(X_actual_history(:, 11:13), X_actual_history(:, 7:10));

ra_camera_actual_history = zeros(n_steps, 1);
dec_camera_actual_history = zeros(n_steps, 1);
pointing_eci_actual_history = zeros(n_steps, 3);
pointing_error_deg_history = zeros(n_steps, 1);

for i = 1:n_steps
    obs_actual = vleo.analysis.state_to_observation(X_actual_history(i, :)', params);
    ra_camera_actual_history(i) = obs_actual.ra;
    dec_camera_actual_history(i) = obs_actual.dec;
    pointing_eci_actual_history(i, :) = obs_actual.pointing_eci';

    r_sat = X_actual_history(i, 1:3)';
    [r_volcano_t, ~] = volcano_state_at_time(tspan(i), r_volcano_0, params.omega_earth);
    apparent_volcano_vec = r_volcano_t - r_sat;
    apparent_volcano_vec = apparent_volcano_vec / norm(apparent_volcano_vec);
    pointing_error_deg_history(i) = rad2deg(acos(clamp_scalar(dot(obs_actual.pointing_eci, apparent_volcano_vec), -1, 1)));
end

tracking_indices = idx_match:idx_end_visible;
ra_error = ra_camera_actual_history(tracking_indices) - ra_volcano_history(tracking_indices);
ra_error = mod(ra_error + 180, 360) - 180;
dec_error = dec_camera_actual_history(tracking_indices) - dec_volcano_history(tracking_indices);
body_rate_tracking_error_deg_s = vecnorm(rad2deg(omega_body_actual_history(tracking_indices, :) - ...
    omega_track_body_history(tracking_indices, :)), 2, 2);

match_attitude_error_deg = rotation_error_angle_deg(X_actual_history(idx_match, 7:10)', q_target);
match_rate_error_deg_s = norm(rad2deg(omega_body_actual_history(idx_match, :)' - omega_target_body));

fprintf('Match-state verification at t = %.2f s:\n', t_match);
fprintf('  Attitude error: %.3e deg\n', match_attitude_error_deg);
fprintf('  Body-rate error: %.3e deg/s\n', match_rate_error_deg_s);
fprintf('Post-match max pointing error: %.3e deg\n', max(pointing_error_deg_history(tracking_indices)));
fprintf('Post-match max RA/Dec error: [%.3e, %.3e] deg\n', max(abs(ra_error)), max(abs(dec_error)));
fprintf('Post-match max body-rate tracking error: %.3e deg/s\n', max(body_rate_tracking_error_deg_s));
fprintf('Max |tau_total|   = %.3e N m\n', max(vecnorm(tau_total_history, 2, 2)));
fprintf('Max |tau_control| = %.3e N m\n', max(vecnorm(tau_control_history, 2, 2)));
fprintf('Max |tau_aero|    = %.3e N m\n\n', max(vecnorm(tau_aero_history, 2, 2)));

% Save animation data for Python rendering
simulations_dir = vleo.util.simulations_dir();

anim_data_file = fullfile(simulations_dir, sprintf('control_test4_anim_data_aero_%s_seed_%s.mat', ...
    on_off_text(params.includeAerodynamics), scenarioSeedText));
save(anim_data_file, 'tspan', 'X_actual_history', 'visibility_history', ...
    'pointing_eci_actual_history', 'r_volcano_0', 'params', 'tau_total_history', ...
    'tau_control_history', 'tau_aero_history', 't_match');
fprintf('Animation data saved to: %s\n', anim_data_file);
fprintf('Simulation outputs directory: %s\n', simulations_dir);
fprintf('Run tools/render_gif.py to generate the GIF from the latest animation data.\n');

torque_csv_path = save_torque_history_csv(tspan, visibility_history, ...
    tau_total_history, tau_control_history, tau_aero_history, params, scenarioSeedText);
fprintf('Torque history CSV saved to: %s\n', torque_csv_path);

attitude_csv_path = save_attitude_history_csv(plotTime, plotDesiredEulerDeg, plotActualEulerDeg, ...
    omega_body_desired_history, omega_body_actual_history, params, scenarioSeedText);
fprintf('Attitude history CSV saved to: %s\n', attitude_csv_path);

if generateVisuals
    animationIndices = 1:numel(tspan);
    if ~isnan(t_start_visible) && ~isnan(t_end_visible)
        animationIndices = idx_start_visible:idx_end_visible;
    end

    guiResults = struct();
    guiResults.t = tspan(animationIndices)';
    guiResults.rs = X_actual_history(animationIndices, 1:3);
    guiResults.vs = X_actual_history(animationIndices, 4:6);
    guiResults.betas = X_actual_history(animationIndices, 7:10);
    guiResults.omegas = X_actual_history(animationIndices, 11:13);
    guiResults.torques = tau_control_history(animationIndices, :);
    guiResults.aeroTorques = tau_aero_history(animationIndices, :);
    guiResults.totalTorques = tau_total_history(animationIndices, :);
    guiResults.eulerDeg = plotActualEulerDeg(animationIndices, :);

    fprintf('Generating GUI-compatible animated results plot...\n');
    vleo.viz.animate_results(struct('results', guiResults));

    fprintf('Generating attitude summary plots...\n');
    figure('Name', 'Euler Angle Tracking', 'WindowStyle', 'normal', 'Position', [120 120 1500 650]);
    axEuler = axes();
    l1 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 1), 'r-', 'LineWidth', 2); hold(axEuler, 'on');
    l2 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 1), 'r--', 'LineWidth', 2);
    l3 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 2), 'b-', 'LineWidth', 2);
    l4 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 2), 'b--', 'LineWidth', 2);
    l5 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 3), 'k-', 'LineWidth', 2);
    l6 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 3), 'k--', 'LineWidth', 2);
    grid(axEuler, 'on');
    xlabel(axEuler, 'Time (minutes)', 'FontSize', 16);
    ylabel(axEuler, 'Angle (deg)', 'FontSize', 16);
    title(axEuler, 'Desired and Actual Euler Angles', 'FontSize', 18);
    legend(axEuler, [l1, l2, l3, l4, l5, l6], ...
        {'Roll desired', 'Roll actual', 'Pitch desired', 'Pitch actual', 'Yaw desired', 'Yaw actual'}, ...
        'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
    mark_phase_lines(axEuler, t_start_visible, t_eruption, t_match, t_end_visible, 13);
    set_visible_domain(axEuler, t_start_visible, t_end_visible);

    figure('Name', 'Body Rate Tracking', 'WindowStyle', 'normal', 'Position', [140 140 1500 650]);
    axRates = axes();
    r1 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 1)), 'r-', 'LineWidth', 2); hold(axRates, 'on');
    r2 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 1)), 'r--', 'LineWidth', 2);
    r3 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 2)), 'b-', 'LineWidth', 2);
    r4 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 2)), 'b--', 'LineWidth', 2);
    r5 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 3)), 'k-', 'LineWidth', 2);
    r6 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 3)), 'k--', 'LineWidth', 2);
    grid(axRates, 'on');
    xlabel(axRates, 'Time (minutes)', 'FontSize', 16);
    ylabel(axRates, 'Body Rate (deg/s)', 'FontSize', 16);
    title(axRates, 'Desired and Actual Body Rates', 'FontSize', 18);
    legend(axRates, [r1, r2, r3, r4, r5, r6], ...
        {'wx desired', 'wx actual', 'wy desired', 'wy actual', 'wz desired', 'wz actual'}, ...
        'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
    mark_phase_lines(axRates, t_start_visible, t_eruption, t_match, t_end_visible, 13);
    set_visible_domain(axRates, t_start_visible, t_end_visible);

    fprintf('Generating torque breakdown plots...\n');
    figure('Name', sprintf('Tracking Torque Breakdown (%s)', on_off_text(params.includeAerodynamics)), ...
        'WindowStyle', 'normal', 'Position', [160 160 1450 900]);
    torqueLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'x', 'y', 'z'};

    for axisIdx = 1:3
        ax = nexttile(torqueLayout, axisIdx);
        plot(ax, tspan / 60, tau_total_history(:, axisIdx), 'k-', 'LineWidth', 2.2); hold(ax, 'on');
        plot(ax, tspan / 60, tau_control_history(:, axisIdx), 'b-', 'LineWidth', 1.8);
        plot(ax, tspan / 60, tau_aero_history(:, axisIdx), 'r--', 'LineWidth', 1.8);
        grid(ax, 'on');
        ylabel(ax, sprintf('\\tau_%s (N\\cdotm)', axisLabels{axisIdx}), 'FontSize', 16);
        title(ax, sprintf('Body-%s Torque Component', upper(axisLabels{axisIdx})), 'FontSize', 17);
        mark_phase_lines(ax, t_start_visible, t_eruption, t_match, t_end_visible, 13);
        set_visible_domain(ax, t_start_visible, t_end_visible);
        if axisIdx == 1
            legend(ax, {'Total required', 'Actuator required', 'Aerodynamic disturbance'}, ...
                'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
        end
    end
    xlabel(torqueLayout, 'Time (minutes)', 'FontSize', 16);
    title(torqueLayout, sprintf('Torque Breakdown with Aerodynamics %s', upper(on_off_text(params.includeAerodynamics))), ...
        'FontSize', 18);

    figure('Name', 'Post-Match Tracking Error', 'WindowStyle', 'normal', 'Position', [180 180 1250 600]);
    plot(tspan(tracking_indices) / 60, pointing_error_deg_history(tracking_indices), 'k-', 'LineWidth', 2); hold on;
    plot(tspan(tracking_indices) / 60, body_rate_tracking_error_deg_s, 'm--', 'LineWidth', 2);
    xlabel('Time (minutes)', 'FontSize', 18);
    ylabel('Error', 'FontSize', 18);
    title('Post-Match Tracking Errors', 'FontSize', 20);
    legend('Pointing error (deg)', 'Body-rate error (deg/s)', 'Location', 'best', 'FontSize', 14);
    grid on;
    xlim([tspan(tracking_indices(1)) / 60, tspan(tracking_indices(end)) / 60]);
else
    fprintf('Skipping GUI animation and figure generation (generateVisuals=false).\n');
end

fprintf('\nSimulation complete!\n');

function I_CB = cuboid_inertia_matrix(mass, bodyDims)
    I_CB = (mass / 12) * diag([
        bodyDims(2)^2 + bodyDims(3)^2, ...
        bodyDims(1)^2 + bodyDims(3)^2, ...
        bodyDims(1)^2 + bodyDims(2)^2]);
end

function [rVolcano, vVolcano] = volcano_state_at_time(t, rVolcano0, omegaEarth)
    cosWt = cos(omegaEarth * t);
    sinWt = sin(omegaEarth * t);
    rVolcano = [cosWt * rVolcano0(1) - sinWt * rVolcano0(2); ...
                sinWt * rVolcano0(1) + cosWt * rVolcano0(2); ...
                rVolcano0(3)];
    vVolcano = cross([0; 0; omegaEarth], rVolcano);
end

function R_B2E = build_tracking_frame(rSat, vSat, rVolcano, vVolcano)
    rRel = rVolcano - rSat;
    vRel = vVolcano - vSat;

    zBody = rRel / norm(rRel);
    dist = norm(rRel);
    zBodyDot = (vRel / dist) - zBody * (dot(rRel, vRel) / dist^2);

    yBody = cross(zBody, zBodyDot);
    if norm(yBody) < 1e-10
        yBody = cross(zBody, vSat);
        if norm(yBody) < 1e-10
            yBody = cross(zBody, rSat);
        end
    end
    yBody = yBody / norm(yBody);

    xBody = cross(yBody, zBody);
    xBody = xBody / norm(xBody);
    yBody = cross(zBody, xBody);
    yBody = yBody / norm(yBody);

    R_B2E = [xBody, yBody, zBody];
end

function omegaBodyHistory = compute_body_rate_history_from_dcm_history(RHist, dt)
    nSteps = size(RHist, 3);
    RDot = zeros(3, 3, nSteps);
    omegaBodyHistory = zeros(nSteps, 3);

    for i = 2:nSteps-1
        RDot(:, :, i) = (RHist(:, :, i+1) - RHist(:, :, i-1)) / (2 * dt);
    end
    RDot(:, :, 1) = (RHist(:, :, 2) - RHist(:, :, 1)) / dt;
    RDot(:, :, nSteps) = (RHist(:, :, nSteps) - RHist(:, :, nSteps-1)) / dt;

    for i = 1:nSteps
        R = RHist(:, :, i);
        skewW = R' * RDot(:, :, i);
        omegaBodyHistory(i, 1) = skewW(3, 2);
        omegaBodyHistory(i, 2) = skewW(1, 3);
        omegaBodyHistory(i, 3) = skewW(2, 1);
    end
end

function alphaBodyHistory = compute_angular_acceleration_history(omegaBodyHistory, dt)
    nSteps = size(omegaBodyHistory, 1);
    alphaBodyHistory = zeros(nSteps, 3);
    for i = 2:nSteps-1
        alphaBodyHistory(i, :) = (omegaBodyHistory(i+1, :) - omegaBodyHistory(i-1, :)) / (2 * dt);
    end
    alphaBodyHistory(1, :) = (omegaBodyHistory(2, :) - omegaBodyHistory(1, :)) / dt;
    alphaBodyHistory(end, :) = (omegaBodyHistory(end, :) - omegaBodyHistory(end-1, :)) / dt;
end

function tauBody = rigid_body_total_torque(omegaBody, alphaBody, inertiaBody)
    tauBody = inertiaBody * alphaBody + cross(omegaBody, inertiaBody * omegaBody);
end

function result = solve_minimax_attitude_slew(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, nIntervals)
    nNodes = nIntervals + 1;
    timeNodes = linspace(0, maneuverDuration, nNodes)';

    [tauPositiveGuess, ~, switchFractionGuess] = bang_bang_initial_guess( ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody, maneuverDuration);

    x0 = [tauPositiveGuess; switchFractionGuess];
    residualFun = @(x) fast_bang_bang_terminal_residual(x, qInitial, omegaInitialBody, ...
        qTarget, omegaTargetBody, inertiaBody, maneuverDuration);

    solverUsed = 'analytic bang-bang approximation';
    exitflag = 0;
    output = struct('message', 'Used analytic bang-bang approximation for fast slew generation.');
    xOpt = x0;

    try
        fsOpts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 400, ...
            'MaxFunctionEvaluations', 4000, 'FunctionTolerance', 1e-10, ...
            'StepTolerance', 1e-14, 'OptimalityTolerance', 1e-10);
        [xSolved, ~, fsFlag, fsOutput] = fsolve(residualFun, x0, fsOpts);
        xSolved(4:6) = min(max(xSolved(4:6), 0.01), 0.99);

        if fsFlag > 0
            xOpt = xSolved;
            solverUsed = sprintf('bang-bang with fsolve refinement (exitflag=%d)', fsFlag);
            exitflag = fsFlag;
            output = fsOutput;
        else
            residualOpt = norm(residualFun(xSolved));
            residualGuess = norm(residualFun(x0));
            if residualOpt < residualGuess
                xOpt = xSolved;
                solverUsed = sprintf('bang-bang with partial fsolve refinement (exitflag=%d)', fsFlag);
                exitflag = fsFlag;
                output = fsOutput;
            end
        end
    catch ME
        warning('control_test4:fsolveUnavailable', ...
            'fsolve unavailable (%s); using analytic guess only.', ME.message);
    end

    tauPositive = xOpt(1:3);
    switchFractions = xOpt(4:6);
    tauNodes = build_bang_bang_torque_nodes(timeNodes, tauPositive, -tauPositive, switchFractions);
    [qNodes, omegaNodes] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, timeNodes);
    ceqOpt = [terminal_quaternion_error_vector(qNodes(end, :)', qTarget); ...
        omegaNodes(end, :)' - omegaTargetBody];
    cOpt = vecnorm(tauNodes, 2, 2) - norm(tauPositive);
    gamma = norm(tauPositive);

    timeVerify = linspace(0, maneuverDuration, max(200, 8 * nIntervals + 1))';
    [qVerify, omegaVerifyBody] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, timeVerify);

    result = struct();
    result.timeNodes = timeNodes;
    result.qNodes = qNodes;
    result.omegaNodesBody = omegaNodes;
    result.tauNodes = tauNodes;
    result.interpMethod = 'previous';
    result.solverUsed = solverUsed;
    result.gamma = gamma;
    result.exitflag = exitflag;
    result.output = output;
    result.message = output.message;
    result.maxTorqueNodeNorm = max(vecnorm(tauNodes, 2, 2));
    result.maxEqualityResidual = max(abs(ceqOpt));
    result.maxInequalityViolation = max(max(cOpt), 0);
    result.timeVerify = timeVerify;
    result.qVerify = qVerify;
    result.omegaVerifyBody = omegaVerifyBody;
    result.terminalAngleErrorDeg = rotation_error_angle_deg(qVerify(end, :)', qTarget);
    result.terminalRateErrorDegPerSec = norm(rad2deg(omegaVerifyBody(end, :)' - omegaTargetBody));
end

function x = pack_bang_bang_decision(tauPositive, tauNegative, switchFractions, gamma)
    x = [tauPositive(:); tauNegative(:); switchFractions(:); gamma];
end

function [tauPositive, tauNegative, switchFractions, gamma] = unpack_bang_bang_decision(x)
    tauPositive = x(1:3);
    tauNegative = x(4:6);
    switchFractions = x(7:9);
    gamma = x(10);
end

function [c, ceq, qNodes, omegaNodes, tauNodes] = bang_bang_shooting_constraints(x, timeNodes, ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody)
    [tauPositive, tauNegative, switchFractions, gamma] = unpack_bang_bang_decision(x);
    tauNodes = build_bang_bang_torque_nodes(timeNodes, tauPositive, tauNegative, switchFractions);
    [qNodes, omegaNodes] = simulate_attitude_maneuver(qInitial, omegaInitialBody, ...
        timeNodes, tauNodes, inertiaBody, timeNodes, 'previous');

    attitudeErrorBody = terminal_quaternion_error_vector(qNodes(end, :)', qTarget);
    rateErrorBody = omegaNodes(end, :)' - omegaTargetBody;

    ceq = [attitudeErrorBody; rateErrorBody];
    c = sum(tauNodes.^2, 2) - gamma^2;
end

function tauNodes = build_bang_bang_torque_nodes(timeNodes, tauPositive, tauNegative, switchFractions)
    nNodes = numel(timeNodes);
    maneuverDuration = timeNodes(end);
    switchTimes = maneuverDuration * switchFractions(:)';
    tauNodes = zeros(nNodes, 3);

    for nodeIdx = 1:nNodes
        t = timeNodes(nodeIdx);
        for axisIdx = 1:3
            if t <= switchTimes(axisIdx)
                tauNodes(nodeIdx, axisIdx) = tauPositive(axisIdx);
            else
                tauNodes(nodeIdx, axisIdx) = tauNegative(axisIdx);
            end
        end
    end
end

function [tauPositiveGuess, tauNegativeGuess, switchFractionGuess] = bang_bang_initial_guess( ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody, maneuverDuration)
    relativeRotationVector = relative_rotation_vector_body(qInitial, qTarget);
    principalInertia = diag(inertiaBody);
    tauPositiveGuess = zeros(3, 1);
    tauNegativeGuess = zeros(3, 1);
    switchFractionGuess = 0.5 * ones(3, 1);

    for axisIdx = 1:3
        [alphaAmplitude, switchTime] = solve_bang_bang_1d(0, omegaInitialBody(axisIdx), ...
            relativeRotationVector(axisIdx), omegaTargetBody(axisIdx), maneuverDuration);
        tauPositiveGuess(axisIdx) = principalInertia(axisIdx) * alphaAmplitude;
        tauNegativeGuess(axisIdx) = -principalInertia(axisIdx) * alphaAmplitude;
        switchFractionGuess(axisIdx) = clamp_scalar(switchTime / maneuverDuration, 0.05, 0.95);
    end
end

function residual = fast_bang_bang_terminal_residual(x, qInitial, omegaInitialBody, ...
        qTarget, omegaTargetBody, inertiaBody, maneuverDuration)
    tauPositive = x(1:3);
    switchFractions = x(4:6);
    [qFinal, omegaFinal] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, maneuverDuration);
    residual = [terminal_quaternion_error_vector(qFinal(end, :)', qTarget); ...
        omegaFinal(end, :)' - omegaTargetBody];
end

function [qHistory, omegaHistory] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, ...
        tauPositive, switchFractions, inertiaBody, maneuverDuration, outputTimes)
    outputTimes = reshape(outputTimes, [], 1);
    outputTimes = sort(outputTimes);
    if isempty(outputTimes) || abs(outputTimes(1)) > 1e-12
        outputTimes = [0; outputTimes];
    end
    if abs(outputTimes(end) - maneuverDuration) > 1e-12
        outputTimes = [outputTimes; maneuverDuration];
    end

    switchTimes = sort(unique(maneuverDuration * switchFractions(:)));
    segmentBoundaries = sort(unique([0; switchTimes; outputTimes; maneuverDuration]));

    qHistory = zeros(numel(outputTimes), 4);
    omegaHistory = zeros(numel(outputTimes), 3);
    state = [qInitial / norm(qInitial); omegaInitialBody];
    outputIdx = 1;
    qHistory(outputIdx, :) = state(1:4)';
    omegaHistory(outputIdx, :) = state(5:7)';
    outputIdx = outputIdx + 1;

    for segIdx = 1:numel(segmentBoundaries) - 1
        t0 = segmentBoundaries(segIdx);
        t1 = segmentBoundaries(segIdx + 1);
        if t1 <= t0
            continue;
        end

        tauSegment = symmetric_bang_bang_torque(0.5 * (t0 + t1), tauPositive, switchFractions, maneuverDuration);
        dtSegment = t1 - t0;
        nSubsteps = max(1, ceil(dtSegment / 0.25));
        h = dtSegment / nSubsteps;

        for stepIdx = 1:nSubsteps
            state = rk4_attitude_step(state, h, tauSegment, inertiaBody);
            state(1:4) = state(1:4) / norm(state(1:4));
        end

        while outputIdx <= numel(outputTimes) && abs(outputTimes(outputIdx) - t1) < 1e-10
            qHistory(outputIdx, :) = state(1:4)';
            omegaHistory(outputIdx, :) = state(5:7)';
            outputIdx = outputIdx + 1;
        end
    end

    qHistory = normalize_quaternion_rows(qHistory);
    qHistory = align_quaternion_signs(qHistory);
end

function tauBody = symmetric_bang_bang_torque(t, tauPositive, switchFractions, maneuverDuration)
    switchTimes = maneuverDuration * switchFractions(:);
    tauBody = tauPositive(:);
    for axisIdx = 1:3
        if t > switchTimes(axisIdx)
            tauBody(axisIdx) = -tauBody(axisIdx);
        end
    end
end

function nextState = rk4_attitude_step(state, h, tauBody, inertiaBody)
    k1 = attitude_state_derivative(state, tauBody, inertiaBody);
    k2 = attitude_state_derivative(state + 0.5 * h * k1, tauBody, inertiaBody);
    k3 = attitude_state_derivative(state + 0.5 * h * k2, tauBody, inertiaBody);
    k4 = attitude_state_derivative(state + h * k3, tauBody, inertiaBody);
    nextState = state + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
end

function stateDot = attitude_state_derivative(state, tauBody, inertiaBody)
    qScalarFirst = state(1:4);
    qScalarFirst = qScalarFirst / norm(qScalarFirst);
    omegaBody = state(5:7);

    stateDot = zeros(7, 1);
    stateDot(1:4) = quaternion_kinematics(qScalarFirst, omegaBody);
    stateDot(5:7) = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody);
end

function [qGuess, omegaGuess, tauGuess, gammaGuess] = initial_attitude_guess( ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody, timeNodes)
    nNodes = numel(timeNodes);
    maneuverDuration = timeNodes(end);
    fineCount = max(201, 20 * (nNodes - 1) + 1);
    timeFine = linspace(0, maneuverDuration, fineCount)';
    qRefFine = zeros(fineCount, 4);
    omegaRefFine = zeros(fineCount, 3);

    for i = 1:fineCount
        u = timeFine(i) / maneuverDuration;
        s = quintic_time_scaling(u);
        qRefFine(i, :) = slerp_quaternion(qInitial, qTarget, s)';
        omegaRefFine(i, :) = ((1 - u) * omegaInitialBody + u * omegaTargetBody)';
    end

    qRefFine = normalize_quaternion_rows(qRefFine);
    qRefFine = align_quaternion_signs(qRefFine);
    alphaRefFine = compute_angular_acceleration_history(omegaRefFine, timeFine(2) - timeFine(1));

    gainPairs = [1e-2, 5e-2; 5e-2, 1e-1; 1e-1, 2e-1; 2e-1, 4e-1];
    bestScore = inf;
    bestGuess = struct();

    for gainIdx = 1:size(gainPairs, 1)
        kp = gainPairs(gainIdx, 1);
        kd = gainPairs(gainIdx, 2);
        [qFine, omegaFine, tauFine, success] = simulate_guided_initial_guess( ...
            qInitial, omegaInitialBody, inertiaBody, timeFine, qRefFine, omegaRefFine, alphaRefFine, kp, kd);
        if ~success
            continue;
        end

        terminalAngleError = rotation_error_angle_deg(qFine(end, :)', qTarget);
        terminalRateError = norm(rad2deg(omegaFine(end, :)' - omegaTargetBody));
        score = terminalAngleError + terminalRateError;
        if score < bestScore
            bestScore = score;
            bestGuess.qFine = qFine;
            bestGuess.omegaFine = omegaFine;
            bestGuess.tauFine = tauFine;
        end
    end

    if isempty(fieldnames(bestGuess))
        qGuess = zeros(nNodes, 4);
        omegaGuess = zeros(nNodes, 3);
        for i = 1:nNodes
            u = timeNodes(i) / maneuverDuration;
            qGuess(i, :) = slerp_quaternion(qInitial, qTarget, u)';
            omegaGuess(i, :) = ((1 - u) * omegaInitialBody + u * omegaTargetBody)';
        end
        alphaGuess = compute_angular_acceleration_history(omegaGuess, timeNodes(2) - timeNodes(1));
        tauGuess = zeros(nNodes, 3);
        for i = 1:nNodes
            tauGuess(i, :) = rigid_body_total_torque(omegaGuess(i, :)', alphaGuess(i, :)', inertiaBody)';
        end
    else
        qGuess = interp1(timeFine, bestGuess.qFine, timeNodes, 'linear');
        omegaGuess = interp1(timeFine, bestGuess.omegaFine, timeNodes, 'pchip');
        tauGuess = interp1(timeFine, bestGuess.tauFine, timeNodes, 'pchip');
    end

    qGuess = normalize_quaternion_rows(qGuess);
    qGuess = align_quaternion_signs(qGuess);
    gammaGuess = max(vecnorm(tauGuess, 2, 2));
    gammaGuess = max(1.1 * gammaGuess, 1e-8);
end

function [qFine, omegaFine, tauFine, success] = simulate_guided_initial_guess( ...
        qInitial, omegaInitialBody, inertiaBody, timeFine, qRefFine, omegaRefFine, alphaRefFine, kp, kd)
    success = false;
    qFine = [];
    omegaFine = [];
    tauFine = [];

    y0 = [qInitial; omegaInitialBody];
    odeOpts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);

    try
        odeFun = @(t, y) guided_attitude_dynamics(t, y, inertiaBody, timeFine, qRefFine, omegaRefFine, alphaRefFine, kp, kd);
        [~, yFine] = ode45(odeFun, timeFine, y0, odeOpts);
        qFine = normalize_quaternion_rows(yFine(:, 1:4));
        qFine = align_quaternion_signs(qFine);
        omegaFine = yFine(:, 5:7);
        tauFine = zeros(size(omegaFine));
        for i = 1:numel(timeFine)
            qRef = interpolate_quaternion_rows(timeFine, qRefFine, timeFine(i))';
            omegaRef = interp1(timeFine, omegaRefFine, timeFine(i), 'pchip')';
            alphaRef = interp1(timeFine, alphaRefFine, timeFine(i), 'pchip')';
            tauFine(i, :) = guided_guess_torque(qFine(i, :)', omegaFine(i, :)', qRef, omegaRef, alphaRef, inertiaBody, kp, kd)';
        end
        success = all(isfinite(qFine), 'all') && all(isfinite(omegaFine), 'all') && all(isfinite(tauFine), 'all');
    catch
        success = false;
    end
end

function yDot = guided_attitude_dynamics(t, y, inertiaBody, timeFine, qRefFine, omegaRefFine, alphaRefFine, kp, kd)
    qCurrent = y(1:4);
    qCurrent = qCurrent / norm(qCurrent);
    omegaBody = y(5:7);

    qRef = interpolate_quaternion_rows(timeFine, qRefFine, t)';
    omegaRef = interp1(timeFine, omegaRefFine, t, 'pchip', 'extrap')';
    alphaRef = interp1(timeFine, alphaRefFine, t, 'pchip', 'extrap')';
    tauBody = guided_guess_torque(qCurrent, omegaBody, qRef, omegaRef, alphaRef, inertiaBody, kp, kd);

    yDot = zeros(7, 1);
    yDot(1:4) = quaternion_kinematics(qCurrent, omegaBody);
    yDot(5:7) = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody);
end

function tauBody = guided_guess_torque(qCurrent, omegaBody, qRef, omegaRef, alphaRef, inertiaBody, kp, kd)
    tauFeedforward = rigid_body_total_torque(omegaRef, alphaRef, inertiaBody);
    attitudeErrorBody = full_attitude_error_body(qCurrent, qRef);
    rateErrorBody = omegaRef - omegaBody;
    tauBody = tauFeedforward + kp * attitudeErrorBody + kd * rateErrorBody;
end

function qInterp = interpolate_quaternion_rows(timeGrid, qRows, queryTime)
    qInterp = interp1(timeGrid, qRows, queryTime, 'linear', 'extrap');
    qInterp = qInterp(:) / norm(qInterp);
end

function attitudeErrorBody = full_attitude_error_body(qCurrent, qTarget)
    R_B2E_Current = quat2dcm(qCurrent')';
    R_B2E_Target = quat2dcm(qTarget')';
    errorInertial = 0.5 * ( ...
        cross(R_B2E_Current(:, 1), R_B2E_Target(:, 1)) + ...
        cross(R_B2E_Current(:, 2), R_B2E_Target(:, 2)) + ...
        cross(R_B2E_Current(:, 3), R_B2E_Target(:, 3)));
    attitudeErrorBody = R_B2E_Current' * errorInertial;
end

function qErrorVector = terminal_quaternion_error_vector(qCurrent, qTarget)
    RCurrent = quat2dcm(qCurrent');
    RTarget = quat2dcm(qTarget');
    RErr = RTarget * RCurrent';
    qErr = dcm2quat(RErr)';
    if qErr(1) < 0
        qErr = -qErr;
    end
    qErrorVector = qErr(2:4);
end

function rotationVector = relative_rotation_vector_body(qInitial, qTarget)
    RInitial = quat2dcm(qInitial');
    RTarget = quat2dcm(qTarget');
    RRel = RTarget * RInitial';
    qRel = dcm2quat(RRel)';
    if qRel(1) < 0
        qRel = -qRel;
    end

    vectorNorm = norm(qRel(2:4));
    if vectorNorm < 1e-12
        rotationVector = [0; 0; 0];
        return;
    end

    angle = 2 * atan2(vectorNorm, qRel(1));
    axis = qRel(2:4) / vectorNorm;
    rotationVector = axis * angle;
end

function [amplitude, switchTime] = solve_bang_bang_1d(x0, v0, xf, vf, tau)
    if abs(vf - v0) < 1e-12
        switchTime = tau / 2;
        amplitude = 4 * (xf - x0 - v0 * tau) / tau^2;
        return;
    end

    candidateGuesses = [0.3, 0.5, 0.7, 0.25, 0.75] * tau;
    amplitude = 0;
    switchTime = tau / 2;

    for guess = candidateGuesses
        try
            root = fzero(@(ts) bang_bang_switch_residual(ts, x0, v0, xf, vf, tau), guess);
            if root > 0 && root < tau
                denominator = 2 * root - tau;
                if abs(denominator) > 1e-12
                    amplitude = (vf - v0) / denominator;
                    switchTime = root;
                    return;
                end
            end
        catch
        end
    end

    denominator = 2 * switchTime - tau;
    if abs(denominator) > 1e-12
        amplitude = (vf - v0) / denominator;
    end
end

function residual = bang_bang_switch_residual(ts, x0, v0, xf, vf, tau)
    denominator = 2 * ts - tau;
    if abs(denominator) < 1e-12
        residual = 1e10;
        return;
    end

    amplitude = (vf - v0) / denominator;
    shapeTerm = 2 * ts * tau - ts^2 - 0.5 * tau^2;
    xTau = x0 + v0 * tau + amplitude * shapeTerm;
    residual = xTau - xf;
end

function value = quintic_time_scaling(u)
    u = clamp_scalar(u, 0, 1);
    value = 10 * u^3 - 15 * u^4 + 6 * u^5;
end

function [c, ceq] = attitude_minimax_constraints(x, nNodes, dtNode, qInitial, ...
        omegaInitialBody, qTarget, omegaTargetBody, inertiaBody)
    [qNodes, omegaNodes, tauNodes, gamma] = unpack_attitude_decision(x, nNodes);

    ceqBoundary = [qNodes(1, :)' - qInitial; ...
                   omegaNodes(1, :)' - omegaInitialBody; ...
                   qNodes(end, :)' - qTarget; ...
                   omegaNodes(end, :)' - omegaTargetBody];

    ceqNorm = sum(qNodes.^2, 2) - 1;
    ceqDynamics = zeros((nNodes - 1) * 7, 1);

    for i = 1:nNodes-1
        q_i = qNodes(i, :)';
        q_ip1 = qNodes(i+1, :)';
        omega_i = omegaNodes(i, :)';
        omega_ip1 = omegaNodes(i+1, :)';
        tau_i = tauNodes(i, :)';
        tau_ip1 = tauNodes(i+1, :)';

        qdot_i = quaternion_kinematics(q_i, omega_i);
        qdot_ip1 = quaternion_kinematics(q_ip1, omega_ip1);
        omegadot_i = rigid_body_omega_dot(omega_i, tau_i, inertiaBody);
        omegadot_ip1 = rigid_body_omega_dot(omega_ip1, tau_ip1, inertiaBody);

        qDefect = q_ip1 - q_i - 0.5 * (qdot_i + qdot_ip1) * dtNode;
        omegaDefect = omega_ip1 - omega_i - 0.5 * (omegadot_i + omegadot_ip1) * dtNode;

        rowIdx = (i - 1) * 7 + (1:7);
        ceqDynamics(rowIdx) = [qDefect; omegaDefect];
    end

    ceq = [ceqBoundary; ceqNorm; ceqDynamics];
    c = sum(tauNodes.^2, 2) - gamma^2;
end

function x = pack_attitude_decision(qNodes, omegaNodes, tauNodes, gamma)
    x = [qNodes(:); omegaNodes(:); tauNodes(:); gamma];
end

function [qNodes, omegaNodes, tauNodes, gamma] = unpack_attitude_decision(x, nNodes)
    nQuat = 4 * nNodes;
    nOmega = 3 * nNodes;
    nTau = 3 * nNodes;

    qNodes = reshape(x(1:nQuat), nNodes, 4);
    omegaNodes = reshape(x(nQuat + (1:nOmega)), nNodes, 3);
    tauNodes = reshape(x(nQuat + nOmega + (1:nTau)), nNodes, 3);
    gamma = x(end);
end

function qDot = quaternion_kinematics(qScalarFirst, omegaBody)
    omegaMatrix = [0, -omegaBody(1), -omegaBody(2), -omegaBody(3); ...
        omegaBody(1), 0, omegaBody(3), -omegaBody(2); ...
        omegaBody(2), -omegaBody(3), 0, omegaBody(1); ...
        omegaBody(3), omegaBody(2), -omegaBody(1), 0];
    qDot = 0.5 * omegaMatrix * qScalarFirst;
end

function omegaDot = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody)
    omegaDot = inertiaBody \ (tauBody - cross(omegaBody, inertiaBody * omegaBody));
end

function [qHistory, omegaBodyHistory] = simulate_attitude_maneuver(qInitial, omegaInitialBody, ...
        timeNodes, tauNodes, inertiaBody, timeEval, interpMethod)
    if nargin < 7 || isempty(interpMethod)
        interpMethod = 'pchip';
    end

    tauInterp = @(t) interp1(timeNodes, tauNodes, t, interpMethod, 'extrap')';
    odeFun = @(t, y) attitude_dynamics_openloop(t, y, inertiaBody, tauInterp);
    y0 = [qInitial; omegaInitialBody];
    odeOpts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [~, yHistory] = ode45(odeFun, timeEval, y0, odeOpts);
    qHistory = normalize_quaternion_rows(yHistory(:, 1:4));
    qHistory = align_quaternion_signs(qHistory);
    omegaBodyHistory = yHistory(:, 5:7);
end

function yDot = attitude_dynamics_openloop(t, y, inertiaBody, tauInterp)
    qScalarFirst = y(1:4);
    qScalarFirst = qScalarFirst / norm(qScalarFirst);
    omegaBody = y(5:7);
    tauBody = tauInterp(t);

    yDot = zeros(7, 1);
    yDot(1:4) = quaternion_kinematics(qScalarFirst, omegaBody);
    yDot(5:7) = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody);
end

function qInterp = slerp_quaternion(q0, q1, s)
    q0 = q0 / norm(q0);
    q1 = q1 / norm(q1);

    dotProd = dot(q0, q1);
    if dotProd < 0
        q1 = -q1;
        dotProd = -dotProd;
    end

    if dotProd > 0.9995
        qInterp = q0 + s * (q1 - q0);
        qInterp = qInterp / norm(qInterp);
        return;
    end

    theta0 = acos(clamp_scalar(dotProd, -1, 1));
    sinTheta0 = sin(theta0);
    w0 = sin((1 - s) * theta0) / sinTheta0;
    w1 = sin(s * theta0) / sinTheta0;
    qInterp = w0 * q0 + w1 * q1;
    qInterp = qInterp / norm(qInterp);
end

function qRows = normalize_quaternion_rows(qRows)
    for rowIdx = 1:size(qRows, 1)
        qNorm = norm(qRows(rowIdx, :));
        if qNorm > 0
            qRows(rowIdx, :) = qRows(rowIdx, :) / qNorm;
        else
            qRows(rowIdx, :) = [1, 0, 0, 0];
        end
    end
end

function qRows = align_quaternion_signs(qRows)
    if isempty(qRows)
        return;
    end

    for rowIdx = 2:size(qRows, 1)
        if dot(qRows(rowIdx - 1, :), qRows(rowIdx, :)) < 0
            qRows(rowIdx, :) = -qRows(rowIdx, :);
        end
    end
end

function RHist = body_to_eci_history_from_quaternion_rows(qRows)
    nSteps = size(qRows, 1);
    RHist = zeros(3, 3, nSteps);
    for i = 1:nSteps
        RHist(:, :, i) = quat2dcm(qRows(i, :))';
    end
end

function omegaBodyHistory = eci_rates_to_body_history(omegaEciHistory, qRows)
    nSteps = size(omegaEciHistory, 1);
    omegaBodyHistory = zeros(nSteps, 3);
    for i = 1:nSteps
        R_ECI_to_Body = quat2dcm(qRows(i, :));
        omegaBodyHistory(i, :) = (R_ECI_to_Body * omegaEciHistory(i, :)')';
    end
end

function angleErrorDeg = rotation_error_angle_deg(qCurrent, qTarget)
    RCurrent = quat2dcm(qCurrent');
    RTarget = quat2dcm(qTarget');
    RErr = RCurrent * RTarget';
    traceVal = clamp_scalar((trace(RErr) - 1) / 2, -1, 1);
    angleErrorDeg = rad2deg(acos(traceVal));
end

function tauAeroHistory = estimate_aero_torque_history(tspan, X_history, R_B2E_hist, params)
    if numel(tspan) < 2
        tauAeroHistory = zeros(numel(tspan), 3);
        return;
    end

    dt = abs(tspan(2) - tspan(1));
    sampleStride = max(1, round(params.aeroSampleStep / dt));
    sampleIndices = unique([1:sampleStride:numel(tspan), numel(tspan)]);
    tauSamples = zeros(numel(sampleIndices), 3);

    for k = 1:numel(sampleIndices)
        idx = sampleIndices(k);
        R_ECI_to_Body = R_B2E_hist(:, :, idx)';
        tauSamples(k, :) = evaluate_aero_torque(tspan(idx), X_history(idx, :)', R_ECI_to_Body, params)';
    end

    if numel(sampleIndices) == 1
        tauAeroHistory = repmat(tauSamples, numel(tspan), 1);
        return;
    end

    tauAeroHistory = zeros(numel(tspan), 3);
    for axisIdx = 1:3
        tauAeroHistory(:, axisIdx) = interp1(tspan(sampleIndices), tauSamples(:, axisIdx), tspan, 'pchip');
    end
end

function tauAero = evaluate_aero_torque(t, X, R_ECI_to_Body, params)
    tauAero = [0; 0; 0];

    r_eci = X(1:3);
    v_eci = X(4:6);
    if ~(all(isfinite(r_eci)) && all(isfinite(v_eci)) && norm(r_eci) > 0)
        return;
    end

    v_body = R_ECI_to_Body * v_eci;
    if norm(v_body) > 1
        attitude.aoa = atan2d(v_body(3), v_body(1));
        attitude.aos = atan2d(v_body(2), v_body(1));
    else
        attitude.aoa = 0;
        attitude.aos = 0;
    end

    location.positionECI = r_eci;

    totalSeconds = params.aero.secondsOfDay + t;
    time.year = params.aero.year;
    time.dayOfYear = params.aero.dayOfYear + floor(totalSeconds / 86400);
    time.UTseconds = mod(totalSeconds, 86400);

    try
        aeroResults = vleo.aero.compute_aero_forces(params.aero.objFilePath, location, attitude, time, params.aero.options);
        tauCandidate = reshape(aeroResults.M_aero, [], 1);
        if isfield(aeroResults, 'M_solar')
            tauCandidate = tauCandidate + reshape(aeroResults.M_solar, [], 1);
        end
        if numel(tauCandidate) == 3 && all(isfinite(tauCandidate))
            tauAero = tauCandidate;
        end
    catch
        tauAero = [0; 0; 0];
    end
end

function csvPath = save_torque_history_csv(tspan, visibilityHistory, ...
        tauTotalHistory, tauControlHistory, tauAeroHistory, params, scenarioSeedText)
    simulationsDir = vleo.util.simulations_dir();

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    csvName = sprintf('control_test4_torque_history_aero_%s_seed_%s_%s.csv', ...
        on_off_text(params.includeAerodynamics), scenarioSeedText, timestamp);
    csvPath = fullfile(simulationsDir, csvName);

    torqueTable = table(tspan(:), tspan(:) / 60, logical(visibilityHistory(:)), ...
        tauTotalHistory(:, 1), tauTotalHistory(:, 2), tauTotalHistory(:, 3), ...
        tauControlHistory(:, 1), tauControlHistory(:, 2), tauControlHistory(:, 3), ...
        tauAeroHistory(:, 1), tauAeroHistory(:, 2), tauAeroHistory(:, 3), ...
        'VariableNames', {'time_s', 'time_min', 'visible', ...
        'tau_total_x', 'tau_total_y', 'tau_total_z', ...
        'tau_actuator_x', 'tau_actuator_y', 'tau_actuator_z', ...
        'tau_aero_x', 'tau_aero_y', 'tau_aero_z'});

    writetable(torqueTable, csvPath);
end

function csvPath = save_attitude_history_csv(plotTime, desiredEulerDeg, actualEulerDeg, ...
        omegaBodyDesiredHistory, omegaBodyActualHistory, params, scenarioSeedText)
    simulationsDir = vleo.util.simulations_dir();

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    csvName = sprintf('control_test4_attitude_history_aero_%s_seed_%s_%s.csv', ...
        on_off_text(params.includeAerodynamics), scenarioSeedText, timestamp);
    csvPath = fullfile(simulationsDir, csvName);

    attitudeTable = table(plotTime(:), plotTime(:) / 60, ...
        desiredEulerDeg(:, 1), desiredEulerDeg(:, 2), desiredEulerDeg(:, 3), ...
        actualEulerDeg(:, 1), actualEulerDeg(:, 2), actualEulerDeg(:, 3), ...
        rad2deg(omegaBodyDesiredHistory(:, 1)), rad2deg(omegaBodyDesiredHistory(:, 2)), rad2deg(omegaBodyDesiredHistory(:, 3)), ...
        rad2deg(omegaBodyActualHistory(:, 1)), rad2deg(omegaBodyActualHistory(:, 2)), rad2deg(omegaBodyActualHistory(:, 3)), ...
        'VariableNames', {'time_s', 'time_min', ...
        'roll_desired_deg', 'pitch_desired_deg', 'yaw_desired_deg', ...
        'roll_actual_deg', 'pitch_actual_deg', 'yaw_actual_deg', ...
        'wx_desired_deg_s', 'wy_desired_deg_s', 'wz_desired_deg_s', ...
        'wx_actual_deg_s', 'wy_actual_deg_s', 'wz_actual_deg_s'});

    writetable(attitudeTable, csvPath);
end

function eulerDeg = rotation_matrix_to_euler_deg(R_B2E)
    eulerDeg = vleo.util.quat_to_euler_deg(dcm2quat(R_B2E'));
end

function mark_phase_lines(ax, tStartVisible, tEruption, tMatch, tEndVisible, fontSize)
    hold(ax, 'on');
    if ~isnan(tStartVisible)
        marker = xline(ax, tStartVisible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', fontSize);
        hide_from_legend(marker);
    end
    marker = xline(ax, tEruption / 60, 'k--', 'LineWidth', 2, 'Label', 'Eruption', 'FontSize', fontSize);
    hide_from_legend(marker);
    marker = xline(ax, tMatch / 60, 'c-.', 'LineWidth', 2, 'Label', 'State Match', 'FontSize', fontSize);
    hide_from_legend(marker);
    if ~isnan(tEndVisible)
        marker = xline(ax, tEndVisible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', fontSize);
        hide_from_legend(marker);
    end
end

function set_visible_domain(ax, tStartVisible, tEndVisible)
    if ~isnan(tStartVisible) && ~isnan(tEndVisible)
        xlim(ax, [tStartVisible, tEndVisible] / 60);
    end
end

function hide_from_legend(graphicsHandle)
    if isprop(graphicsHandle, 'HandleVisibility')
        graphicsHandle.HandleVisibility = 'off';
    end

    try
        graphicsHandle.Annotation.LegendInformation.IconDisplayStyle = 'off';
    catch
    end
end

function Xd = sat_dynamics_openloop(t, X, params, tauControlInterp, tauAeroInterp)
    mu = params.mu;
    I_CB = params.I_CB;

    r_eci = X(1:3);
    v_eci = X(4:6);
    q_eci_to_body = X(7:10);
    omega_eci = X(11:13);

    qNorm = norm(q_eci_to_body);
    if qNorm < 1e-10 || ~isfinite(qNorm)
        Xd = zeros(13, 1);
        return;
    end
    q_eci_to_body = q_eci_to_body / qNorm;

    R_ECI_to_Body = quat2dcm(q_eci_to_body');
    R_Body_to_ECI = R_ECI_to_Body';

    omega_body = R_ECI_to_Body * omega_eci;

    Xd = zeros(13, 1);
    Xd(1:3) = v_eci;
    r_norm = norm(r_eci);
    Xd(4:6) = -mu * r_eci / r_norm^3;

    omegaMatrix = [0, -omega_body(1), -omega_body(2), -omega_body(3); ...
        omega_body(1), 0, omega_body(3), -omega_body(2); ...
        omega_body(2), -omega_body(3), 0, omega_body(1); ...
        omega_body(3), omega_body(2), -omega_body(1), 0];
    Xd(7:10) = 0.5 * omegaMatrix * q_eci_to_body;

    tau_body = tauControlInterp(t) + tauAeroInterp(t);
    alpha_body = I_CB \ (tau_body - cross(omega_body, I_CB * omega_body));
    Xd(11:13) = R_Body_to_ECI * alpha_body;
end

function textValue = on_off_text(flag)
    if flag
        textValue = 'on';
    else
        textValue = 'off';
    end
end

function clampedValue = clamp_scalar(value, lowerBound, upperBound)
    clampedValue = min(max(value, lowerBound), upperBound);
end
