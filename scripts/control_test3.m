% Control Test Earth - Volcano Observation Scenario
% Simulates a 6U CubeSat observing a volcano eruption with Earth rotation.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('includeAerodynamics', 'var')
    includeAerodynamics = true;
end

if ~exist('scenarioSeed', 'var')
    scenarioSeed = [];
end

clearvars -except includeAerodynamics scenarioSeed
clc;

if ~(isscalar(includeAerodynamics) && (islogical(includeAerodynamics) || isnumeric(includeAerodynamics)))
    error('control_test3:InvalidAeroToggle', ...
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
    error('control_test3:InvalidScenarioSeed', ...
        'scenarioSeed must be empty or a nonnegative integer scalar.');
end

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

fprintf('=== VOLCANO OBSERVATION SCENARIO ===\n');
fprintf('Spacecraft: 6U CubeSat (10 x 20 x 30 cm, 12 kg, uniform density)\n');
fprintf('Aerodynamic disturbance torque: %s\n', on_off_text(params.includeAerodynamics));
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
zenith_distance_0 = rad2deg(acos(cos_zenith_0));
fprintf('Zenith distance at t=0: %.2f degrees (should be < 90)\n\n', zenith_distance_0);

% Simulation window: half period before to half period after eruption
t_start = -orbital_period / 2;
t_end = orbital_period / 2;
t_eruption = 0;
dt = 1;
tspan = t_start:dt:t_end;

fprintf('Simulation from t = %.2f to %.2f seconds\n', t_start, t_end);
fprintf('Total simulation time: %.2f minutes\n\n', (t_end - t_start) / 60);

% Propagate the satellite state across the entire time window
fprintf('Propagating satellite trajectory...\n');
ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
ode_fun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, z_body_eci, params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false);

tspan_fwd = 0:dt:t_end;
[~, X_fwd] = ode45(ode_fun, tspan_fwd, X0, ode_opts);

tspan_bwd = 0:-dt:t_start;
[~, X_bwd] = ode45(ode_fun, tspan_bwd, X0, ode_opts);

X_history = [flipud(X_bwd(2:end, :)); X_fwd];
tspan = [fliplr(tspan_bwd(2:end)), tspan_fwd];
n_steps = size(X_history, 1);

% Pre-allocate output arrays
ra_camera_history = zeros(n_steps, 1);
dec_camera_history = zeros(n_steps, 1);
ra_volcano_history = zeros(n_steps, 1);
dec_volcano_history = zeros(n_steps, 1);
visibility_history = zeros(n_steps, 1);
zenith_distance_history = zeros(n_steps, 1);
pointing_eci_history = zeros(n_steps, 3);
R_B2E_hist = zeros(3, 3, n_steps);

% Post-process all time steps
for i = 1:n_steps
    t = tspan(i);

    obs = vleo.analysis.state_to_observation(X_history(i, :)', params);
    ra_camera_history(i) = obs.ra;
    dec_camera_history(i) = obs.dec;
    pointing_eci_history(i, :) = obs.pointing_eci';

    cos_wt = cos(params.omega_earth * t);
    sin_wt = sin(params.omega_earth * t);
    r_volcano_t = [cos_wt * r_volcano_0(1) - sin_wt * r_volcano_0(2);
                   sin_wt * r_volcano_0(1) + cos_wt * r_volcano_0(2);
                   r_volcano_0(3)];

    r_sat = X_history(i, 1:3)';
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
    zenith_distance_history(i) = rad2deg(acos(cos_zenith));

    omega_earth_vec = [0; 0; params.omega_earth];
    v_volcano_t = cross(omega_earth_vec, r_volcano_t);
    v_sat = X_history(i, 4:6)';
    r_rel = r_volcano_t - r_sat;
    v_rel = v_volcano_t - v_sat;

    z_b = r_rel / norm(r_rel);
    dist = norm(r_rel);
    z_b_dot = (v_rel / dist) - z_b * (dot(r_rel, v_rel) / dist^2);

    y_b = cross(z_b, z_b_dot);
    if norm(y_b) < 1e-10
        y_b = cross(z_b, v_sat);
        if norm(y_b) < 1e-10
            y_b = cross(z_b, r_sat);
        end
    end
    y_b = y_b / norm(y_b);

    x_b = cross(y_b, z_b);
    R_B2E_hist(:, :, i) = [x_b, y_b, z_b];
end

% Compute angular velocity from the DCM history
omega_body_history = zeros(n_steps, 3);
R_dot = zeros(3, 3, n_steps);

for i = 2:n_steps-1
    R_dot(:, :, i) = (R_B2E_hist(:, :, i+1) - R_B2E_hist(:, :, i-1)) / (2 * dt);
end
R_dot(:, :, 1) = (R_B2E_hist(:, :, 2) - R_B2E_hist(:, :, 1)) / dt;
R_dot(:, :, n_steps) = (R_B2E_hist(:, :, n_steps) - R_B2E_hist(:, :, n_steps-1)) / dt;

for i = 1:n_steps
    R = R_B2E_hist(:, :, i);
    skew_w = R' * R_dot(:, :, i);
    omega_body_history(i, 1) = skew_w(3, 2);
    omega_body_history(i, 2) = skew_w(1, 3);
    omega_body_history(i, 3) = skew_w(2, 1);
end

% Compute angular acceleration and torques
alpha_body_history = zeros(n_steps, 3);
for i = 2:n_steps-1
    alpha_body_history(i, :) = (omega_body_history(i+1, :) - omega_body_history(i-1, :)) / (2 * dt);
end
alpha_body_history(1, :) = (omega_body_history(2, :) - omega_body_history(1, :)) / dt;
alpha_body_history(end, :) = (omega_body_history(end, :) - omega_body_history(end-1, :)) / dt;

tau_total_history = zeros(n_steps, 3);
euler_desired_history_deg = zeros(n_steps, 3);
for i = 1:n_steps
    omega_body = omega_body_history(i, :)';
    alpha_body = alpha_body_history(i, :)';
    tau_total = params.I_CB * alpha_body + cross(omega_body, params.I_CB * omega_body);
    tau_total_history(i, :) = tau_total';
    euler_desired_history_deg(i, :) = rotation_matrix_to_euler_deg(R_B2E_hist(:, :, i));
end

if params.includeAerodynamics
    fprintf('Evaluating aerodynamic torque history every %.0f seconds...\n', params.aeroSampleStep);
    tau_aero_history = estimate_aero_torque_history(tspan, X_history, R_B2E_hist, params);
else
    fprintf('Aerodynamic torque disabled; using zero disturbance torque.\n');
    tau_aero_history = zeros(n_steps, 3);
end
tau_control_history = tau_total_history - tau_aero_history;

fprintf('Max |tau_total|   = %.3e N m\n', max(vecnorm(tau_total_history, 2, 2)));
fprintf('Max |tau_control| = %.3e N m\n', max(vecnorm(tau_control_history, 2, 2)));
fprintf('Max |tau_aero|    = %.3e N m\n\n', max(vecnorm(tau_aero_history, 2, 2)));

% Find visibility events
visibility_changes = diff([0; visibility_history]);
idx_start_visible = find(visibility_changes == 1);
idx_end_visible = find(visibility_changes == -1);

if ~isempty(idx_start_visible)
    t_start_visible = tspan(idx_start_visible(1));
    fprintf('Satellite starts seeing volcano at t = %.2f s (%.2f min)\n', ...
        t_start_visible, t_start_visible / 60);
else
    t_start_visible = NaN;
    if visibility_history(1) > 0
        fprintf('Volcano visible at simulation start\n');
    else
        fprintf('Volcano not visible at simulation start\n');
    end
end

fprintf('Volcano erupts at t = %.2f s (%.2f min)\n', t_eruption, t_eruption / 60);

if ~isempty(idx_end_visible)
    t_end_visible = tspan(idx_end_visible(1));
    fprintf('Volcano goes out of sight at t = %.2f s (%.2f min)\n\n', ...
        t_end_visible, t_end_visible / 60);
else
    t_end_visible = NaN;
    if visibility_history(end) > 0
        fprintf('Volcano still visible at simulation end\n\n');
    else
        fprintf('Volcano not visible at simulation end\n\n');
    end
end

% Save animation data for Python rendering
simulations_dir = vleo.util.simulations_dir();

anim_data_file = fullfile(simulations_dir, sprintf('control_test3_anim_data_aero_%s_seed_%s.mat', ...
    on_off_text(params.includeAerodynamics), scenarioSeedText));
save(anim_data_file, 'tspan', 'X_history', 'visibility_history', ...
    'pointing_eci_history', 'r_volcano_0', 'params', 'tau_total_history', ...
    'tau_control_history', 'tau_aero_history');
fprintf('Animation data saved to: %s\n', anim_data_file);
fprintf('Simulation outputs directory: %s\n', simulations_dir);
fprintf('Run tools/render_gif.py to generate the GIF from the latest animation data.\n');

torque_csv_path = save_torque_history_csv(tspan, visibility_history, ...
    tau_total_history, tau_control_history, tau_aero_history, params, scenarioSeedText);
fprintf('Torque history CSV saved to: %s\n', torque_csv_path);

plotTime = tspan;
plotDesiredEulerDeg = euler_desired_history_deg;
plotActualEulerDeg = vleo.util.quat_history_to_euler_deg(X_history(:, 7:10));
ra_error = [];
dec_error = [];
tspan_verif = [];

% Verification Simulation
if ~isnan(t_start_visible) && ~isnan(t_end_visible)
    fprintf('Running Verification Simulation...\n');
    idx_start = idx_start_visible(1);
    idx_end = idx_end_visible(1);

    t0_verif = tspan(idx_start);
    r_sat0 = X_history(idx_start, 1:3)';
    v_sat0 = X_history(idx_start, 4:6)';

    r_volcano_t0 = [cos(params.omega_earth * t0_verif) * r_volcano_0(1) - sin(params.omega_earth * t0_verif) * r_volcano_0(2);
                    sin(params.omega_earth * t0_verif) * r_volcano_0(1) + cos(params.omega_earth * t0_verif) * r_volcano_0(2);
                    r_volcano_0(3)];
    omega_earth_vec = [0; 0; params.omega_earth];
    v_volcano_t0 = cross(omega_earth_vec, r_volcano_t0);

    r_rel0 = r_volcano_t0 - r_sat0;
    v_rel0 = v_volcano_t0 - v_sat0;

    z_b0 = r_rel0 / norm(r_rel0);
    dist0 = norm(r_rel0);
    z_b_dot0 = (v_rel0 / dist0) - z_b0 * (dot(r_rel0, v_rel0) / dist0^2);

    y_b0 = cross(z_b0, z_b_dot0);
    if norm(y_b0) < 1e-10
        y_b0 = cross(z_b0, v_sat0);
        if norm(y_b0) < 1e-10
            y_b0 = cross(z_b0, r_sat0);
        end
    end
    y_b0 = y_b0 / norm(y_b0);

    x_b0 = cross(y_b0, z_b0);
    R_ECI_to_Body0 = [x_b0, y_b0, z_b0]';
    q0_verif = dcm2quat(R_ECI_to_Body0)';

    omega0_verif_body = omega_body_history(idx_start, :)';
    omega0_verif_eci = R_ECI_to_Body0' * omega0_verif_body;
    X0_verif = [r_sat0; v_sat0; q0_verif; omega0_verif_eci];

    tau_control_interp = @(t) interp1(tspan, tau_control_history, t, 'pchip', 'extrap')';
    tau_aero_interp = @(t) interp1(tspan, tau_aero_history, t, 'pchip', 'extrap')';
    ode_fun_verif = @(t, X) sat_dynamics_openloop(t, X, params, tau_control_interp, tau_aero_interp);

    tspan_verif = tspan(idx_start:idx_end);
    [~, X_verif] = ode45(ode_fun_verif, tspan_verif, X0_verif, ode_opts);

    plotTime = tspan_verif;
    plotDesiredEulerDeg = euler_desired_history_deg(idx_start:idx_end, :);
    plotActualEulerDeg = vleo.util.quat_history_to_euler_deg(X_verif(:, 7:10));

    ra_verif_history = zeros(length(tspan_verif), 1);
    dec_verif_history = zeros(length(tspan_verif), 1);
    for k = 1:length(tspan_verif)
        obs_verif = vleo.analysis.state_to_observation(X_verif(k, :)', params);
        ra_verif_history(k) = obs_verif.ra;
        dec_verif_history(k) = obs_verif.dec;
    end

    ra_target_history = ra_volcano_history(idx_start:idx_end);
    dec_target_history = dec_volcano_history(idx_start:idx_end);

    ra_error = ra_verif_history - ra_target_history;
    ra_error = mod(ra_error + 180, 360) - 180;
    dec_error = dec_verif_history - dec_target_history;

    fprintf('Max RA Error: %.2e deg\n', max(abs(ra_error)));
    fprintf('Max Dec Error: %.2e deg\n', max(abs(dec_error)));
end

plotDesiredEulerDeg = vleo.util.unwrap_angle_history_deg(plotDesiredEulerDeg);
plotActualEulerDeg = vleo.util.unwrap_angle_history_deg(plotActualEulerDeg);

euler_csv_path = save_euler_history_csv(plotTime, plotDesiredEulerDeg, plotActualEulerDeg, params, scenarioSeedText);
fprintf('Euler history CSV saved to: %s\n', euler_csv_path);

animationIndices = 1:numel(tspan);
if ~isnan(t_start_visible) && ~isnan(t_end_visible)
    animationIndices = idx_start_visible(1):idx_end_visible(1);
end

animationTime = tspan(animationIndices);
animationStateHistory = X_history(animationIndices, :);
animationControlHistory = tau_control_history(animationIndices, :);
animationAeroHistory = tau_aero_history(animationIndices, :);

if ~isempty(tspan_verif)
    animationTime = tspan_verif;
    animationStateHistory = X_verif;
end

guiResults = struct();
guiResults.t = animationTime(:);
guiResults.rs = animationStateHistory(:, 1:3);
guiResults.vs = animationStateHistory(:, 4:6);
guiResults.betas = animationStateHistory(:, 7:10);
guiResults.omegas = animationStateHistory(:, 11:13);
guiResults.torques = animationControlHistory;
guiResults.aeroTorques = animationAeroHistory;
guiResults.totalTorques = animationControlHistory + animationAeroHistory;
guiResults.eulerDeg = plotActualEulerDeg;

fprintf('Generating GUI-compatible animated results plot...\n');
vleo.viz.animate_results(struct('results', guiResults));

% Attitude summary plot
fprintf('Generating attitude summary plot...\n');

figure('Name', 'Euler Angle Tracking', 'WindowStyle', 'normal', 'Position', [120 120 1500 650]);
attitudeLayout = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axEuler = nexttile(attitudeLayout, 1);
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
    {'Roll \phi desired', 'Roll \phi actual', 'Pitch \theta desired', 'Pitch \theta actual', ...
    'Yaw \psi desired', 'Yaw \psi actual'}, 'Location', 'best', 'FontSize', 14, 'AutoUpdate', 'off');
mark_visibility_window(axEuler, t_start_visible, t_eruption, t_end_visible, 13);
set_visible_domain(axEuler, t_start_visible, t_end_visible);

% Torque breakdown plots
fprintf('Generating torque breakdown plots...\n');
figure('Name', sprintf('Tracking Torque Breakdown (%s)', on_off_text(params.includeAerodynamics)), ...
    'WindowStyle', 'normal', 'Position', [140 140 1450 900]);
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
    mark_visibility_window(ax, t_start_visible, t_eruption, t_end_visible, 13);
    set_visible_domain(ax, t_start_visible, t_end_visible);
    if axisIdx == 1
        legend(ax, {'Total required', 'Actuator required', 'Aerodynamic disturbance'}, ...
            'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
    end
end
xlabel(torqueLayout, 'Time (minutes)', 'FontSize', 16);
title(torqueLayout, sprintf('Torque Breakdown with Aerodynamics %s', upper(on_off_text(params.includeAerodynamics))), ...
    'FontSize', 18);

if ~isempty(ra_error)
    figure('Name', 'Verification Pointing Error', 'WindowStyle', 'normal', 'Position', [150 150 1200 600]);
    plot(tspan_verif / 60, ra_error, 'r-', 'LineWidth', 2); hold on;
    plot(tspan_verif / 60, dec_error, 'b--', 'LineWidth', 2);
    xlabel('Time (minutes)', 'FontSize', 20);
    ylabel('Pointing Error (degrees)', 'FontSize', 20);
    title(sprintf('Verification: Pointing Error with Aerodynamics %s', upper(on_off_text(params.includeAerodynamics))), ...
        'FontSize', 24);
    legend('RA Error', 'Dec Error', 'Location', 'best', 'FontSize', 16);
    set(gca, 'FontSize', 16);
    grid on;
    xlim([tspan_verif(1) / 60, tspan_verif(end) / 60]);
end

fprintf('\nSimulation complete!\n');

function I_CB = cuboid_inertia_matrix(mass, bodyDims)
    I_CB = (mass / 12) * diag([
        bodyDims(2)^2 + bodyDims(3)^2, ...
        bodyDims(1)^2 + bodyDims(3)^2, ...
        bodyDims(1)^2 + bodyDims(2)^2]);
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
    csvName = sprintf('control_test3_torque_history_aero_%s_seed_%s_%s.csv', ...
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

function csvPath = save_euler_history_csv(plotTime, desiredEulerDeg, actualEulerDeg, params, scenarioSeedText)
    simulationsDir = vleo.util.simulations_dir();

    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    csvName = sprintf('control_test3_euler_history_aero_%s_seed_%s_%s.csv', ...
        on_off_text(params.includeAerodynamics), scenarioSeedText, timestamp);
    csvPath = fullfile(simulationsDir, csvName);

    eulerTable = table(plotTime(:), plotTime(:) / 60, ...
        desiredEulerDeg(:, 1), desiredEulerDeg(:, 2), desiredEulerDeg(:, 3), ...
        actualEulerDeg(:, 1), actualEulerDeg(:, 2), actualEulerDeg(:, 3), ...
        'VariableNames', {'time_s', 'time_min', ...
        'roll_desired_deg', 'pitch_desired_deg', 'yaw_desired_deg', ...
        'roll_actual_deg', 'pitch_actual_deg', 'yaw_actual_deg'});

    writetable(eulerTable, csvPath);
end

function eulerDeg = rotation_matrix_to_euler_deg(R_B2E)
    eulerDeg = vleo.util.quat_to_euler_deg(dcm2quat(R_B2E'));
end

function mark_visibility_window(ax, tStartVisible, tEruption, tEndVisible, fontSize)
    hold(ax, 'on');
    if ~isnan(tStartVisible)
        marker = xline(ax, tStartVisible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', fontSize);
        hide_from_legend(marker);
    end
    marker = xline(ax, tEruption / 60, 'k--', 'LineWidth', 2, 'Label', 'Eruption (t=0)', 'FontSize', fontSize);
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

    q_eci_to_body = q_eci_to_body / norm(q_eci_to_body);

    R_ECI_to_Body = quat2dcm(q_eci_to_body');
    R_Body_to_ECI = R_ECI_to_Body';

    omega_body = R_ECI_to_Body * omega_eci;

    Xd = zeros(13, 1);
    Xd(1:3) = v_eci;
    r_norm = norm(r_eci);
    Xd(4:6) = -mu * r_eci / r_norm^3;

    Omega_matrix = [0, -omega_body(1), -omega_body(2), -omega_body(3);
        omega_body(1), 0, omega_body(3), -omega_body(2);
        omega_body(2), -omega_body(3), 0, omega_body(1);
        omega_body(3), omega_body(2), -omega_body(1), 0];
    Xd(7:10) = 0.5 * Omega_matrix * q_eci_to_body;

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
