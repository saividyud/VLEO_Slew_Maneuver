% Control Test Earth - Volcano Observation Scenario
% Simulates satellite observation of a volcano eruption with Earth rotation
% Volcano erupts at t=0 and is guaranteed to be visible at that time

clear; close all; clc;
rng('shuffle');

% Parameters
params.mu = 3.986004418e14; % Earth gravitational parameter [m^3/s^2]
params.R_e = 6378137; % Earth radius [m]
params.mass = 83.6; % Satellite mass [kg]
params.radius = 0.58/2; % Satellite radius [m]
params.J2 = 1.08263e-3; % J2 coefficient
params.I_CB = 2/5 * params.mass * params.radius^2 * eye(3); % Moment of inertia [kg*m^2]
params.Kp_att = 10; % Proportional control gain [N*m]
params.omega_earth = 7.2921159e-5; % Earth rotation rate [rad/s]

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
% Use Aerospace Toolbox dcm2quat (which takes DCM mapping ECI->Body and outputs scalar-first [w x y z])
q0_aerospace = dcm2quat(R_Body_to_ECI');
q0 = [q0_aerospace(2:4), q0_aerospace(1)]'; % Convert to scalar-last [x y z w]

omega_orbit_mag = v_circ / r_orbit;
omega_eci = omega_orbit_mag * orbit_normal_eci;

X0 = [r_eci; v_eci; q0; omega_eci];

% Calculate orbital period
orbital_period = 2 * pi * sqrt(r_orbit^3 / params.mu); % seconds

% === Place volcano so it's visible at t=0 (eruption time) ===
% Project satellite position onto Earth's surface
r_sat_0 = X0(1:3);
r_projected = r_sat_0 / norm(r_sat_0) * params.R_e;

% Maximum angle from nadir to horizon (as seen from Earth's center)
theta_max = acos(params.R_e / norm(r_sat_0));

% Random placement within visible cone
theta = rand() * theta_max;  % Angle from nadir
phi = rand() * 2 * pi;        % Azimuthal angle

% First rotation: rotate r_projected by theta around v_eci (axis orthogonal to r_sat_0)
axis1 = v_eci / norm(v_eci);  % Use velocity as rotation axis

% Rodrigues rotation formula: v_rot = v*cos(theta) + (k × v)*sin(theta) + k*(k·v)*(1-cos(theta))
% where k is the unit rotation axis
r_temp = r_projected * cos(theta) + ...
         cross(axis1, r_projected) * sin(theta) + ...
         axis1 * dot(axis1, r_projected) * (1 - cos(theta));

% Second rotation: rotate r_temp by phi around r_projected (or r_eci)
axis2 = r_projected / norm(r_projected);  % Radial direction as rotation axis

r_volcano_0 = r_temp * cos(phi) + ...
              cross(axis2, r_temp) * sin(phi) + ...
              axis2 * dot(axis2, r_temp) * (1 - cos(phi));

% Normalize to Earth radius (should already be close, but ensure precision)
r_volcano_0 = r_volcano_0 / norm(r_volcano_0) * params.R_e;

% Convert to lat/lon for display
volcano_direction = r_volcano_0 / norm(r_volcano_0);
volcano_lat = rad2deg(asin(volcano_direction(3)));
volcano_lon = rad2deg(atan2(volcano_direction(2), volcano_direction(1)));

fprintf('=== VOLCANO OBSERVATION SCENARIO ===\n');
fprintf('Volcano erupts at t = 0\n');
fprintf('Volcano Location (at t=0): Lat=%.2f deg, Lon=%.2f deg\n', volcano_lat, volcano_lon);
fprintf('Orbital Period: %.2f minutes\n', orbital_period / 60);
fprintf('Orbital Altitude: %.2f km\n', altitude / 1000);
fprintf('Max horizon angle (theta_max): %.2f degrees\n', rad2deg(theta_max));
fprintf('Selected theta: %.2f degrees, phi: %.2f degrees\n\n', rad2deg(theta), rad2deg(phi));

% Verify visibility at t=0
visibility_at_eruption = dot(r_volcano_0, r_sat_0 - r_volcano_0);
fprintf('Visibility check at t=0: %.2e (should be > 0)\n', visibility_at_eruption);

% Calculate zenith distance at t=0
zenith_at_volcano_0 = r_volcano_0 / norm(r_volcano_0);
volcano_to_sat_0 = r_sat_0 - r_volcano_0;
volcano_to_sat_0_unit = volcano_to_sat_0 / norm(volcano_to_sat_0);
cos_zenith_0 = dot(zenith_at_volcano_0, volcano_to_sat_0_unit);
zenith_distance_0 = rad2deg(acos(cos_zenith_0));
fprintf('Zenith distance at t=0: %.2f degrees (should be < 90)\n\n', zenith_distance_0);

% Simulation window: half period before to half period after eruption
t_start = -orbital_period / 2;
t_end = orbital_period / 2;
t_eruption = 0; % Eruption at t=0
dt = 1; % Decrease time step to 1 second for better numerical differentiation and torque tracking
tspan = t_start:dt:t_end;
n_steps = length(tspan);

fprintf('Simulation from t = %.2f to %.2f seconds\n', t_start, t_end);
fprintf('Total simulation time: %.2f minutes\n\n', (t_end - t_start) / 60);

% Propagate the satellite state across the entire time window
fprintf('Propagating satellite trajectory...\n');

% Two ODE45 calls instead of one per time step
ode_opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
ode_fun = @(t,X) Sat_template_control(t, X, z_body_eci, params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false);

% Forward propagation: t=0 to t_end
tspan_fwd = 0:dt:t_end;
[~, X_fwd] = ode45(ode_fun, tspan_fwd, X0, ode_opts);

% Backward propagation: t=0 to t_start
tspan_bwd = 0:-dt:t_start;
[~, X_bwd] = ode45(ode_fun, tspan_bwd, X0, ode_opts);

% Combine: flip backward (exclude t=0 duplicate), then forward
X_history = [flipud(X_bwd(2:end,:)); X_fwd];

% The concatenated time arrays might have a slightly different size due to length(tspan) 
% vs length(X_history). Make sure they match exactly.
n_steps = size(X_history, 1);
tspan = [fliplr(tspan_bwd(2:end)), tspan_fwd];

% Pre-allocate output arrays
ra_camera_history = zeros(n_steps, 1);
dec_camera_history = zeros(n_steps, 1);
ra_volcano_history = zeros(n_steps, 1);
dec_volcano_history = zeros(n_steps, 1);
visibility_history = zeros(n_steps, 1);
zenith_distance_history = zeros(n_steps, 1);
pointing_eci_history = zeros(n_steps, 3);
omega_req_eci_history = zeros(n_steps, 3);
tau_body_history = zeros(n_steps, 3);

% Post-process all time steps (no ODE integration needed)
for i = 1:n_steps
    t = tspan(i);

    % Get camera pointing direction (RA/Dec)
    obs = state_to_observation(X_history(i, :)', params);
    ra_camera_history(i) = obs.ra;
    dec_camera_history(i) = obs.dec;
    pointing_eci_history(i,:) = obs.pointing_eci';

    % Volcano position at time t (rotates with Earth)
    cos_wt = cos(params.omega_earth * t);
    sin_wt = sin(params.omega_earth * t);
    r_volcano_t = [cos_wt * r_volcano_0(1) - sin_wt * r_volcano_0(2);
                   sin_wt * r_volcano_0(1) + cos_wt * r_volcano_0(2);
                   r_volcano_0(3)];

    % Apparent volcano direction from satellite
    r_sat = X_history(i, 1:3)';
    apparent_volcano_vec = r_volcano_t - r_sat;
    apparent_volcano_vec = apparent_volcano_vec / norm(apparent_volcano_vec);

    % Convert to RA/Dec
    dec_rad = asin(apparent_volcano_vec(3));
    ra_rad = atan2(apparent_volcano_vec(2), apparent_volcano_vec(1));
    if ra_rad < 0
        ra_rad = ra_rad + 2*pi;
    end
    ra_volcano_history(i) = rad2deg(ra_rad);
    dec_volcano_history(i) = rad2deg(dec_rad);

    % Visibility check: r_volcano · (r_sat - r_volcano) > 0
    visibility_check = dot(r_volcano_t, r_sat - r_volcano_t);
    visibility_history(i) = visibility_check > 0;

    % Zenith distance
    zenith_at_volcano = r_volcano_t / norm(r_volcano_t);
    volcano_to_sat = r_sat - r_volcano_t;
    volcano_to_sat_unit = volcano_to_sat / norm(volcano_to_sat);
    cos_zenith = dot(zenith_at_volcano, volcano_to_sat_unit);
    zenith_distance_history(i) = rad2deg(acos(cos_zenith));

    % Required angular velocity for tracking
    omega_earth_vec = [0; 0; params.omega_earth];
    v_volcano_t = cross(omega_earth_vec, r_volcano_t);
    v_sat = X_history(i, 4:6)';
    r_rel = r_volcano_t - r_sat;
    v_rel = v_volcano_t - v_sat;
    
    % Use kinematic approach with DCM to guarantee correct body frame omega
    % Define desired Body Frame
    % Z points to volcano
    z_b = r_rel / norm(r_rel); 
    
    % Minimum-roll tracking frame: Y is perpendicular to the instantaneous 
    % plane of motion of the pointing vector (Z cross Z_dot)
    dist = norm(r_rel);
    z_b_dot = (v_rel / dist) - z_b * (dot(r_rel, v_rel) / dist^2);
    
    y_b = cross(z_b, z_b_dot);
    if norm(y_b) < 1e-10
        % Fallback if looking directly along the velocity vector
        y_b = cross(z_b, v_sat); 
        if norm(y_b) < 1e-10
            y_b = cross(z_b, r_sat);
        end
    end
    y_b = y_b / norm(y_b);
    
    % X completes the triad
    x_b = cross(y_b, z_b);
    
    % This is the rotation matrix from Body to ECI
    R_B2E_hist(:,:,i) = [x_b, y_b, z_b];
end

% Compute angular velocity by numerically differentiating the DCM
% R_dot = R * skew(omega_body) -> skew(omega_body) = R' * R_dot
omega_body_history = zeros(n_steps, 3);
R_dot = zeros(3, 3, n_steps);

for i = 2:n_steps-1
    R_dot(:,:,i) = (R_B2E_hist(:,:,i+1) - R_B2E_hist(:,:,i-1)) / (2*dt);
end
R_dot(:,:,1) = (R_B2E_hist(:,:,2) - R_B2E_hist(:,:,1)) / dt;
R_dot(:,:,n_steps) = (R_B2E_hist(:,:,n_steps) - R_B2E_hist(:,:,n_steps-1)) / dt;

for i = 1:n_steps
    R = R_B2E_hist(:,:,i);
    R_d = R_dot(:,:,i);
    skew_w = R' * R_d;
    
    % Extract omega from skew-symmetric matrix:
    % [  0   -w3   w2]
    % [ w3    0   -w1]
    % [-w2   w1    0 ]
    omega_body_history(i, 1) = skew_w(3, 2);
    omega_body_history(i, 2) = skew_w(1, 3);
    omega_body_history(i, 3) = skew_w(2, 1);
end

% Calculate angular acceleration in body frame
alpha_body_history = zeros(n_steps, 3);
for i = 2:n_steps-1
    alpha_body_history(i, :) = (omega_body_history(i+1, :) - omega_body_history(i-1, :)) / (2 * dt);
end
alpha_body_history(1, :) = (omega_body_history(2, :) - omega_body_history(1, :)) / dt;
alpha_body_history(end, :) = (omega_body_history(end, :) - omega_body_history(end-1, :)) / dt;

% Calculate torque in body frame
for i = 1:n_steps
    omega_body = omega_body_history(i, :)';
    alpha_body = alpha_body_history(i, :)';
    
    tau_body = params.I_CB * alpha_body + cross(omega_body, params.I_CB * omega_body);
    tau_body_history(i, :) = tau_body';
end

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

% === Save animation data for Python rendering ===
anim_data_file = fullfile(fileparts(mfilename('fullpath')), 'anim_data.mat');
save(anim_data_file, 'tspan', 'X_history', 'visibility_history', ...
    'pointing_eci_history', 'r_volcano_0', 'params');
fprintf('Animation data saved to: %s\n', anim_data_file);
fprintf('Run render_gif.py in the same folder to generate the GIF.\n');

% === RA and Dec vs Time Plots ===
fprintf('Generating RA and Dec plots...\n');
figure('Name', 'RA and Dec vs Time', 'Position', [100 100 1400 800]);

% RA plot
subplot(2, 1, 1);
plot(tspan / 60, ra_camera_history, 'b-', 'LineWidth', 2); hold on;
plot(tspan / 60, ra_volcano_history, 'r--', 'LineWidth', 2);
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Right Ascension (degrees)', 'FontSize', 16);
title('Right Ascension vs Time', 'FontSize', 18);
legend('Camera Pointing', 'Apparent Volcano Position', 'Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', 13);
end
xline(t_eruption / 60, 'k-', 'LineWidth', 2, 'Label', 'Eruption (t=0)', 'FontSize', 13);
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', 13);
end

% Dec plot
subplot(2, 1, 2);
plot(tspan / 60, dec_camera_history, 'b-', 'LineWidth', 2); hold on;
plot(tspan / 60, dec_volcano_history, 'r--', 'LineWidth', 2);
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Declination (degrees)', 'FontSize', 16);
title('Declination vs Time', 'FontSize', 18);
legend('Camera Pointing', 'Apparent Volcano Position', 'Location', 'best', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', 13);
end
xline(t_eruption / 60, 'k-', 'LineWidth', 2, 'Label', 'Eruption (t=0)', 'FontSize', 13);
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', 13);
end

% === Zenith Distance Plot ===
fprintf('Generating Zenith Distance plot...\n');
figure('Name', 'Zenith Distance vs Time', 'Position', [100 100 1400 500]);

plot(tspan / 60, zenith_distance_history, 'k-', 'LineWidth', 2);
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Zenith Distance (degrees)', 'FontSize', 16);
title('Zenith Distance of Satellite from Volcano vs Time', 'FontSize', 18);
set(gca, 'FontSize', 14);
grid on;
hold on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', 13);
end
xline(t_eruption / 60, 'k--', 'LineWidth', 2, 'Label', 'Eruption (t=0)', 'FontSize', 13);
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', 13);
end

% Add horizontal line at 90 degrees (horizon)
yline(90, 'r--', 'LineWidth', 1.5, 'Label', 'Horizon (90°)', 'FontSize', 13);

% Shade the visible region (zenith distance < 90 degrees)
ylim([0 180]);
legend('Zenith Distance', 'Location', 'best', 'FontSize', 14);

% === Required Torque Plot ===
fprintf('Generating Required Torque plot...\n');
figure('Name', 'Required Torque (Body Frame) vs Time', 'Position', [100 100 1400 500]);

plot(tspan / 60, tau_body_history(:, 1), 'r-', 'LineWidth', 2.5); hold on;
plot(tspan / 60, tau_body_history(:, 2), 'g-', 'LineWidth', 2.5);
plot(tspan / 60, tau_body_history(:, 3), 'b-', 'LineWidth', 2.5);
xlabel('Time (minutes)', 'FontSize', 24);
ylabel('Torque (N\cdotm)', 'FontSize', 24);
title('Required Tracking Torque (Body Frame) vs Time', 'FontSize', 28);
legend('\tau_x', '\tau_y', '\tau_z', 'Location', 'best', 'FontSize', 20);
set(gca, 'FontSize', 20);
grid on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 2, 'Label', 'Start Visible', 'FontSize', 18);
end
xline(t_eruption / 60, 'k--', 'LineWidth', 2.5, 'Label', 'Eruption (t=0)', 'FontSize', 18);
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 2, 'Label', 'Out of Sight', 'FontSize', 18);
end

% Set the x-axis limits to the visibility window
if ~isnan(t_start_visible) && ~isnan(t_end_visible)
    xlim([t_start_visible / 60, t_end_visible / 60]);
end

% === Verification Simulation ===
if ~isnan(t_start_visible) && ~isnan(t_end_visible)
    fprintf('Running Verification Simulation...\n');
    idx_start = idx_start_visible(1);
    idx_end = idx_end_visible(1);
    
    % Initial state for verification: Perfectly tracking the volcano
    t0_verif = tspan(idx_start);
    r_sat0 = X_history(idx_start, 1:3)';
    v_sat0 = X_history(idx_start, 4:6)';
    
    % Reconstruct the initial attitude that points at the volcano using Minimum Roll frame
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
    R_ECI_to_Body0 = [x_b0, y_b0, z_b0]'; % 3x3 matrix where rows are body axes in ECI
    
    % Use Aerospace Toolbox dcm2quat (scalar first [w x y z])
    q0_aerospace = dcm2quat(R_ECI_to_Body0);
    q0_verif = [q0_aerospace(2:4), q0_aerospace(1)]'; % convert to scalar-last [x y z w]
    
    % Body frame angular velocity (calculated directly from DCM derivative earlier)
    omega0_verif_body = omega_body_history(idx_start, :)';
    % State expects ECI frame angular velocity
    omega0_verif_eci = R_ECI_to_Body0' * omega0_verif_body;
    
    X0_verif = [r_sat0; v_sat0; q0_verif; omega0_verif_eci];
    
    % Interpolator for open-loop torque (Body Frame) using spline for smoother and more accurate torque curves
    tau_interp = @(t) interp1(tspan, tau_body_history, t, 'spline')';
    
    % ODE Function for verification (open-loop torque)
    ode_fun_verif = @(t, X) sat_dynamics_openloop(t, X, params, tau_interp);
    
    tspan_verif = tspan(idx_start:idx_end);
    [~, X_verif] = ode45(ode_fun_verif, tspan_verif, X0_verif, ode_opts);
    
    % Calculate pointing RA/Dec and compare with apparent Volcano RA/Dec
    ra_verif_history = zeros(length(tspan_verif), 1);
    dec_verif_history = zeros(length(tspan_verif), 1);
    
    for k = 1:length(tspan_verif)
        obs_verif = state_to_observation(X_verif(k, :)', params);
        ra_verif_history(k) = obs_verif.ra;
        dec_verif_history(k) = obs_verif.dec;
    end
    
    % Get the apparent volcano positions for the same window
    ra_target_history = ra_volcano_history(idx_start:idx_end);
    dec_target_history = dec_volcano_history(idx_start:idx_end);
    
    % We need to handle 360 to 0 wrap around for RA error
    ra_error = ra_verif_history - ra_target_history;
    ra_error = mod(ra_error + 180, 360) - 180; % Wrap to [-180, 180]
    dec_error = dec_verif_history - dec_target_history;
    
    figure('Name', 'Verification Pointing Error', 'Position', [150 150 1200 600]);
    plot(tspan_verif / 60, ra_error, 'r-', 'LineWidth', 2); hold on;
    plot(tspan_verif / 60, dec_error, 'b--', 'LineWidth', 2);
    xlabel('Time (minutes)', 'FontSize', 20);
    ylabel('Pointing Error (degrees)', 'FontSize', 20);
    title('Verification: Pointing Error with Open-Loop Torque', 'FontSize', 24);
    legend('RA Error', 'Dec Error', 'Location', 'best', 'FontSize', 16);
    set(gca, 'FontSize', 16);
    grid on;
    xlim([tspan_verif(1)/60, tspan_verif(end)/60]);
    
    fprintf('Max RA Error: %.2e deg\n', max(abs(ra_error)));
    fprintf('Max Dec Error: %.2e deg\n', max(abs(dec_error)));
end

fprintf('\nSimulation complete!\n');

% --- Helper Function for Verification ---
function Xd = sat_dynamics_openloop(t, X, params, tau_interp)
    % Extract parameters
    mu = params.mu;
    I_CB = params.I_CB; 
    
    % Extract state components
    r_eci = X(1:3);
    v_eci = X(4:6);
    q_eci_to_body = X(7:10);
    omega_eci = X(11:13);
    
    q_eci_to_body = q_eci_to_body / norm(q_eci_to_body);
    
    % DCMs using Aerospace Toolbox quat2dcm
    qx = q_eci_to_body(1); qy = q_eci_to_body(2); qz = q_eci_to_body(3); qw = q_eci_to_body(4);
    q_aerospace = [qw, qx, qy, qz]; % scalar-first
    R_ECI_to_Body = quat2dcm(q_aerospace);
    R_Body_to_ECI = R_ECI_to_Body';
    
    omega_body = R_ECI_to_Body * omega_eci;
    
    Xd = zeros(13, 1);
    Xd(1:3) = v_eci;
    r_norm = norm(r_eci);
    Xd(4:6) = -mu * r_eci / r_norm^3;
    
    Omega_matrix = [ 0, omega_body(3), -omega_body(2), omega_body(1);
                    -omega_body(3), 0, omega_body(1), omega_body(2);
                     omega_body(2), -omega_body(1), 0, omega_body(3);
                    -omega_body(1), -omega_body(2), -omega_body(3), 0 ];
    Xd(7:10) = 0.5 * Omega_matrix * q_eci_to_body;
    
    tau_body = tau_interp(t);
    
    alpha_body = I_CB \ (tau_body - cross(omega_body, I_CB * omega_body));
    Xd(11:13) = R_Body_to_ECI * alpha_body;
end
