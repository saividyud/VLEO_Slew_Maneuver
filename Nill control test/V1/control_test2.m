% Control Test Earth - Volcano Observation Scenario
% Simulates satellite observation of a volcano eruption with Earth rotation
% Volcano erupts at t=0 and is guaranteed to be visible at that time

clear; close all; clc;

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
q0 = dcm_to_quaternion(R_Body_to_ECI);

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
dt = 10; % 10 second time steps
tspan = t_start:dt:t_end;
n_steps = length(tspan);

fprintf('Simulation from t = %.2f to %.2f seconds\n', t_start, t_end);
fprintf('Total simulation time: %.2f minutes\n\n', (t_end - t_start) / 60);

% Propagate the satellite state across the entire time window
fprintf('Propagating satellite trajectory...\n');
X_history = zeros(n_steps, 13);
ra_camera_history = zeros(n_steps, 1);
dec_camera_history = zeros(n_steps, 1);
ra_volcano_history = zeros(n_steps, 1);
dec_volcano_history = zeros(n_steps, 1);
visibility_history = zeros(n_steps, 1);
zenith_distance_history = zeros(n_steps, 1);

for i = 1:n_steps
    t = tspan(i);

    % Propagate satellite state from t=0
    if t == 0
        X_history(i, :) = X0';
    else
        [~, X_seg] = ode45(@(t,X) Sat_template(t, X, z_body_eci, params, ...
            'useJ2', false, 'useAtmDrag', false, 'useControl', false), ...
            [0, t], X0, odeset('RelTol',1e-8, 'AbsTol',1e-10));
        X_history(i, :) = X_seg(end, :);
    end

    % Get camera pointing direction (RA/Dec)
    obs = state_to_observation(X_history(i, :)', params);
    ra_camera_history(i) = obs.ra;
    dec_camera_history(i) = obs.dec;

    % Volcano position at time t (rotates with Earth)
    Rotz_t = [cos(params.omega_earth * t), -sin(params.omega_earth * t), 0;
              sin(params.omega_earth * t),  cos(params.omega_earth * t), 0;
              0, 0, 1];
    r_volcano_t = Rotz_t * r_volcano_0;

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

    % Calculate zenith distance
    % Zenith direction from volcano (pointing up from volcano)
    zenith_at_volcano = r_volcano_t / norm(r_volcano_t);

    % Direction from volcano to satellite
    volcano_to_sat = r_sat - r_volcano_t;
    volcano_to_sat_unit = volcano_to_sat / norm(volcano_to_sat);

    % Zenith distance is the angle between zenith and satellite direction
    cos_zenith = dot(zenith_at_volcano, volcano_to_sat_unit);
    zenith_distance_deg = rad2deg(acos(cos_zenith));
    zenith_distance_history(i) = zenith_distance_deg;
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

% === 3D Animated Plot ===
fprintf('Generating 3D animated trajectory...\n');
figure('Name', '3D Satellite and Volcano Trajectory', 'Position', [100 100 1200 800]);

% Earth sphere
[xs, ys, zs] = sphere(50);
xs = xs * params.R_e;
ys = ys * params.R_e;
zs = zs * params.R_e;

% Animation
n_frames = n_steps;
frame_skip = max(1, floor(n_frames / 25));

for i = 1:frame_skip:n_frames
    clf;

    % Draw Earth
    surf(xs, ys, zs, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.3 0.5 0.8]);
    hold on;
    axis equal;
    grid on;
    xlabel('X (ECI) [m]'); ylabel('Y (ECI) [m]'); zlabel('Z (ECI) [m]');

    t = tspan(i);
    title(sprintf('Volcano Observation - Time: %.2f min (Eruption at t=0)', t/60));

    % Satellite position
    r_sat = X_history(i, 1:3);
    plot3(r_sat(1), r_sat(2), r_sat(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    % Satellite trajectory (past)
    if i > 1
        plot3(X_history(1:i, 1), X_history(1:i, 2), X_history(1:i, 3), ...
            'r-', 'LineWidth', 1.5);
    end

    % Volcano position at time t
    Rotz_t = [cos(params.omega_earth * t), -sin(params.omega_earth * t), 0;
              sin(params.omega_earth * t),  cos(params.omega_earth * t), 0;
              0, 0, 1];
    r_volcano_t = Rotz_t * r_volcano_0;
    plot3(r_volcano_t(1), r_volcano_t(2), r_volcano_t(3), ...
        'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'y');

    % Line from volcano to satellite with dynamic color based on visibility
    if visibility_history(i) > 0
        % Visible - Green line
        plot3([r_volcano_t(1), r_sat(1)], [r_volcano_t(2), r_sat(2)], ...
            [r_volcano_t(3), r_sat(3)], 'g-', 'LineWidth', 2.5);
        line_label = 'Line (Visible)';
    else
        % Not visible - Red line
        plot3([r_volcano_t(1), r_sat(1)], [r_volcano_t(2), r_sat(2)], ...
            [r_volcano_t(3), r_sat(3)], 'r--', 'LineWidth', 1.5);
        line_label = 'Line (Not Visible)';
    end

    % Camera pointing direction
    obs = state_to_observation(X_history(i, :)', params);
    pointing_vec = obs.pointing_eci * 2e6; % Scale for visibility
    quiver3(r_sat(1), r_sat(2), r_sat(3), ...
        pointing_vec(1), pointing_vec(2), pointing_vec(3), ...
        'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);

    legend('Earth', 'Satellite', 'Trajectory', 'Volcano', line_label, ...
        'Camera Pointing', 'Location', 'best');
    view(45, 30);

    drawnow;
end

% === RA and Dec vs Time Plots ===
fprintf('Generating RA and Dec plots...\n');
figure('Name', 'RA and Dec vs Time', 'Position', [100 100 1400 800]);

% RA plot
subplot(2, 1, 1);
plot(tspan / 60, ra_camera_history, 'b-', 'LineWidth', 2); hold on;
plot(tspan / 60, ra_volcano_history, 'r--', 'LineWidth', 2);
xlabel('Time (minutes)');
ylabel('Right Ascension (degrees)');
title('Right Ascension vs Time');
legend('Camera Pointing', 'Apparent Volcano Position', 'Location', 'best');
grid on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible');
end
xline(t_eruption / 60, 'k-', 'LineWidth', 2, 'Label', 'Eruption (t=0)');
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight');
end

% Dec plot
subplot(2, 1, 2);
plot(tspan / 60, dec_camera_history, 'b-', 'LineWidth', 2); hold on;
plot(tspan / 60, dec_volcano_history, 'r--', 'LineWidth', 2);
xlabel('Time (minutes)');
ylabel('Declination (degrees)');
title('Declination vs Time');
legend('Camera Pointing', 'Apparent Volcano Position', 'Location', 'best');
grid on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible');
end
xline(t_eruption / 60, 'k-', 'LineWidth', 2, 'Label', 'Eruption (t=0)');
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight');
end

% === Zenith Distance Plot ===
fprintf('Generating Zenith Distance plot...\n');
figure('Name', 'Zenith Distance vs Time', 'Position', [100 100 1400 500]);

plot(tspan / 60, zenith_distance_history, 'k-', 'LineWidth', 2);
xlabel('Time (minutes)');
ylabel('Zenith Distance (degrees)');
title('Zenith Distance of Satellite from Volcano vs Time');
grid on;
hold on;

% Mark events
if ~isnan(t_start_visible)
    xline(t_start_visible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible');
end
xline(t_eruption / 60, 'k--', 'LineWidth', 2, 'Label', 'Eruption (t=0)');
if ~isnan(t_end_visible)
    xline(t_end_visible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight');
end

% Add horizontal line at 90 degrees (horizon)
yline(90, 'r--', 'LineWidth', 1.5, 'Label', 'Horizon (90°)');

% Shade the visible region (zenith distance < 90 degrees)
ylim([0 180]);
legend('Zenith Distance', 'Location', 'best');

fprintf('\nSimulation complete!\n');
