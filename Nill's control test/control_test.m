% Modified Supernova Slew Maneuver with Smooth Proportional Control and 3D Animation
% Mission: Point to supernova in 5 minutes, observe behavior for 6 minutes total

clear; close all; clc;

% Parameters
params.mu = 3.986004418e14; % Earth gravitational parameter [m^3/s^2]
params.R_e = 6378137; % Earth radius [m]
params.mass = 83.6; % Satellite mass [kg]
params.radius = 0.58/2; % Satellite radius [m] (for thrust calculation)
params.J2 = 1.08263e-3; % J2 coefficient
params.I_CB = 2/5 * params.mass * params.radius^2 * eye(3); % Moment of inertia [kg*m^2]

% Proportional control gain
params.Kp_att = 0.5; % [N*m]

% Initial Conditions
altitude = 200e3; % 200 km
r_orbit = params.R_e + altitude;
v_circ = sqrt(params.mu / r_orbit);

% Orbital plane random orientation
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

% Initial attitude - nadir pointing basis
orbit_normal_eci = cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
z_body_eci = -r_eci / norm(r_eci);
x_body_eci = v_eci / norm(v_eci);
y_body_eci = cross(z_body_eci, x_body_eci);
R_Body_to_ECI = [x_body_eci, y_body_eci, z_body_eci];
q0 = dcm_to_quaternion(R_Body_to_ECI);
omega_orbit_mag = v_circ / r_orbit;
omega_eci = omega_orbit_mag * orbit_normal_eci;
X0 = [r_eci; v_eci; q0; omega_eci];

% Supernova event
supernova_ra = rand() * 360; % RA degrees
supernova_dec = (rand() - 0.5) * 180; % Dec degrees
ra_rad = deg2rad(supernova_ra);
dec_rad = deg2rad(supernova_dec);
actual_target = [cos(dec_rad)*cos(ra_rad); cos(dec_rad)*sin(ra_rad); sin(dec_rad)];
actual_target = actual_target / norm(actual_target);

% Original pointing vector (at t=0) - camera aligned +Z body axis
original_point_vec = R_Body_to_ECI(:,3); % +Z in ECI

fprintf('=== SUPERNOVA SLEW MANEUVER ===\n');
fprintf('Supernova Location: RA=%.2f deg, Dec=%.2f deg\n', supernova_ra, supernova_dec);
fprintf('Mission: Reach target in 5 minutes\n');
fprintf('Simulation: Run for 6 minutes to observe behavior\n\n');

total_angle = acos(max(min(dot(original_point_vec, actual_target), 1), -1));
rotation_axis = cross(original_point_vec, actual_target);
if norm(rotation_axis) < 1e-6
    rotation_axis = [0; 0; 1];
else
    rotation_axis = rotation_axis / norm(rotation_axis);
end

% Timing parameters
maneuver_time = 5 * 60; % 5 minutes to reach target
total_sim_time = 6 * 60; % 6 minutes total observation
sim_time = 20; % 20 seconds wall-clock time
speedup_factor = total_sim_time / sim_time;

dt = 0.1;
tspan = 0:dt:sim_time;
n_steps = length(tspan);

X_history = zeros(n_steps, 13);
pointing_error_history = zeros(n_steps, 1);
tau_history = zeros(n_steps, 3);
thrust_history = zeros(n_steps, 3);
omega_mag_history = zeros(n_steps, 1);
alpha_mag_history = zeros(n_steps, 1);
desired_point_history = zeros(n_steps, 3);

X_history(1,:) = X0';
event_acquired = false;
t_acquired = NaN;

fprintf('Starting controlled slew with smooth target trajectory...\n');

for i = 2:n_steps
    t_real = tspan(i) * speedup_factor; % Real mission time

    % Desired pointing vector based on 5-minute maneuver time
    if t_real < maneuver_time
        alpha = t_real / maneuver_time;
    else
        alpha = 1; % After 5 minutes, stay at target
    end
    desired_point_vec = (1-alpha)*original_point_vec + alpha*actual_target;
    desired_point_vec = desired_point_vec / norm(desired_point_vec);
    desired_point_history(i, :) = desired_point_vec';

    % Propagate one step with desired pointing
    [~, X_seg] = ode45(@(t,X) Sat_template(t, X, desired_point_vec, params, ...
        'useJ2', false, 'useAtmDrag', false, 'useControl', true), ...
        [tspan(i-1), tspan(i)], X_history(i-1,:)', ...
        odeset('RelTol',1e-8, 'AbsTol',1e-10));
    X_history(i,:) = X_seg(end,:);

    % Compute pointing error to actual target
    obs = state_to_observation(X_history(i,:)', params);
    current_pointing_eci = obs.pointing_eci;
    dot_prod = dot(current_pointing_eci, actual_target);
    dot_prod = min(max(dot_prod, -1), 1); % Clamp
    pointing_error_deg = rad2deg(acos(dot_prod));
    pointing_error_history(i) = pointing_error_deg;

    % Torque logging
    q = X_history(i, 7:10)';
    q = q / norm(q);
    R_Body_to_ECI = quaternion_to_dcm(q);
    R_ECI_to_Body = R_Body_to_ECI';
    pointing_error_eci = cross(current_pointing_eci, desired_point_vec);
    pointing_error_body = R_ECI_to_Body * pointing_error_eci;
    tau_body = params.Kp_att * pointing_error_body;
    tau_eci = R_Body_to_ECI * tau_body;
    tau_history(i, :) = tau_eci';

    % Thrust logging
    thrust_history(i, :) = tau_eci' / (2 * params.radius);

    omega_eci = X_history(i, 11:13)';
    omega_mag_history(i) = norm(omega_eci);

    % Angular acceleration approx
    if i > 2
        alpha_eci = (X_history(i, 11:13) - X_history(i-1, 11:13)) / dt;
        alpha_mag_history(i) = norm(alpha_eci);
    end

    % Acquisition check (only report once)
    if pointing_error_deg < 0.5 && ~event_acquired
        event_acquired = true;
        t_acquired = tspan(i);
        fprintf('*** TARGET ACQUIRED at t=%.2f s (real time %.2f min) ***\n', ...
            t_acquired, t_real / 60);
        fprintf('Pointing error: %.4f degrees\n\n', pointing_error_deg);
    end
end

% Store initial desired pointing
desired_point_history(1, :) = original_point_vec';

if ~event_acquired
    fprintf('*** MISSION FAILED to acquire target within 6 minutes ***\n');
    fprintf('Final pointing error: %.4f degrees\n\n', pointing_error_history(end));
end

% === 3D Animated Plot ===
fprintf('Generating 3D animated plot...\n');

figure('Name', '3D Attitude Visualization', 'Position', [100 100 1200 800]);

[xs, ys, zs] = sphere(50);

n_frames = length(tspan);
frame_skip = max(1, floor(n_frames / 100));

for i = 1:frame_skip:n_frames
    clf;
    surf(xs, ys, zs, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);
    hold on;
    axis equal;
    grid on;
    xlabel('X (ECI)'); ylabel('Y (ECI)'); zlabel('Z (ECI)');
    t_real = tspan(i) * speedup_factor;
    title(sprintf('Satellite Attitude Control - Time: %.2f s (Real: %.2f min)', ...
        tspan(i), t_real/60));

    % Current pointing vector (blue)
    obs = state_to_observation(X_history(i,:)', params);
    current_pointing = obs.pointing_eci;
    quiver3(0, 0, 0, current_pointing(1), current_pointing(2), current_pointing(3), ...
        1.2, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5);

    % Actual target pointing vector (red)
    quiver3(0, 0, 0, actual_target(1), actual_target(2), actual_target(3), ...
        1.2, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5);

    % Desired pointing vector (yellow/gold) - time-varying
    desired_point_vec = desired_point_history(i, :)';
    quiver3(0, 0, 0, desired_point_vec(1), desired_point_vec(2), desired_point_vec(3), ...
        1.2, 'Color', [1 0.6 0], 'LineWidth', 3, 'MaxHeadSize', 0.5);

    % Angular velocity vector (green) - normalized
    omega_eci = X_history(i, 11:13)';
    if norm(omega_eci) > 1e-6
        omega_unit = omega_eci / norm(omega_eci);
        quiver3(0, 0, 0, omega_unit(1), omega_unit(2), omega_unit(3), ...
            1.0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    end

    % Angular acceleration vector (magenta) - normalized
    if i > 1
        alpha_eci = (X_history(i, 11:13) - X_history(i-1, 11:13)) / dt;
        if norm(alpha_eci) > 1e-6
            alpha_unit = alpha_eci / norm(alpha_eci);
            quiver3(0, 0, 0, alpha_unit(1), alpha_unit(2), alpha_unit(3), ...
                0.9, 'm', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        end
    end

    % Control torque vector (cyan) - normalized
    tau_eci = tau_history(i, :)';
    if norm(tau_eci) > 1e-6
        tau_unit = tau_eci / norm(tau_eci);
        quiver3(0, 0, 0, tau_unit(1), tau_unit(2), tau_unit(3), ...
            0.8, 'c', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    end

    legend('Unit Sphere', 'Current Pointing', 'Target Pointing (Final)', ...
        'Desired Pointing (Traj)', 'Angular Velocity', 'Angular Acceleration', ...
        'Control Torque', 'Location', 'best');

    view(45, 30);
    xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);

    drawnow;
    pause(0.05);
end

% === Pointing Error Plot ===
figure('Name', 'Pointing Error');
plot(tspan * speedup_factor / 60, pointing_error_history, 'b-', 'LineWidth', 2);
hold on;
yline(0.5, 'r--', 'LineWidth', 1.5, 'Label', 'FOV Tolerance (0.5°)');
xline(5, 'k--', 'LineWidth', 1.5, 'Label', '5 min (Target Time)');
xlabel('Real Time (minutes)');
ylabel('Pointing Error (degrees)');
title('Pointing Error to Final Target Over Time');
grid on;
xlim([0, 6]);

% === Torque / Thrust Plot ===
figure('Name', 'Control Effort', 'Position', [100 100 1200 800]);

subplot(2, 1, 1);
time_min = tspan * speedup_factor / 60;
plot(time_min, tau_history(:, 1), 'r-', 'LineWidth', 1.5); hold on;
plot(time_min, tau_history(:, 2), 'g-', 'LineWidth', 1.5);
plot(time_min, tau_history(:, 3), 'b-', 'LineWidth', 1.5);
plot(time_min, vecnorm(tau_history, 2, 2), 'k--', 'LineWidth', 2);
xline(5, 'k--', 'LineWidth', 1, 'Label', '5 min');
xlabel('Real Time (minutes)'); 
ylabel('Torque (N·m)');
title('Control Torque vs Time');
legend('\tau_x', '\tau_y', '\tau_z', '||\tau||', 'Location', 'best');
grid on;
xlim([0, 6]);

subplot(2, 1, 2);
plot(time_min, thrust_history(:, 1), 'r-', 'LineWidth', 1.5); hold on;
plot(time_min, thrust_history(:, 2), 'g-', 'LineWidth', 1.5);
plot(time_min, thrust_history(:, 3), 'b-', 'LineWidth', 1.5);
plot(time_min, vecnorm(thrust_history, 2, 2), 'k--', 'LineWidth', 2);
xline(5, 'k--', 'LineWidth', 1, 'Label', '5 min');
xlabel('Real Time (minutes)'); 
ylabel('Thrust (N)');
title('Control Thrust vs Time (Thrust = Torque / (2×radius))');
legend('F_x', 'F_y', 'F_z', '||F||', 'Location', 'best');
grid on;
xlim([0, 6]);

fprintf('\nSimulation complete!\n');

% === Helper Functions ===
function R = quaternion_to_dcm(q)
    qx = q(1); qy = q(2); qz = q(3); qw = q(4);
    R = [1-2*(qy^2+qz^2), 2*(qx*qy-qw*qz), 2*(qx*qz+qw*qy); 
         2*(qx*qy+qw*qz), 1-2*(qx^2+qz^2), 2*(qy*qz-qw*qx); 
         2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx), 1-2*(qx^2+qy^2)];
end
