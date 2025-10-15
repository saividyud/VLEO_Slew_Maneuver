%% Quaternion Test Using Sat_template - Azimuth/Elevation Tracking
% This test verifies that body-fixed vectors rotate correctly according to
% R(omega*t) when integrated through Sat_template's quaternion kinematics.

clear; clc; close all;

%% Physical Parameters

params.mass = 100; % kg
params.mu = 3.986004418e14; % m^3/s^2
params.R_e = 6378137; % m
params.I_CB = diag([100, 100, 100]); % Symmetric inertia for simplicity

%% Define Constant Angular Velocity

omega_magnitude = 0.1; % rad/s (about 5.7 deg/s)
omega_axis = [1; 1; 1]; % Arbitrary axis
omega_axis = omega_axis / norm(omega_axis); % Normalize
omega = omega_magnitude * omega_axis;

fprintf('=== Quaternion Kinematics Test with Sat_template ===\n\n');
fprintf('Angular Velocity (Body Frame):\n');
fprintf('  ω = [%.6f, %.6f, %.6f] rad/s\n', omega(1), omega(2), omega(3));
fprintf('  |ω| = %.6f rad/s (%.2f deg/s)\n', omega_magnitude, rad2deg(omega_magnitude));
fprintf('  Rotation axis: [%.4f, %.4f, %.4f]\n\n', omega_axis(1), omega_axis(2), omega_axis(3));

%% Initial Conditions - Stationary in Space

% Position: doesn't matter for this test, use high orbit to minimize gravity gradient
altitude = 1000e3; % 1000 km
r0 = params.R_e + altitude;
r_vec0 = [r0; 0; 0];

% Velocity: minimal (essentially stationary in inertial frame)
v_vec0 = [0; 0; 0];

% Initial quaternion: Identity (body frame = inertial frame initially)
q0 = [0; 0; 0; 1]; % [qx, qy, qz, qw]

% Initial angular velocity: our test value
omega0 = omega;

% Assemble state
X0 = [r_vec0; v_vec0; q0; omega0];

%% Define Body-Fixed Reference Vectors

% These vectors are fixed in the body frame
% We'll track how they appear in the inertial frame as the body rotates

v_body1 = [1; 0; 0]; % Body X-axis
v_body2 = [0; 1; 0]; % Body Y-axis  
v_body3 = [0; 0; 1]; % Body Z-axis

fprintf('Body-Fixed Reference Vectors:\n');
fprintf('  v1 = [%.0f, %.0f, %.0f] (Body X)\n', v_body1(1), v_body1(2), v_body1(3));
fprintf('  v2 = [%.0f, %.0f, %.0f] (Body Y)\n', v_body2(1), v_body2(2), v_body2(3));
fprintf('  v3 = [%.0f, %.0f, %.0f] (Body Z)\n\n', v_body3(1), v_body3(2), v_body3(3));

%% Simulation Setup

% Simulate for one complete rotation
period = 2 * pi / omega_magnitude;
t_span = linspace(0, period, 200);

fprintf('Simulating for %.2f seconds (1 rotation period)...\n', period);

%% Run Simulation (Torque-Free, No External Forces)

desired_state = X0; % No control needed

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

[t, X] = ode45(@(t, X) Sat_template(t, X, desired_state, params, ...
                                     'useJ2', false, ...
                                     'useAtmDrag', false, ...
                                     'useControl', false, ...
                                     'useSRP', false), ...
               t_span, X0, options);

fprintf('Simulation complete.\n\n');

%% Extract Results

q_numeric = X(:, 7:10); % [qx, qy, qz, qw]
omega_numeric = X(:, 11:13);

%% Compute Analytical Solution

fprintf('Computing analytical solution...\n');

% For constant omega, the quaternion evolves as rotation about omega axis
q_analytic = zeros(length(t), 4);
v1_numeric_inertial = zeros(length(t), 3);
v1_analytic_inertial = zeros(length(t), 3);
v2_numeric_inertial = zeros(length(t), 3);
v3_numeric_inertial = zeros(length(t), 3);

for i = 1:length(t)
    % Analytical quaternion at time t
    theta = omega_magnitude * t(i);
    q_analytic(i, :) = axis_angle_to_quat(omega_axis, theta)';

    % Rotate body-fixed vectors to inertial frame using numeric quaternion
    R_numeric = quat_to_rotation_matrix(q_numeric(i, :)');
    v1_numeric_inertial(i, :) = (R_numeric * v_body1)';
    v2_numeric_inertial(i, :) = (R_numeric * v_body2)';
    v3_numeric_inertial(i, :) = (R_numeric * v_body3)';

    % Analytical rotation: R(ω*t)
    R_analytic = rodrigues_rotation(omega_axis, theta);
    v1_analytic_inertial(i, :) = (R_analytic * v_body1)';
end

fprintf('Analytical computation complete.\n\n');

%% Compute Azimuth and Elevation for v1

% Azimuth: angle in XY plane from X-axis
% Elevation: angle from XY plane toward Z-axis

az_numeric = atan2(v1_numeric_inertial(:, 2), v1_numeric_inertial(:, 1));
el_numeric = asin(v1_numeric_inertial(:, 3));

az_analytic = atan2(v1_analytic_inertial(:, 2), v1_analytic_inertial(:, 1));
el_analytic = asin(v1_analytic_inertial(:, 3));

%% Error Analysis

% Quaternion error
q_error = q_numeric - q_analytic;
q_error_norm = vecnorm(q_error, 2, 2);

% Vector rotation error for v1
v1_error = v1_numeric_inertial - v1_analytic_inertial;
v1_error_norm = vecnorm(v1_error, 2, 2);

% Azimuth/Elevation error
az_error = rad2deg(wrapToPi(az_numeric - az_analytic));
el_error = rad2deg(el_numeric - el_analytic);

% Angular velocity constancy check
omega_deviation = omega_numeric - repmat(omega0', length(t), 1);
omega_dev_norm = vecnorm(omega_deviation, 2, 2);

fprintf('=== Error Analysis ===\n');
fprintf('Quaternion Error:\n');
fprintf('  Max norm error: %.3e\n', max(q_error_norm));
fprintf('  Mean norm error: %.3e\n\n', mean(q_error_norm));

fprintf('Vector Rotation Error (v1):\n');
fprintf('  Max norm error: %.3e\n', max(v1_error_norm));
fprintf('  Mean norm error: %.3e\n\n', mean(v1_error_norm));

fprintf('Azimuth/Elevation Error (v1):\n');
fprintf('  Max azimuth error: %.4f deg\n', max(abs(az_error)));
fprintf('  Max elevation error: %.4f deg\n\n', max(abs(el_error)));

fprintf('Angular Velocity Constancy:\n');
fprintf('  Max deviation from ω₀: %.3e rad/s\n', max(omega_dev_norm));
fprintf('  Mean deviation: %.3e rad/s\n\n', mean(omega_dev_norm));

%% Pass/Fail Criteria

fprintf('=== Verification Results ===\n');

tolerance_quat = 1e-4;
tolerance_vector = 1e-4;
tolerance_angle = 0.01; % degrees
tolerance_omega = 1e-8;

pass_count = 0;
total_tests = 4;

if max(q_error_norm) < tolerance_quat
    fprintf('✓ PASS: Quaternion integration (error: %.3e)\n', max(q_error_norm));
    pass_count = pass_count + 1;
else
    fprintf('✗ FAIL: Quaternion integration (error: %.3e > %.3e)\n', max(q_error_norm), tolerance_quat);
end

if max(v1_error_norm) < tolerance_vector
    fprintf('✓ PASS: Vector rotation (error: %.3e)\n', max(v1_error_norm));
    pass_count = pass_count + 1;
else
    fprintf('✗ FAIL: Vector rotation (error: %.3e > %.3e)\n', max(v1_error_norm), tolerance_vector);
end

if max(abs(az_error)) < tolerance_angle && max(abs(el_error)) < tolerance_angle
    fprintf('✓ PASS: Az/El tracking (az: %.4f°, el: %.4f°)\n', max(abs(az_error)), max(abs(el_error)));
    pass_count = pass_count + 1;
else
    fprintf('✗ FAIL: Az/El tracking exceeds %.2f degrees\n', tolerance_angle);
end

if max(omega_dev_norm) < tolerance_omega
    fprintf('✓ PASS: Angular velocity constancy (dev: %.3e)\n', max(omega_dev_norm));
    pass_count = pass_count + 1;
else
    fprintf('✗ FAIL: Angular velocity not constant (dev: %.3e)\n', max(omega_dev_norm));
end

fprintf('\n=== Overall: %d/%d tests passed ===\n\n', pass_count, total_tests);

%% Visualization

% Plot 1: 3D Trajectory of Body X-axis in Inertial Frame
figure('Position', [100, 100, 900, 700]);
plot3(v1_numeric_inertial(:, 1), v1_numeric_inertial(:, 2), v1_numeric_inertial(:, 3), ...
      'b-', 'LineWidth', 2); hold on;
plot3(v1_analytic_inertial(:, 1), v1_analytic_inertial(:, 2), v1_analytic_inertial(:, 3), ...
      'r--', 'LineWidth', 2);

% Mark start/end
plot3(v1_numeric_inertial(1, 1), v1_numeric_inertial(1, 2), v1_numeric_inertial(1, 3), ...
      'go', 'MarkerSize', 12, 'LineWidth', 2);
plot3(v1_numeric_inertial(end, 1), v1_numeric_inertial(end, 2), v1_numeric_inertial(end, 3), ...
      'ro', 'MarkerSize', 12, 'LineWidth', 2);

% Draw rotation axis
quiver3(0, 0, 0, omega_axis(1), omega_axis(2), omega_axis(3), ...
        'k', 'LineWidth', 3, 'MaxHeadSize', 0.5);

% Unit sphere
[X_sp, Y_sp, Z_sp] = sphere(40);
surf(X_sp, Y_sp, Z_sp, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

xlabel('X (Inertial)', 'FontSize', 12);
ylabel('Y (Inertial)', 'FontSize', 12);
zlabel('Z (Inertial)', 'FontSize', 12);
title('Body X-axis Trajectory in Inertial Frame', 'FontSize', 14, 'FontWeight', 'bold');
legend('Sat\_template', 'Analytic R(\omega t)', 'Start', 'End', 'Rotation Axis', ...
       'Location', 'best');
grid on;
axis equal;
view(45, 30);

% Plot 2: Azimuth and Elevation vs Time
figure('Position', [150, 150, 1200, 500]);

subplot(1, 2, 1);
plot(t, rad2deg(az_numeric), 'b-', 'LineWidth', 1.5); hold on;
plot(t, rad2deg(az_analytic), 'r--', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Azimuth (deg)', 'FontSize', 12);
title('Azimuth Angle of Body X-axis', 'FontSize', 13);
legend('Sat\_template', 'Analytic', 'Location', 'best');
grid on;

subplot(1, 2, 2);
plot(t, rad2deg(el_numeric), 'b-', 'LineWidth', 1.5); hold on;
plot(t, rad2deg(el_analytic), 'r--', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Elevation (deg)', 'FontSize', 12);
title('Elevation Angle of Body X-axis', 'FontSize', 13);
legend('Sat\_template', 'Analytic', 'Location', 'best');
grid on;

sgtitle('Azimuth/Elevation Tracking', 'FontSize', 15, 'FontWeight', 'bold');

% Plot 3: Az/El Errors
figure('Position', [200, 200, 1200, 500]);

subplot(1, 2, 1);
plot(t, az_error, 'r-', 'LineWidth', 1.5); hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Error (deg)', 'FontSize', 12);
title('Azimuth Error', 'FontSize', 13);
grid on;

subplot(1, 2, 2);
plot(t, el_error, 'b-', 'LineWidth', 1.5); hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Error (deg)', 'FontSize', 12);
title('Elevation Error', 'FontSize', 13);
grid on;

sgtitle('Az/El Tracking Errors', 'FontSize', 15, 'FontWeight', 'bold');

% Plot 4: Body Frame Axes in Inertial Frame at Selected Times
figure('Position', [250, 250, 900, 700]);

% Select 8 time points evenly distributed
n_snapshots = 8;
snapshot_indices = round(linspace(1, length(t), n_snapshots));

colors = lines(n_snapshots);

for i = 1:n_snapshots
    idx = snapshot_indices(i);
    R = quat_to_rotation_matrix(q_numeric(idx, :)');

    % Plot body axes
    scale = 0.3;
    quiver3(0, 0, 0, R(1,1), R(2,1), R(3,1), scale, ...
            'Color', colors(i,:), 'LineWidth', 2, 'MaxHeadSize', 0.3); hold on;
    quiver3(0, 0, 0, R(1,2), R(2,2), R(3,2), scale, ...
            'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', '--', 'MaxHeadSize', 0.3);
    quiver3(0, 0, 0, R(1,3), R(2,3), R(3,3), scale, ...
            'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', ':', 'MaxHeadSize', 0.3);
end

% Draw rotation axis
quiver3(0, 0, 0, omega_axis(1), omega_axis(2), omega_axis(3), ...
        'k', 'LineWidth', 3, 'MaxHeadSize', 0.5);

xlabel('X', 'FontSize', 12);
ylabel('Y', 'FontSize', 12);
zlabel('Z', 'FontSize', 12);
title('Body Frame Orientation at Different Times', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis equal;
axis([-1 1 -1 1 -1 1]);
view(45, 30);

%% Helper Functions

function q = axis_angle_to_quat(axis, angle)
    % Convert axis-angle to quaternion (scalar-last: qx, qy, qz, qw)
    axis = axis / norm(axis);
    half_angle = angle / 2;
    q = [axis * sin(half_angle); cos(half_angle)];
end

function R = rodrigues_rotation(axis, angle)
    % Rodrigues' rotation formula: R = I + sin(θ)K + (1-cos(θ))K²
    axis = axis / norm(axis);
    K = [0,        -axis(3),  axis(2);
         axis(3),   0,       -axis(1);
        -axis(2),   axis(1),  0];
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
end

function R = quat_to_rotation_matrix(q)
    % Convert quaternion to rotation matrix (scalar-last convention)
    % q = [qx, qy, qz, qw]
    qx = q(1); qy = q(2); qz = q(3); qw = q(4);

    R = [1 - 2*(qy^2 + qz^2),     2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw);
         2*(qx*qy + qz*qw),       1 - 2*(qx^2 + qz^2),   2*(qy*qz - qx*qw);
         2*(qx*qz - qy*qw),       2*(qy*qz + qx*qw),     1 - 2*(qx^2 + qy^2)];
end
