%% Quaternion Kinematics Verification Test
% Test: For constant angular velocity omega, verify that the rotation 
% applied to body-fixed vectors matches the analytical solution R(omega*t)

clear; clc; close all;

%% Test Parameters

% Define constant angular velocity (rad/s)
omega_magnitude = 0.1; % rad/s
omega_axis = [1; 2; 3]; % Arbitrary rotation axis
omega_axis = omega_axis / norm(omega_axis); % Normalize
omega = omega_magnitude * omega_axis; % Angular velocity vector

fprintf('=== Quaternion Kinematics Verification Test ===\n\n');
fprintf('Angular Velocity:\n');
fprintf('  omega = [%.6f, %.6f, %.6f] rad/s\n', omega(1), omega(2), omega(3));
fprintf('  |omega| = %.6f rad/s\n', norm(omega));
fprintf('  Rotation axis: [%.6f, %.6f, %.6f]\n\n', omega_axis(1), omega_axis(2), omega_axis(3));

%% Define Test Vectors in Body Frame

% Choose 3 orthogonal unit vectors for testing
v1_body = [1; 0; 0]; % Body X-axis
v2_body = [0; 1; 0]; % Body Y-axis
v3_body = [0; 0; 1]; % Body Z-axis

% Additional test: arbitrary vector
v4_body = [1; 1; 1] / sqrt(3); % Normalized diagonal vector

fprintf('Test vectors in body frame:\n');
fprintf('  v1 = [%.2f, %.2f, %.2f] (Body X-axis)\n', v1_body(1), v1_body(2), v1_body(3));
fprintf('  v2 = [%.2f, %.2f, %.2f] (Body Y-axis)\n', v2_body(1), v2_body(2), v2_body(3));
fprintf('  v3 = [%.2f, %.2f, %.2f] (Body Z-axis)\n', v3_body(1), v3_body(2), v3_body(3));
fprintf('  v4 = [%.4f, %.4f, %.4f] (Diagonal)\n\n', v4_body(1), v4_body(2), v4_body(3));

%% Initial Quaternion (Identity - Body frame aligned with inertial frame)

q0 = [0; 0; 0; 1]; % Identity quaternion (no initial rotation)

fprintf('Initial quaternion (identity):\n');
fprintf('  q0 = [%.2f, %.2f, %.2f, %.2f]\n\n', q0(1), q0(2), q0(3), q0(4));

%% Setup Physical Parameters (Minimal for this test)

params.mass = 100; % kg
params.mu = 3.986004418e14; % m^3/s^2 (not used, but required by template)
params.R_e = 6378137; % m (not used)

% Moment of inertia (symmetric for simple test)
params.I_CB = 100 * eye(3); % kg*m^2

%% Simulation Parameters

test_duration = 2 * pi / omega_magnitude; % One full rotation period
t_span = linspace(0, test_duration, 100); % Time vector

%% Initialize Storage Arrays

q_numeric = zeros(length(t_span), 4);
q_analytic = zeros(length(t_span), 4);

% Store rotated vectors (for v1 as example)
v1_numeric = zeros(length(t_span), 3);
v1_analytic = zeros(length(t_span), 3);

%% Time-stepping Integration

fprintf('Running numerical integration...\n');

% Use a simple fixed-step integration for clarity
dt = t_span(2) - t_span(1);
q_current = q0;

for i = 1:length(t_span)
    t = t_span(i);

    % Store current quaternion (numeric)
    q_numeric(i, :) = q_current';

    % Compute analytic quaternion at this time
    theta = omega_magnitude * t; % Total rotation angle
    q_analytic(i, :) = axis_angle_to_quaternion(omega_axis, theta)';

    % Rotate v1 using numeric quaternion
    v1_numeric(i, :) = quat_rotate_vector(q_current, v1_body)';

    % Rotate v1 using analytic rotation matrix
    R_analytic = axis_angle_to_rotation_matrix(omega_axis, theta);
    v1_analytic(i, :) = (R_analytic * v1_body)';

    % Integrate quaternion for next step (if not last step)
    if i < length(t_span)
        % Quaternion derivative: dq/dt = 0.5 * Omega(omega) * q
        Omega = [0,        -omega(3),  omega(2),  omega(1);
                 omega(3),  0,        -omega(1),  omega(2);
                -omega(2),  omega(1),  0,         omega(3);
                -omega(1), -omega(2), -omega(3),  0];

        dqdt = 0.5 * Omega * q_current;

        % Euler integration (simple for demonstration)
        q_current = q_current + dqdt * dt;

        % Renormalize quaternion to prevent drift
        q_current = q_current / norm(q_current);
    end
end

fprintf('Integration complete.\n\n');

%% Error Analysis

% Quaternion error
q_error = q_numeric - q_analytic;
q_error_norm = vecnorm(q_error, 2, 2);

% Vector rotation error
v1_error = v1_numeric - v1_analytic;
v1_error_norm = vecnorm(v1_error, 2, 2);

fprintf('=== Error Analysis ===\n');
fprintf('Maximum quaternion error: %.3e\n', max(q_error_norm));
fprintf('Mean quaternion error: %.3e\n', mean(q_error_norm));
fprintf('Maximum vector rotation error: %.3e\n', max(v1_error_norm));
fprintf('Mean vector rotation error: %.3e\n\n', mean(v1_error_norm));

% Check if errors are within tolerance
tolerance_quaternion = 1e-3;
tolerance_vector = 1e-6;

fprintf('=== Verification Results ===\n');
if max(q_error_norm) < tolerance_quaternion
    fprintf('PASS: Quaternion integration matches analytic solution\n');
    fprintf('      (max error: %.3e < tolerance: %.3e)\n', max(q_error_norm), tolerance_quaternion);
else
    fprintf('FAIL: Quaternion integration error exceeds tolerance\n');
    fprintf('      (max error: %.3e > tolerance: %.3e)\n', max(q_error_norm), tolerance_quaternion);
end

if max(v1_error_norm) < tolerance_vector
    fprintf('PASS: Vector rotation matches R(omega*t)\n');
    fprintf('      (max error: %.3e < tolerance: %.3e)\n', max(v1_error_norm), tolerance_vector);
else
    fprintf('FAIL: Vector rotation error exceeds tolerance\n');
    fprintf('      (max error: %.3e > tolerance: %.3e)\n', max(v1_error_norm), tolerance_vector);
end

%% Visualization

% Plot 1: Quaternion Components vs Time
figure('Position', [100, 100, 1200, 800]);

subplot(2, 2, 1);
plot(t_span, q_numeric(:, 1), 'r-', 'LineWidth', 1.5); hold on;
plot(t_span, q_analytic(:, 1), 'r--', 'LineWidth', 2);
ylabel('q_x', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Component q_x', 'FontSize', 13);
legend('Numeric', 'Analytic', 'Location', 'best');
grid on;

subplot(2, 2, 2);
plot(t_span, q_numeric(:, 2), 'g-', 'LineWidth', 1.5); hold on;
plot(t_span, q_analytic(:, 2), 'g--', 'LineWidth', 2);
ylabel('q_y', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Component q_y', 'FontSize', 13);
legend('Numeric', 'Analytic', 'Location', 'best');
grid on;

subplot(2, 2, 3);
plot(t_span, q_numeric(:, 3), 'b-', 'LineWidth', 1.5); hold on;
plot(t_span, q_analytic(:, 3), 'b--', 'LineWidth', 2);
ylabel('q_z', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Component q_z', 'FontSize', 13);
legend('Numeric', 'Analytic', 'Location', 'best');
grid on;

subplot(2, 2, 4);
plot(t_span, q_numeric(:, 4), 'm-', 'LineWidth', 1.5); hold on;
plot(t_span, q_analytic(:, 4), 'm--', 'LineWidth', 2);
ylabel('q_w', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Component q_w (scalar)', 'FontSize', 13);
legend('Numeric', 'Analytic', 'Location', 'best');
grid on;

sgtitle('Quaternion Components: Numeric vs Analytic', 'FontSize', 15, 'FontWeight', 'bold');

% Plot 2: Error in Quaternion Components
figure('Position', [150, 150, 1200, 600]);

subplot(1, 2, 1);
plot(t_span, q_error(:, 1), 'r-', 'LineWidth', 1.5); hold on;
plot(t_span, q_error(:, 2), 'g-', 'LineWidth', 1.5);
plot(t_span, q_error(:, 3), 'b-', 'LineWidth', 1.5);
plot(t_span, q_error(:, 4), 'm-', 'LineWidth', 1.5);
plot(t_span, zeros(size(t_span)), 'k--', 'LineWidth', 1);
ylabel('Error', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Component Errors', 'FontSize', 13);
legend('q_x error', 'q_y error', 'q_z error', 'q_w error', 'Zero', 'Location', 'best');
grid on;

subplot(1, 2, 2);
plot(t_span, q_error_norm, 'k-', 'LineWidth', 2);
ylabel('2-Norm Error', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Quaternion Error Magnitude', 'FontSize', 13);
grid on;

sgtitle('Quaternion Integration Error', 'FontSize', 15, 'FontWeight', 'bold');

% Plot 3: 3D Trajectory of Rotated Vector v1
figure('Position', [200, 200, 800, 800]);
plot3(v1_numeric(:, 1), v1_numeric(:, 2), v1_numeric(:, 3), 'b-', 'LineWidth', 2); hold on;
plot3(v1_analytic(:, 1), v1_analytic(:, 2), v1_analytic(:, 3), 'r--', 'LineWidth', 2);
plot3(v1_numeric(1, 1), v1_numeric(1, 2), v1_numeric(1, 3), 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(v1_numeric(end, 1), v1_numeric(end, 2), v1_numeric(end, 3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% Draw rotation axis
quiver3(0, 0, 0, omega_axis(1), omega_axis(2), omega_axis(3), 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);

% Draw unit sphere
[X, Y, Z] = sphere(30);
surf(X, Y, Z, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

xlabel('X', 'FontSize', 12);
ylabel('Y', 'FontSize', 12);
zlabel('Z', 'FontSize', 12);
title('Trajectory of Rotated Vector v1 on Unit Sphere', 'FontSize', 14);
legend('Numeric', 'Analytic', 'Start', 'End', 'Rotation Axis', 'Location', 'best');
grid on;
axis equal;
view(3);

%% Helper Functions

function q = axis_angle_to_quaternion(axis, angle)
    % Convert axis-angle representation to quaternion (scalar-last)
    % q = [qx, qy, qz, qw]
    axis = axis / norm(axis); % Ensure unit vector
    half_angle = angle / 2;
    q = [axis * sin(half_angle); cos(half_angle)];
end

function R = axis_angle_to_rotation_matrix(axis, angle)
    % Rodrigues' rotation formula
    axis = axis / norm(axis);
    K = [0,        -axis(3),  axis(2);
         axis(3),   0,       -axis(1);
        -axis(2),   axis(1),  0];
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K * K);
end

function v_rotated = quat_rotate_vector(q, v)
    % Rotate vector v by quaternion q (scalar-last convention)
    % v_rotated = q * v * q^*

    % Convert vector to quaternion form [v; 0]
    v_quat = [v; 0];

    % Quaternion multiplication: q * v_quat
    q_conj = [-q(1:3); q(4)];

    % First multiplication: q * v_quat
    temp = quatmultiply_helper(q, v_quat);

    % Second multiplication: temp * q_conj
    result = quatmultiply_helper(temp, q_conj);

    % Extract vector part
    v_rotated = result(1:3);
end

function q_result = quatmultiply_helper(q1, q2)
    % Quaternion multiplication (scalar-last convention)
    v1 = q1(1:3);
    s1 = q1(4);
    v2 = q2(1:3);
    s2 = q2(4);

    v_result = s1*v2 + s2*v1 + cross(v1, v2);
    s_result = s1*s2 - dot(v1, v2);

    q_result = [v_result; s_result];
end
