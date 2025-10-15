%% Test Script: Verify Constant Angular Velocity in Torque-Free Environment

% Hypothesis: With only gravity (no J2, drag, or torques), omega should remain constant

clear; clc; close all;

%% Setup Physical Parameters

params.mass = 100; % kg

params.mu = 3.986004418e14; % m^3/s^2

params.R_e = 6378137; % m

% Moment of inertia tensor (asymmetric for realism)
% params.I_CB = [50, 0, 0;
%                 0, 75, 0;
%                 0, 0, 100]; % kg*m^2

%% Initial Conditions for Circular Orbit

altitude = 500e3; % 500 km altitude

r0 = params.R_e + altitude; % orbital radius

% Circular orbit velocity
v_circ = sqrt(params.mu / r0);

% Initial position (start on equatorial plane)
r_vec0 = [r0; 0; 0]; % m

% Initial velocity (circular orbit in XY plane)
v_vec0 = [0; v_circ; 0]; % m/s

%% Random Initial Angular Velocity

% Generate random angular velocity vector (rad/s)
omega0 = randn(3, 1) * 0.05; % Random omega with magnitude 0.05 rad/s

% Initial quaternion (random orientation)
% Using Shoemake's method for uniform random quaternion
u1 = rand();
u2 = rand();
u3 = rand();

q0 = [sqrt(1-u1)*sin(2*pi*u2);   % qx
      sqrt(1-u1)*cos(2*pi*u2);   % qy
      sqrt(u1)*sin(2*pi*u3);     % qz
      sqrt(u1)*cos(2*pi*u3)];    % qw (scalar)

%% Assemble Initial State Vector

X0 = [r_vec0; v_vec0; q0; omega0];

%% Simulation Time

orbital_period = 2 * pi * sqrt(r0^3 / params.mu); % seconds

t_span = [0, 3 * orbital_period]; % Simulate for 3 orbits

%% Desired State (Set to Initial State for Torque-Free Test)

desired_state = X0; % No desired state change

%% ODE Options

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);

%% Run Simulation with Only Gravity (No J2, No Drag, No Control)

fprintf('Running simulation with only gravitational force...\n');
fprintf('Initial omega: [%.6f, %.6f, %.6f] rad/s\n', omega0(1), omega0(2), omega0(3));

[t, X] = ode45(@(t, X) Sat_template(t, X, desired_state, params, ...
                                     'useJ2', false, ...
                                     'useAtmDrag', false, ...
                                     'useControl', false, ...
                                     'useSRP', false), ...
               t_span, X0, options);

%% Extract Angular Velocity Components

omega_x = X(:, 11);
omega_y = X(:, 12);
omega_z = X(:, 13);

%% Extract Position for Orbit Plot

r_x = X(:, 1);
r_y = X(:, 2);
r_z = X(:, 3);

%% Calculate Deviation from Initial Values

delta_omega_x = omega_x - omega0(1);
delta_omega_y = omega_y - omega0(2);
delta_omega_z = omega_z - omega0(3);

% Calculate magnitude of angular velocity
omega_mag = sqrt(omega_x.^2 + omega_y.^2 + omega_z.^2);

%% Statistical Analysis

fprintf('\n=== Angular Velocity Stability Analysis ===\n');
fprintf('Simulation duration: %.2f hours (%.2f orbits)\n', ...
        t(end)/3600, t(end)/orbital_period);

fprintf('\nInitial omega: [%.8f, %.8f, %.8f] rad/s\n', omega0(1), omega0(2), omega0(3));
fprintf('Final omega:   [%.8f, %.8f, %.8f] rad/s\n', omega_x(end), omega_y(end), omega_z(end));

fprintf('\nMaximum deviation from initial:\n');
fprintf('  omega_x: %.3e rad/s\n', max(abs(delta_omega_x)));
fprintf('  omega_y: %.3e rad/s\n', max(abs(delta_omega_y)));
fprintf('  omega_z: %.3e rad/s\n', max(abs(delta_omega_z)));

fprintf('\nStandard deviation:\n');
fprintf('  omega_x: %.3e rad/s\n', std(omega_x));
fprintf('  omega_y: %.3e rad/s\n', std(omega_y));
fprintf('  omega_z: %.3e rad/s\n', std(omega_z));

fprintf('\nInitial omega magnitude: %.8f rad/s\n', norm(omega0));
fprintf('Final omega magnitude:   %.8f rad/s\n', omega_mag(end));
fprintf('Change in magnitude:     %.3e rad/s\n', omega_mag(end) - norm(omega0));

%% 3D Orbit Plot

figure;
plot3(r_x/1e3, r_y/1e3, r_z/1e3, 'b-', 'LineWidth', 1.5);
hold on;

% Plot Earth
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth*params.R_e/1e3, Y_earth*params.R_e/1e3, Z_earth*params.R_e/1e3, ...
     'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Mark start and end points
plot3(r_x(1)/1e3, r_y(1)/1e3, r_z(1)/1e3, 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot3(r_x(end)/1e3, r_y(end)/1e3, r_z(end)/1e3, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('X (km)', 'FontSize', 12);
ylabel('Y (km)', 'FontSize', 12);
zlabel('Z (km)', 'FontSize', 12);
title('3D Orbital Trajectory', 'FontSize', 14);
legend('Orbit', 'Earth', 'Start', 'End', 'Location', 'best');
grid on;
axis equal;
view(3);

%% Improved Separate Plots for Each Omega Component

% Convert time to hours for better readability
t_hours = t / 3600;

% Figure: Omega_x vs Time
figure;
plot(t_hours, omega_x, 'r-', 'LineWidth', 1.5); hold on;
plot(t_hours, omega0(1) * ones(size(t_hours)), 'r--', 'LineWidth', 2);
ylabel('\omega_x (rad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Omega_x vs Time', 'FontSize', 14);
legend('\omega_x(t)', 'Initial \omega_x', 'Location', 'best');
grid on;

% Figure: Omega_y vs Time
figure;
plot(t_hours, omega_y, 'g-', 'LineWidth', 1.5); hold on;
plot(t_hours, omega0(2) * ones(size(t_hours)), 'g--', 'LineWidth', 2);
ylabel('\omega_y (rad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Omega_y vs Time', 'FontSize', 14);
legend('\omega_y(t)', 'Initial \omega_y', 'Location', 'best');
grid on;

% Figure: Omega_z vs Time
figure;
plot(t_hours, omega_z, 'b-', 'LineWidth', 1.5); hold on;
plot(t_hours, omega0(3) * ones(size(t_hours)), 'b--', 'LineWidth', 2);
ylabel('\omega_z (rad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Omega_z vs Time', 'FontSize', 14);
legend('\omega_z(t)', 'Initial \omega_z', 'Location', 'best');
grid on;

% Figure: Delta Omega_x vs Time
figure;
plot(t_hours, delta_omega_x * 1e6, 'r-', 'LineWidth', 1.5); hold on;
plot(t_hours, zeros(size(t_hours)), 'r--', 'LineWidth', 2);
ylabel('\Delta\omega_x (\murad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Deviation in Omega_x vs Time', 'FontSize', 14);
legend('\Delta\omega_x(t)', 'Zero Reference', 'Location', 'best');
grid on;

% Figure: Delta Omega_y vs Time
figure;
plot(t_hours, delta_omega_y * 1e6, 'g-', 'LineWidth', 1.5); hold on;
plot(t_hours, zeros(size(t_hours)), 'g--', 'LineWidth', 2);
ylabel('\Delta\omega_y (\murad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Deviation in Omega_y vs Time', 'FontSize', 14);
legend('\Delta\omega_y(t)', 'Zero Reference', 'Location', 'best');
grid on;

% Figure: Delta Omega_z vs Time
figure;
plot(t_hours, delta_omega_z * 1e6, 'b-', 'LineWidth', 1.5); hold on;
plot(t_hours, zeros(size(t_hours)), 'b--', 'LineWidth', 2);
ylabel('\Delta\omega_z (\murad/s)', 'FontSize', 12);
xlabel('Time (hours)', 'FontSize', 12);
title('Deviation in Omega_z vs Time', 'FontSize', 14);
legend('\Delta\omega_z(t)', 'Zero Reference', 'Location', 'best');
grid on;

%% Verification Result

tolerance = 1e-8; % rad/s
max_deviation = max([max(abs(delta_omega_x)), ...
                     max(abs(delta_omega_y)), ...
                     max(abs(delta_omega_z))]);

fprintf('\n=== VERIFICATION RESULT ===\n');
if max_deviation < tolerance
    fprintf('PASS: Angular velocity remained constant (max deviation: %.3e rad/s)\n', max_deviation);
else
    fprintf('FAIL: Angular velocity changed significantly (max deviation: %.3e rad/s)\n', max_deviation);
    fprintf('This suggests an error in the quaternion kinematics implementation.\n');
end
fprintf('Tolerance: %.3e rad/s\n', tolerance);
