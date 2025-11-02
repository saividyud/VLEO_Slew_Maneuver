clear
clc
close all

%% Simplified Benchmark Testing: Numerical 2-Body vs Analytical
% This script compares the Sat_template numerical propagator (2-body only)
% with the analytical Kepler propagator to assess numerical accuracy

fprintf('================================================\n');
fprintf('  2-Body Propagator Benchmark Test\n');
fprintf('  Numerical vs Analytical Comparison\n');
fprintf('================================================\n\n');

%% Test Configuration
% Single test case: LEO Circular Orbit
altitude = 400e3;           % Altitude [m]
eccentricity = 0.01;         % Circular orbit
propagation_time = 3600 * 24;  % 24 hours [s]
time_step = 60;             % 1 minute [s]
time_vec = 0:time_step:propagation_time;

fprintf('Test Configuration:\n');
fprintf('  Orbit: LEO Circular\n');
fprintf('  Altitude: %.0f km\n', altitude/1e3);
fprintf('  Eccentricity: %.3f\n', eccentricity);
fprintf('  Propagation time: %.1f hours\n', propagation_time/3600);
fprintf('  Time step: %.0f seconds\n\n', time_step);

%% Physical Parameters
params.mu = 3.986004418e14;     % Earth gravitational parameter [m^3/s^2]
params.R_e = 6378137;           % Earth equatorial radius [m]
params.J2 = 1.08263e-3;         % J2 coefficient (not used in 2-body)
params.J3 = -2.53e-6;           % J3 coefficient (not used in 2-body)
params.mass = 100;              % Satellite mass [kg]
params.I_CB = params.mass * diag([1, 1, 1]);  % Moment of inertia [kg*m^2]

%% Generate Initial State
r_perigee = params.R_e + altitude;
a = r_perigee / (1 - eccentricity);  % Semi-major axis
v_perigee = sqrt(params.mu * (1 + eccentricity) / a);

X0 = zeros(13, 1);
X0(1:3) = [r_perigee; 0; 0];        % Position at perigee [m]
X0(4:6) = [0; v_perigee; 0];        % Velocity at perigee [m/s]
X0(7:10) = [0; 0; 0; 1];            % Quaternion (identity)
X0(11:13) = [0; 0; 0];              % Angular velocity [rad/s]

fprintf('Initial State:\n');
fprintf('  Position: [%.3f, %.3f, %.3f] km\n', X0(1)/1e3, X0(2)/1e3, X0(3)/1e3);
fprintf('  Velocity: [%.3f, %.3f, %.3f] m/s\n', X0(4), X0(5), X0(6));
fprintf('  Orbital period: %.2f minutes\n\n', 2*pi*sqrt(a^3/params.mu)/60);

%% Analytical Propagation (Baseline)
fprintf('Running analytical propagation...\n');
tic;
X_analytical = Analytical_Propagator(X0, time_vec, params);
t_analytical = toc;
fprintf('  ✓ Analytical propagation completed in %.4f seconds\n\n', t_analytical);

%% Numerical Propagation (Two-Body Only)
fprintf('Running numerical propagation (2-body only)...\n');
X_numerical = zeros(6, length(time_vec));
X_numerical(:, 1) = X0(1:6);

Xi = X0;
u = zeros(6, 1);  % No control input

tic;
for i = 2:length(time_vec)
    dt = time_step;

    % Integrate using ode45 with 2-body dynamics only
    [~, X_traj] = ode45(@(t, X) Sat_template(t, X, u, params, ...
                                 'useJ2', false, 'useAtmDrag', false, ...
                                 'useControl', false), ...
                        [0, dt], Xi);

    Xi = X_traj(end, :)';
    Xi(7:10) = Xi(7:10) / norm(Xi(7:10));  % Normalize quaternion
    X_numerical(:, i) = Xi(1:6);

end
t_numerical = toc;
fprintf('  ✓ Numerical propagation completed in %.4f seconds\n\n', t_numerical);

%% Calculate Errors
err_pos = zeros(1, length(time_vec));
err_vel = zeros(1, length(time_vec));

for i = 1:length(time_vec)
    err_pos(i) = norm(X_numerical(1:3, i) - X_analytical(1:3, i));
    err_vel(i) = norm(X_numerical(4:6, i) - X_analytical(4:6, i));
end

%% Print Statistics
fprintf('================================================\n');
fprintf('  RESULTS SUMMARY\n');
fprintf('================================================\n\n');

fprintf('Computation Time:\n');
fprintf('  Analytical: %.4f seconds\n', t_analytical);
fprintf('  Numerical:  %.4f seconds\n', t_numerical);
fprintf('  Speedup:    %.1fx (Analytical is %.1fx faster)\n\n', ...
        t_numerical/t_analytical, t_numerical/t_analytical);

fprintf('Position Error Statistics (Numerical vs Analytical):\n');
fprintf('  Mean error:   %.6f m  (%.6e m)\n', mean(err_pos), mean(err_pos));
fprintf('  Max error:    %.6f m  (%.6e m)\n', max(err_pos), max(err_pos));
fprintf('  Final error:  %.6f m  (%.6e m)\n\n', err_pos(end), err_pos(end));

fprintf('Velocity Error Statistics (Numerical vs Analytical):\n');
fprintf('  Mean error:   %.9f m/s  (%.6e m/s)\n', mean(err_vel), mean(err_vel));
fprintf('  Max error:    %.9f m/s  (%.6e m/s)\n', max(err_vel), max(err_vel));
fprintf('  Final error:  %.9f m/s  (%.6e m/s)\n\n', err_vel(end), err_vel(end));

fprintf('Relative Errors:\n');
r_mean = mean(sqrt(sum(X_analytical(1:3,:).^2, 1)));
v_mean = mean(sqrt(sum(X_analytical(4:6,:).^2, 1)));
fprintf('  Position (relative): %.6e\n', mean(err_pos)/r_mean);
fprintf('  Velocity (relative): %.6e\n\n', mean(err_vel)/v_mean);

%% Generate Plots
fprintf('Generating plots...\n');

figure('Name', '2-Body vs Analytical Comparison', 'Position', [100, 100, 1400, 900]);

% Plot 1: 3D Orbit Comparison
subplot(2, 3, 1);
hold on; grid on; axis equal;
plot3(X_analytical(1,:)/1e3, X_analytical(2,:)/1e3, X_analytical(3,:)/1e3, ...
      'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
plot3(X_numerical(1,:)/1e3, X_numerical(2,:)/1e3, X_numerical(3,:)/1e3, ...
      'r--', 'LineWidth', 1.5, 'DisplayName', 'Numerical 2-body');

% Draw Earth
[X_e, Y_e, Z_e] = sphere(30);
surf(X_e*params.R_e/1e3, Y_e*params.R_e/1e3, Z_e*params.R_e/1e3, ...
     'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Orbit Comparison');
legend('Location', 'best');
view(45, 30);

% Plot 2: Position Error vs Time
subplot(2, 3, 2);
plot(time_vec/3600, err_pos, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time [hours]');
ylabel('Position Error [m]');
title('Position Error: Numerical vs Analytical');

% Plot 3: Velocity Error vs Time
subplot(2, 3, 3);
plot(time_vec/3600, err_vel*1000, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time [hours]');
ylabel('Velocity Error [mm/s]');
title('Velocity Error: Numerical vs Analytical');

% Plot 4: Position Error Components
subplot(2, 3, 4);
hold on; grid on;
for i = 1:length(time_vec)
    err_x(i) = X_numerical(1,i) - X_analytical(1,i);
    err_y(i) = X_numerical(2,i) - X_analytical(2,i);
    err_z(i) = X_numerical(3,i) - X_analytical(3,i);
end
plot(time_vec/3600, err_x, 'r-', 'DisplayName', 'X Error');
plot(time_vec/3600, err_y, 'g-', 'DisplayName', 'Y Error');
plot(time_vec/3600, err_z, 'b-', 'DisplayName', 'Z Error');
xlabel('Time [hours]');
ylabel('Position Error [m]');
title('Position Error Components');
legend('Location', 'best');

% Plot 5: Velocity Error Components
subplot(2, 3, 5);
hold on; grid on;
for i = 1:length(time_vec)
    err_vx(i) = X_numerical(4,i) - X_analytical(4,i);
    err_vy(i) = X_numerical(5,i) - X_analytical(5,i);
    err_vz(i) = X_numerical(6,i) - X_analytical(6,i);
end
plot(time_vec/3600, err_vx*1000, 'r-', 'DisplayName', 'Vx Error');
plot(time_vec/3600, err_vy*1000, 'g-', 'DisplayName', 'Vy Error');
plot(time_vec/3600, err_vz*1000, 'b-', 'DisplayName', 'Vz Error');
xlabel('Time [hours]');
ylabel('Velocity Error [mm/s]');
title('Velocity Error Components');
legend('Location', 'best');

% Plot 6: Log-scale Error Growth
subplot(2, 3, 6);
semilogy(time_vec/3600, err_pos, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Position');
hold on; grid on;
semilogy(time_vec/3600, err_vel*1000, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Velocity (mm/s)');
xlabel('Time [hours]');
ylabel('Error [m or mm/s]');
title('Error Growth (Log Scale)');
legend('Location', 'best');

fprintf('\n================================================\n');
fprintf('Benchmark testing complete!\n');
fprintf('================================================\n');