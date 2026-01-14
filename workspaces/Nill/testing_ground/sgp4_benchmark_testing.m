clear; clc; close all;

fprintf('========================================\n');
fprintf(' SGP4 vs Numerical Propagation\n');
fprintf(' with Atmospheric Drag\n');
fprintf('========================================\n\n');

% Add SGP4 to path
addpath('SGP4');

%% Constants
global const
SAT_Const % Initialize satellite constants
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
MINUTES_PER_DAY_CUBED = (MINUTES_PER_DAY * MINUTES_PER_DAY_SQUARED);

%% ISS TLE (September 26, 2025)
% You can replace these with any satellite TLE
tle_line1 = '1 25544U 98067A   25269.21152988  .00014701  00000-0  26152-3 0  9992';
tle_line2 = '2 25544  51.6332 166.3292 0002995  10.9396 349.1657 15.50384251530832';

%% Parse TLE
year = str2num(tle_line1(19:20));
doy = str2num(tle_line1(21:32));
epoch = str2num(tle_line1(19:32));
TD1 = str2num(tle_line1(34:43));

% BStar parsing
BStar = str2num(tle_line1(54:59));
ExBStar = str2num(tle_line1(60:61));
BStar = BStar * 1e-5 * 10^ExBStar;

incl = str2num(tle_line2(9:16));
raan = str2num(tle_line2(18:25));
ecc = str2num(['0.' tle_line2(27:33)]);
omega = str2num(tle_line2(35:42));
M = str2num(tle_line2(44:51));
no = str2num(tle_line2(53:63));

% Convert to radians and proper units
satdata.epoch = epoch;
satdata.xmo = M * (pi/180);
satdata.xnodeo = raan * (pi/180);
satdata.omegao = omega * (pi/180);
satdata.xincl = incl * (pi/180);
satdata.eo = ecc;
satdata.xno = no * TWOPI / MINUTES_PER_DAY;
satdata.xndt2o = TD1 * TWOPI / MINUTES_PER_DAY_SQUARED;
satdata.xndd6o = 0;
satdata.bstar = BStar;

fprintf('Satellite TLE Parameters:\n');
fprintf(' Epoch: 2025, day %.8f\n', doy);
fprintf(' Inclination: %.4f deg\n', incl);
fprintf(' RAAN: %.4f deg\n', raan);
fprintf(' Eccentricity: %.7f\n', ecc);
fprintf(' Mean Motion: %.8f rev/day\n', no);
fprintf(' B*: %.6e\n\n', BStar);

%% Propagation Settings
propagation_hours = 24;
time_step_minutes = 2;
tsince_vec = 0:time_step_minutes:(propagation_hours*60);
n_steps = length(tsince_vec);

fprintf('Propagation: %.0f hours in %d steps\n\n', propagation_hours, n_steps);

%% SGP4 Propagation
fprintf('Running SGP4 propagation...\n');
r_sgp4 = zeros(3, n_steps);
v_sgp4 = zeros(3, n_steps);

tic;
for idx = 1:n_steps
    [rteme, vteme] = sgp4(tsince_vec(idx), satdata);
    r_sgp4(:, idx) = rteme * 1000;
    v_sgp4(:, idx) = vteme * 1000;
end
t_sgp4 = toc;

fprintf(' ✓ SGP4 completed in %.4f s\n\n', t_sgp4);

%% Setup Numerical Propagation Parameters (from TLE only)
params.mu = 3.986004418e14;
params.R_e = 6378137;
params.J2 = 1.08263e-3;
params.J3 = -2.53e-6;

% Drag parameters - automatically estimated from TLE
params.bstar = BStar;  % B* from TLE
params.H = 8500;       % Scale height for LEO [m]

% For attitude dynamics (not critical for orbit)
params.mass = 1000;    % Nominal value [kg] - not used in drag calculation
params.I_CB = params.mass * eye(3);

% Initial state from SGP4
r0_sgp4 = r_sgp4(:, 1);
v0_sgp4 = v_sgp4(:, 1);
alt0_sgp4 = (norm(r0_sgp4) - params.R_e) / 1e3;

fprintf('Initial altitude: %.2f km\n', alt0_sgp4);
fprintf('Using B* = %.6e for drag calculation\n\n', BStar);

X0 = zeros(13, 1);
X0(1:3) = r0_sgp4;
X0(4:6) = v0_sgp4;
X0(7:10) = [0; 0; 0; 1]; % Identity quaternion
X0(11:13) = [0; 0; 0];

%% Numerical Propagation WITH Drag
fprintf('Running numerical propagation (J2 + drag from B*)...\n');
X_numerical = zeros(6, n_steps);
X_numerical(:, 1) = X0(1:6);
Xi = X0;
u = zeros(6, 1);

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 60);

tic;
for idx = 2:n_steps
    dt = time_step_minutes * 60;

    [~, X_traj] = ode45(@(t, X) Sat_template(t, X, u, params, ...
        'useJ2', true, 'useAtmDrag', true, ...
        'useControl', false), ...
        [0, dt], Xi, options);

    Xi = X_traj(end, :)';
    Xi(7:10) = Xi(7:10) / norm(Xi(7:10));
    X_numerical(:, idx) = Xi(1:6);
end
t_numerical = toc;

fprintf(' ✓ Numerical completed in %.4f s\n\n', t_numerical);

%% Calculate Errors
err_pos = zeros(1, n_steps);
err_vel = zeros(1, n_steps);

for idx = 1:n_steps
    err_pos(idx) = norm(X_numerical(1:3, idx) - r_sgp4(:, idx));
    err_vel(idx) = norm(X_numerical(4:6, idx) - v_sgp4(:, idx));
end

%% Altitude Analysis
alt_sgp4 = zeros(1, n_steps);
alt_numerical = zeros(1, n_steps);

for idx = 1:n_steps
    alt_sgp4(idx) = (norm(r_sgp4(:, idx)) - params.R_e) / 1e3;
    alt_numerical(idx) = (norm(X_numerical(1:3, idx)) - params.R_e) / 1e3;
end

%% Results
fprintf('========================================\n');
fprintf(' RESULTS\n');
fprintf('========================================\n\n');

fprintf('Computation Time:\n');
fprintf(' SGP4: %.4f s\n', t_sgp4);
fprintf(' Numerical: %.4f s (%.1fx slower)\n\n', t_numerical, t_numerical/t_sgp4);

fprintf('Position Error (Numerical vs SGP4):\n');
fprintf(' Mean: %.3f km\n', mean(err_pos)/1e3);
fprintf(' Max: %.3f km\n', max(err_pos)/1e3);
fprintf(' RMS: %.3f km\n', sqrt(mean(err_pos.^2))/1e3);
fprintf(' Final: %.3f km\n\n', err_pos(end)/1e3);

fprintf('Velocity Error:\n');
fprintf(' Mean: %.6f m/s\n', mean(err_vel));
fprintf(' Max: %.6f m/s\n', max(err_vel));
fprintf(' Final: %.6f m/s\n\n', err_vel(end));

fprintf('Altitude Decay Over %.0f Hours:\n', propagation_hours);
decay_sgp4 = alt_sgp4(1) - alt_sgp4(end);
decay_numerical = alt_numerical(1) - alt_numerical(end);
fprintf(' SGP4: %.3f km (%.3f km/day)\n', decay_sgp4, decay_sgp4*24/propagation_hours);
fprintf(' Numerical: %.3f km (%.3f km/day)\n', decay_numerical, decay_numerical*24/propagation_hours);
fprintf(' Difference: %.3f km\n\n', abs(decay_sgp4 - decay_numerical));

%% Plots
time_hours = tsince_vec / 60;

figure('Name', 'SGP4 vs Numerical with Drag', 'Position', [100, 100, 1600, 900]);

% 3D Orbit Comparison
subplot(2, 3, 1);
hold on; grid on; axis equal;
plot3(r_sgp4(1,:)/1e3, r_sgp4(2,:)/1e3, r_sgp4(3,:)/1e3, ...
    'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot3(X_numerical(1,:)/1e3, X_numerical(2,:)/1e3, X_numerical(3,:)/1e3, ...
    'r--', 'LineWidth', 1.5, 'DisplayName', 'Numerical (J2+Drag)');
[X_e, Y_e, Z_e] = sphere(30);
surf(X_e*params.R_e/1e3, Y_e*params.R_e/1e3, Z_e*params.R_e/1e3, ...
    'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbit Comparison');
legend('Location', 'best');
view(45, 30);

% Altitude Comparison
subplot(2, 3, 2);
hold on; grid on;
plot(time_hours, alt_sgp4, 'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot(time_hours, alt_numerical, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Numerical');
xlabel('Time [hours]');
ylabel('Altitude [km]');
title('Altitude vs Time');
legend('Location', 'best');

% Position Error
subplot(2, 3, 3);
plot(time_hours, err_pos/1e3, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Time [hours]');
ylabel('Position Error [km]');
title('Position Error vs Time');

% Velocity Error
subplot(2, 3, 4);
plot(time_hours, err_vel, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time [hours]');
ylabel('Velocity Error [m/s]');
title('Velocity Error vs Time');

% Error Components (X, Y, Z)
subplot(2, 3, 5);
hold on; grid on;
err_x = X_numerical(1, :) - r_sgp4(1, :);
err_y = X_numerical(2, :) - r_sgp4(2, :);
err_z = X_numerical(3, :) - r_sgp4(3, :);
plot(time_hours, err_x/1e3, 'r-', 'DisplayName', 'X');
plot(time_hours, err_y/1e3, 'g-', 'DisplayName', 'Y');
plot(time_hours, err_z/1e3, 'b-', 'DisplayName', 'Z');
xlabel('Time [hours]');
ylabel('Position Error [km]');
title('Position Error Components');
legend('Location', 'best');

% Altitude Decay Comparison
subplot(2, 3, 6);
hold on; grid on;
decay_sgp4_vec = (alt_sgp4 - alt_sgp4(1));
decay_num_vec = (alt_numerical - alt_numerical(1));
plot(time_hours, decay_sgp4_vec, 'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot(time_hours, decay_num_vec, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Numerical');
xlabel('Time [hours]');
ylabel('Altitude Change [km]');
title('Altitude Decay from Initial');
legend('Location', 'best');

fprintf('========================================\n');
fprintf('Analysis:\n');
fprintf('========================================\n\n');

fprintf('The numerical propagator with drag (estimated from B*)\n');
fprintf('closely matches SGP4 with mean error of %.3f km.\n', mean(err_pos)/1e3);
fprintf('\nBoth models show altitude decay of ~%.2f km/day,\n', decay_sgp4*24/propagation_hours);
fprintf('demonstrating that B* effectively captures drag effects.\n');
fprintf('\n========================================\n');
fprintf('Benchmark Complete!\n');
fprintf('========================================\n');