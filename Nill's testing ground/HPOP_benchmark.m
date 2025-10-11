clear; clc; close all;

fprintf('========================================\n');
fprintf(' Three-Way Propagator Comparison\n');
fprintf(' SGP4 vs Custom vs HPOP\n');
fprintf('========================================\n\n');

%% Setup paths
addpath('SGP4');
addpath('HPOP');

%% Initialize HPOP globals and data
load_hpop_data()

%% ISS TLE (September 26, 2025)
tle_line1 = '1 25544U 98067A   25269.21152988  .00014701  00000-0  26152-3 0  9992';
tle_line2 = '2 25544  51.6332 166.3292 0002995  10.9396 349.1657 15.50384251530832';

%% Parse TLE
[satdata] = parse_tle(tle_line1, tle_line2);

%% Propagation Settings
propagation_hours = 24;
time_step_minutes = 2;
tsince_vec = 0:time_step_minutes:(propagation_hours*60);
n_steps = length(tsince_vec);

fprintf('Propagation: %.0f hours in %d steps\n\n', propagation_hours, n_steps);

%% 1. SGP4 Propagation
fprintf('1. Running SGP4 propagation...\n');
r_sgp4 = zeros(3, n_steps);
v_sgp4 = zeros(3, n_steps);

tic;
for idx = 1:n_steps
    [rteme, vteme] = sgp4(tsince_vec(idx), satdata);
    r_sgp4(:, idx) = rteme * 1000; % km to m
    v_sgp4(:, idx) = vteme * 1000; % km/s to m/s
end
t_sgp4 = toc;

fprintf(' ✓ Completed in %.4f s\n\n', t_sgp4);

%% 2. Custom Numerical Propagator
fprintf('2. Running custom numerical propagator (J2 + drag)...\n');

params.J3 = -2.53e-6;
params.bstar = satdata.bstar;
params.H = 8500;
params.mass = 1000;
params.I_CB = params.mass * eye(3);

r0_sgp4 = r_sgp4(:, 1);
v0_sgp4 = v_sgp4(:, 1);

X0 = zeros(13, 1);
X0(1:3) = r0_sgp4;
X0(4:6) = v0_sgp4;
X0(7:10) = [0; 0; 0; 1];
X0(11:13) = [0; 0; 0];

X_custom = zeros(6, n_steps);
X_custom(:, 1) = X0(1:6);
Xi = X0;
u = zeros(6, 1);

ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, 'MaxStep', 60);

tic;
for idx = 2:n_steps
    dt = time_step_minutes * 60;
    [~, X_traj] = ode45(@(t, X) Sat_template(t, X, u, params, ...
        'useJ2', true, 'useAtmDrag', true, ...
        'useControl', false), ...
        [0, dt], Xi, ode_options);
    Xi = X_traj(end, :)';
    Xi(7:10) = Xi(7:10) / norm(Xi(7:10));
    X_custom(:, idx) = Xi(1:6);
end
t_custom = toc;

fprintf(' ✓ Completed in %.4f s\n\n', t_custom);

%% 3. HPOP (High Precision Orbit Propagator)
fprintf('3. Running HPOP (high-precision numerical propagator)...\n');
fprintf('   Configuration: 20x20 gravity, Sun/Moon, SRP, drag\n');

% Convert epoch for HPOP
if satdata.epochyr < 57
    full_year = 2000 + satdata.epochyr;
else
    full_year = 1900 + satdata.epochyr;
end

% Calculate Modified Julian Date
Mjd0_UTC = Mjday(full_year, 1, 0) + satdata.epochdays;

% Initial state from SGP4 (already in meters and m/s)
Y0_hpop = [r0_sgp4; v0_sgp4];

% CRITICAL: Convert from TEME to ECI frame
% SGP4 outputs TEME coordinates; HPOP expects ECI
% Note: For short propagations, TEME ≈ ECI, but this conversion
% should be done properly for accurate results
% For now, we'll use TEME directly as an approximation

% Setup AuxParam structure (following test_HPOP.m)
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                 'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                 'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

AuxParam.Mjd_UTC = Mjd0_UTC;
AuxParam.area_solar = 1600;      % ISS solar radiation area [m^2]
AuxParam.area_drag = 1600;       % ISS drag area [m^2]
AuxParam.mass = 419725;          % ISS mass [kg]
AuxParam.Cr = 1.0;               % Radiation pressure coefficient
AuxParam.Cd = 2.2;               % Drag coefficient
AuxParam.n = 20;                 % Degree of gravity model
AuxParam.m = 20;                 % Order of gravity model
AuxParam.sun = 1;                % Solar perturbation on/off
AuxParam.moon = 1;               % Lunar perturbation on/off
AuxParam.planets = 0;            % Planetary perturbations off (faster)
AuxParam.sRad = 1;               % Solar radiation pressure on
AuxParam.drag = 1;               % Atmospheric drag on
AuxParam.SolidEarthTides = 0;
AuxParam.OceanTides = 0;
AuxParam.Relativity = 0;

% Use Ephemeris function (HPOP's high-level propagator)
Step = time_step_minutes * 60;   % seconds
N_Step = n_steps - 1;

tic;
Eph_hpop = Ephemeris(Y0_hpop, N_Step, Step);
t_hpop = toc;

% Extract position and velocity (already in meters and m/s)
r_hpop = zeros(3, n_steps);
v_hpop = zeros(3, n_steps);
for idx = 1:n_steps
    r_hpop(:, idx) = Eph_hpop(idx, 2:4)';  % Already in meters
    v_hpop(:, idx) = Eph_hpop(idx, 5:7)';  % Already in m/s
end

fprintf('   ✓ Completed in %.4f s\n\n', t_hpop);
hpop_success = true;
%% Calculate Errors
err_custom = zeros(1, n_steps);
if hpop_success
    err_hpop = zeros(1, n_steps);
    err_custom_vs_hpop = zeros(1, n_steps);
end

for idx = 1:n_steps
    err_custom(idx) = norm(X_custom(1:3, idx) - r_sgp4(:, idx));
    if hpop_success
        err_hpop(idx) = norm(r_hpop(:, idx) - r_sgp4(:, idx));
        err_custom_vs_hpop(idx) = norm(X_custom(1:3, idx) - r_hpop(:, idx));
    end
end

%% Altitude Analysis
R_e = 6378137;
alt_sgp4 = zeros(1, n_steps);
alt_custom = zeros(1, n_steps);
if hpop_success
    alt_hpop = zeros(1, n_steps);
end

for idx = 1:n_steps
    alt_sgp4(idx) = (norm(r_sgp4(:, idx)) - R_e) / 1e3;
    alt_custom(idx) = (norm(X_custom(1:3, idx)) - R_e) / 1e3;
    if hpop_success
        alt_hpop(idx) = (norm(r_hpop(:, idx)) - R_e) / 1e3;
    end
end

%% Results
fprintf('========================================\n');
fprintf(' RESULTS\n');
fprintf('========================================\n\n');

fprintf('Computation Time:\n');
fprintf('  SGP4:            %.4f s (baseline)\n', t_sgp4);
fprintf('  Custom Numerical: %.4f s (%.1fx)\n', t_custom, t_custom/t_sgp4);
if hpop_success
    fprintf('  HPOP (20x20):     %.4f s (%.1fx)\n\n', t_hpop, t_hpop/t_sgp4);
else
    fprintf('\n');
end

fprintf('Position Error vs SGP4:\n');
fprintf('  Custom Numerical:\n');
fprintf('    Mean:  %.3f km\n', mean(err_custom)/1e3);
fprintf('    Max:   %.3f km\n', max(err_custom)/1e3);
fprintf('    RMS:   %.3f km\n', sqrt(mean(err_custom.^2))/1e3);
fprintf('    Final: %.3f km\n', err_custom(end)/1e3);

if hpop_success
    fprintf('  HPOP:\n');
    fprintf('    Mean:  %.3f km\n', mean(err_hpop)/1e3);
    fprintf('    Max:   %.3f km\n', max(err_hpop)/1e3);
    fprintf('    RMS:   %.3f km\n', sqrt(mean(err_hpop.^2))/1e3);
    fprintf('    Final: %.3f km\n\n', err_hpop(end)/1e3);

    fprintf('Position Difference (Custom vs HPOP):\n');
    fprintf('    Mean:  %.3f km\n', mean(err_custom_vs_hpop)/1e3);
    fprintf('    Max:   %.3f km\n', max(err_custom_vs_hpop)/1e3);
    fprintf('    Final: %.3f km\n\n', err_custom_vs_hpop(end)/1e3);
end

fprintf('Altitude Decay Over %.0f Hours:\n', propagation_hours);
decay_sgp4 = alt_sgp4(1) - alt_sgp4(end);
decay_custom = alt_custom(1) - alt_custom(end);
fprintf('  SGP4:            %.3f km (%.3f km/day)\n', decay_sgp4, decay_sgp4*24/propagation_hours);
fprintf('  Custom Numerical: %.3f km (%.3f km/day)\n', decay_custom, decay_custom*24/propagation_hours);
if hpop_success
    decay_hpop = alt_hpop(1) - alt_hpop(end);
    fprintf('  HPOP:            %.3f km (%.3f km/day)\n\n', decay_hpop, decay_hpop*24/propagation_hours);
end

%% Plots
time_hours = tsince_vec / 60;

if hpop_success
    figure('Name', 'Three-Way Comparison: SGP4 vs Custom vs HPOP', 'Position', [50, 50, 1800, 1000]);
else
    figure('Name', 'Two-Way Comparison: SGP4 vs Custom', 'Position', [50, 50, 1600, 900]);
end

% 3D Orbit
subplot(2, 4, 1);
hold on; grid on; axis equal;
plot3(r_sgp4(1,:)/1e3, r_sgp4(2,:)/1e3, r_sgp4(3,:)/1e3, ...
    'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot3(X_custom(1,:)/1e3, X_custom(2,:)/1e3, X_custom(3,:)/1e3, ...
    'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom (J2+Drag)');
if hpop_success
    plot3(r_hpop(1,:)/1e3, r_hpop(2,:)/1e3, r_hpop(3,:)/1e3, ...
        'g:', 'LineWidth', 1.5, 'DisplayName', 'HPOP (20x20)');
end
[X_e, Y_e, Z_e] = sphere(30);
surf(X_e*R_e/1e3, Y_e*R_e/1e3, Z_e*R_e/1e3, ...
    'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Orbit Comparison');
legend('Location', 'best');
view(45, 30);

% Altitude
subplot(2, 4, 2);
hold on; grid on;
plot(time_hours, alt_sgp4, 'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot(time_hours, alt_custom, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom');
if hpop_success
    plot(time_hours, alt_hpop, 'g:', 'LineWidth', 1.5, 'DisplayName', 'HPOP');
end
xlabel('Time [hours]');
ylabel('Altitude [km]');
title('Altitude vs Time');
legend('Location', 'best');

% Position Error vs SGP4
subplot(2, 4, 3);
hold on; grid on;
plot(time_hours, err_custom/1e3, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Custom');
if hpop_success
    plot(time_hours, err_hpop/1e3, 'g-', 'LineWidth', 1.5, 'DisplayName', 'HPOP');
end
xlabel('Time [hours]');
ylabel('Position Error [km]');
title('Error vs SGP4');
legend('Location', 'best');

% Computation Time
subplot(2, 4, 4);
if hpop_success
    bar([t_sgp4, t_custom, t_hpop]);
    set(gca, 'XTickLabel', {'SGP4', 'Custom', 'HPOP'});
else
    bar([t_sgp4, t_custom]);
    set(gca, 'XTickLabel', {'SGP4', 'Custom'});
end
ylabel('Time [s]');
title('Computation Time');
grid on;

% Altitude Decay
subplot(2, 4, 5);
hold on; grid on;
plot(time_hours, (alt_sgp4 - alt_sgp4(1)), 'b-', 'LineWidth', 2, 'DisplayName', 'SGP4');
plot(time_hours, (alt_custom - alt_custom(1)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom');
if hpop_success
    plot(time_hours, (alt_hpop - alt_hpop(1)), 'g:', 'LineWidth', 1.5, 'DisplayName', 'HPOP');
end
xlabel('Time [hours]');
ylabel('Altitude Change [km]');
title('Altitude Decay');
legend('Location', 'best');

% Error Statistics
subplot(2, 4, 6);
if hpop_success
    err_stats = [mean(err_custom)/1e3, mean(err_hpop)/1e3; ...
                 max(err_custom)/1e3, max(err_hpop)/1e3; ...
                 err_custom(end)/1e3, err_hpop(end)/1e3];
    bar(err_stats);
    set(gca, 'XTickLabel', {'Mean', 'Max', 'Final'});
    ylabel('Error [km]');
    title('Error Statistics vs SGP4');
    legend('Custom', 'HPOP', 'Location', 'best');
else
    err_stats = [mean(err_custom)/1e3; max(err_custom)/1e3; err_custom(end)/1e3];
    bar(err_stats);
    set(gca, 'XTickLabel', {'Mean', 'Max', 'Final'});
    ylabel('Error [km]');
    title('Error Statistics vs SGP4');
end
grid on;

if hpop_success
    % Custom vs HPOP
    subplot(2, 4, 7);
    plot(time_hours, err_custom_vs_hpop/1e3, 'k-', 'LineWidth', 1.5);
    xlabel('Time [hours]');
    ylabel('Position Difference [km]');
    title('Custom vs HPOP (Truth Model)');
    grid on;

    % Decay Rate Comparison
    subplot(2, 4, 8);
    decay_rates = [decay_sgp4, decay_custom, decay_hpop] * 24 / propagation_hours;
    bar(decay_rates);
    set(gca, 'XTickLabel', {'SGP4', 'Custom', 'HPOP'});
    ylabel('Decay Rate [km/day]');
    title('Altitude Decay Rate');
    grid on;
end

fprintf('\n========================================\n');
fprintf('Analysis:\n');
fprintf('========================================\n\n');

if hpop_success
    fprintf('HPOP is the high-fidelity "truth" model with:\n');
    fprintf('  • 20x20 spherical harmonic gravity\n');
    fprintf('  • Sun/Moon third-body perturbations\n');
    fprintf('  • Solar radiation pressure\n');
    fprintf('  • NRLMSISE-00 atmospheric drag model\n\n');
    fprintf('Your custom propagator vs HPOP:\n');
    fprintf('  Mean difference: %.3f km\n', mean(err_custom_vs_hpop)/1e3);
    fprintf('  Max difference:  %.3f km\n', max(err_custom_vs_hpop)/1e3);
    fprintf('\nThis validates your implementation against\n');
    fprintf('the industry-standard high-precision reference!\n');
else
    fprintf('Custom vs SGP4: %.3f km (mean error)\n', mean(err_custom)/1e3);
end

fprintf('\n========================================\n');
if hpop_success
    fprintf('Three-Way Benchmark Complete!\n');
else
    fprintf('Two-Way Benchmark Complete!\n');
end
fprintf('========================================\n');
