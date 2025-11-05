%% MINIMUM-TIME OPTIMAL CONTROL WITH PONTRYAGIN'S MAXIMUM PRINCIPLE
% Fixed version - properly initializes X0 and params from control_test2 logic
% Objective: Reach target {RA, Dec, RA_dot, Dec_dot} in minimum time
% with constraint: |tau| <= 0.1 N per body axis

clear; close all; clc;

% ===================================================================
% STEP 1: Initialize Parameters and Scenario (from control_test2)
% ===================================================================

fprintf('\n=== MINIMUM-TIME OPTIMAL CONTROL ===\n');
fprintf('Using Pontryagin''s Maximum Principle\n\n');

% Parameters
params.mu = 3.986004418e14;
params.R_e = 6378137;
params.mass = 83.6;
params.radius = 0.58/2;
params.J2 = 1.08263e-3;
params.I_CB = 2/5 * params.mass * params.radius^2 * eye(3);
params.Kp_att = 10;
params.omega_earth = 7.2921159e-5;

% Orbital parameters
altitude = 200e3;
r_orbit = params.R_e + altitude;
v_circ = sqrt(params.mu / r_orbit);

% Initial orbital state 
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
orbital_period = 2 * pi * sqrt(r_orbit^3 / params.mu);

% Volcano position at t=0
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

fprintf('Scenario initialized:\n');
fprintf('  Orbital period: %.2f minutes\n', orbital_period/60);
fprintf('  Initial state X0 defined (13x1)\n\n');

% ===================================================================
% STEP 2: Find Target State (at visibility cutoff) - FAST METHOD
% ===================================================================

fprintf('Finding target state at visibility cutoff (zenith angle = 90°)...\n');

% FAST METHOD: Pre-compute trajectory once, evaluate visibility at all points
t_max_search = orbital_period * 2;  % Search up to 2 orbital periods
n_points = 500;
dt_search = t_max_search / (n_points - 1);
t_search_vec = 0:dt_search:t_max_search;

% Pre-propagate the satellite trajectory once (much faster!)
fprintf('Pre-propagating satellite trajectory...\n');
[t_traj, X_traj_full] = ode45(@(t,X) Sat_template(t, X, zeros(3,1), params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false), ...
    t_search_vec, X0, odeset('RelTol', 1e-8, 'AbsTol', 1e-10));

% Now evaluate visibility at all time points WITHOUT calling ode45 again
zenith_distance_array = zeros(length(t_search_vec), 1);

for idx = 1:length(t_search_vec)
    t = t_search_vec(idx);
    
    % Get pre-computed state (already from trajectory)
    [~, nearest_idx] = min(abs(t_traj - t));
    X_t = X_traj_full(nearest_idx, :)';
    
    % Current satellite position
    r_sat = X_t(1:3);
    
    % Volcano position at time t (rotates with Earth)
    Rotz_t = [cos(params.omega_earth * t), -sin(params.omega_earth * t), 0;
              sin(params.omega_earth * t), cos(params.omega_earth * t), 0;
              0, 0, 1];
    r_volcano_t = Rotz_t * r_volcano_0;
    
    % Calculate zenith distance
    zenith_dir = r_volcano_t / norm(r_volcano_t);
    sat_dir = (r_sat - r_volcano_t) / norm(r_sat - r_volcano_t);
    cos_zenith = dot(zenith_dir, sat_dir);
    zenith_distance_array(idx) = rad2deg(acos(max(-1, min(1, cos_zenith))));
end

% Find index where zenith distance crosses 90 degrees
idx_crossover = find(zenith_distance_array >= 90, 1, 'first');

if isempty(idx_crossover)
    error('Volcano does not go out of sight in search window. Increase t_max_search.');
end

t_target = t_search_vec(idx_crossover);

fprintf('Visibility cutoff at t_target = %.2f seconds (%.2f minutes)\n', t_target, t_target/60);

% Get target state from pre-computed trajectory
[~, target_idx] = min(abs(t_traj - t_target));
X_target_ref = X_traj_full(target_idx, :)';

% Extract target observation state
obs_target = state_to_observation(X_target_ref, params);
ra_target = deg2rad(obs_target.ra);
dec_target = deg2rad(obs_target.dec);

% Calculate target angular rates (numerical differentiation)
dt_rate = 1.0;  % 1 second
[~, X_traj_plus] = ode45(@(t,X) Sat_template(t, X, zeros(3,1), params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false), ...
    [t_target, t_target+dt_rate], X_target_ref, ...
    odeset('RelTol', 1e-8, 'AbsTol', 1e-10));

X_plus = X_traj_plus(end, :)';
obs_plus = state_to_observation(X_plus, params);
ra_dot_target = (deg2rad(obs_plus.ra) - ra_target) / dt_rate;
dec_dot_target = (deg2rad(obs_plus.dec) - dec_target) / dt_rate;

fprintf('\nTarget Observation State:\n');
fprintf(' RA = %.6f degrees\n', obs_target.ra);
fprintf(' Dec = %.6f degrees\n', obs_target.dec);
fprintf(' RA_dot = %.8f rad/s\n', ra_dot_target);
fprintf(' Dec_dot = %.8f rad/s\n\n', dec_dot_target);

% ===================================================================
% STEP 3: Set up Optimization Problem
% ===================================================================

fprintf('Setting up TPBVP optimization...\n');

% Discretize control into N intervals with piecewise constant torque
n_intervals = 25;
t_intervals = linspace(0, t_target, n_intervals+1);

% Optimization variables: [t_f; tau_interval_1; ...; tau_interval_N]
% where each tau_interval is 3D torque vector
n_vars = 1 + n_intervals * 3;

% Objective: minimize final time
objfun = @(x) x(1);

% Constraint: terminal state must match target
% This function will be called by fmincon
confun = @(x) constraint_terminal_state(x, X0, r_volcano_0, params, ...
    ra_target, dec_target, ra_dot_target, dec_dot_target, n_intervals);

% Initial guess
x0 = zeros(n_vars, 1);
x0(1) = t_target * 1.2;  % Start with 20% longer time
% tau = 0 (no control initially)

% Bounds
NEW_TAU_MAX = 0.001;  % N (change this value)
lb = [10; -NEW_TAU_MAX*ones(n_intervals*3, 1)];
ub = [t_target * 2; NEW_TAU_MAX*ones(n_intervals*3, 1)];

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 2000, ...
    'TolFun', 1e-5, ...
    'TolX', 1e-5, ...
    'StepTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-5, ...
    'UseParallel', false, ...
    'FiniteDifferenceType', 'central');

% ===================================================================
% STEP 4: Solve Optimization
% ===================================================================

fprintf('\nSolving with fmincon (SQP method)...\n');
fprintf('This may take several minutes...\n\n');

[x_opt, fval, exitflag, output, lambda, grad, hessian] = fmincon(...
    objfun, x0, [], [], [], [], lb, ub, confun, options);

t_f_optimal = x_opt(1);
tau_optimal = reshape(x_opt(2:end), n_intervals, 3);

fprintf('\n============================================\n');
fprintf('OPTIMAL SOLUTION FOUND\n');
fprintf('============================================\n');
fprintf('Minimum time: %.2f seconds (%.2f minutes)\n', t_f_optimal, t_f_optimal/60);
fprintf('Time improvement: %.1f%% faster than uncontrolled\n', ...
    (t_target - t_f_optimal)/t_target * 100);
fprintf('Exit flag: %d\n\n', exitflag);

% ===================================================================
% STEP 5: Simulate Optimal Trajectory
% ===================================================================

fprintf('Simulating optimal trajectory...\n');

n_sim_steps = 300;
t_sim = linspace(0, t_f_optimal, n_sim_steps);
X_opt_traj = zeros(n_sim_steps, 13);
obs_opt_traj = struct('ra', [], 'dec', []);
tau_applied_traj = zeros(n_sim_steps, 3);

X_current = X0;

for i = 1:n_sim_steps
    t_current = t_sim(i);
    X_opt_traj(i, :) = X_current';
    
    % Get observation
    obs = state_to_observation(X_current, params);
    obs_opt_traj(i).ra = obs.ra;
    obs_opt_traj(i).dec = obs.dec;
    
    % Interpolate control at current time
    interval_fraction = (t_current / t_f_optimal) * n_intervals;
    interval_idx = max(1, min(n_intervals, ceil(interval_fraction)));
    tau_current = tau_optimal(interval_idx, :)';
    tau_applied_traj(i, :) = tau_current';
    
    % Propagate to next time step
    if i < n_sim_steps
        dt_step = t_sim(i+1) - t_sim(i);
        [~, X_seg] = ode45(@(t,X) Sat_template(t, X, tau_current, params), ...
            [t_current, t_current+dt_step], X_current, ...
            odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
        X_current = X_seg(end, :)';
    end
end

% Final observation
obs_final = state_to_observation(X_current, params);

fprintf('\n============================================\n');
fprintf('FINAL STATE ACHIEVED\n');
fprintf('============================================\n');
fprintf('Final RA:  %.6f deg (target: %.6f deg) - Error: %.6f deg\n', ...
    obs_final.ra, obs_target.ra, abs(obs_final.ra - obs_target.ra));
fprintf('Final Dec: %.6f deg (target: %.6f deg) - Error: %.6f deg\n\n', ...
    obs_final.dec, obs_target.dec, abs(obs_final.dec - obs_target.dec));

% ===================================================================
% STEP 6: Plotting (Using Pre-computed Trajectory Data)
% ===================================================================

fprintf('Generating plots from pre-computed trajectory...\n\n');

% ===================================================================
% Plot 1: RA and Dec convergence
% ===================================================================

figure('Position', [100 100 1400 600]);

subplot(1, 2, 1);
plot(t_sim/60, [obs_opt_traj.ra], 'b-', 'LineWidth', 2.5); hold on;
yline(obs_target.ra, 'r--', 'LineWidth', 2, 'Label', 'Target');
ylabel('RA (degrees)'); xlabel('Time (minutes)');
title('Right Ascension vs Time (Optimal Control)');
grid on; legend('Trajectory', 'Target', 'Location', 'best');

subplot(1, 2, 2);
plot(t_sim/60, [obs_opt_traj.dec], 'b-', 'LineWidth', 2.5); hold on;
yline(obs_target.dec, 'r--', 'LineWidth', 2, 'Label', 'Target');
ylabel('Dec (degrees)'); xlabel('Time (minutes)');
title('Declination vs Time (Optimal Control)');
grid on; legend('Trajectory', 'Target', 'Location', 'best');

% ===================================================================
% Plot 2: Control torques
% ===================================================================

figure('Position', [100 100 1400 600]);

subplot(2, 1, 1);
plot(t_sim/60, tau_applied_traj(:,1), 'r-', 'LineWidth', 1.5); hold on;
plot(t_sim/60, tau_applied_traj(:,2), 'g-', 'LineWidth', 1.5);
plot(t_sim/60, tau_applied_traj(:,3), 'b-', 'LineWidth', 1.5);
yline(NEW_TAU_MAX, 'k--', 'LineWidth', 1);
yline(-NEW_TAU_MAX, 'k--', 'LineWidth', 1);
ylabel('Torque (N)'); title('Optimal Control Torques');
legend('\tau_x', '\tau_y', '\tau_z', 'Constraints', 'Location', 'best');
grid on;

subplot(2, 1, 2);
tau_mag = sqrt(sum(tau_applied_traj.^2, 2));
plot(t_sim/60, tau_mag, 'k-', 'LineWidth', 2);
ylabel('||τ|| (N)'); xlabel('Time (minutes)');
title('Control Magnitude vs Time');
grid on;

% ===================================================================
% Plot 3: Observation error convergence
% ===================================================================

figure('Position', [100 100 1000 600]);

ra_error = [obs_opt_traj.ra] - obs_target.ra;
dec_error = [obs_opt_traj.dec] - obs_target.dec;
total_error = sqrt(ra_error.^2 + dec_error.^2);

subplot(2, 1, 1);
plot(t_sim/60, ra_error, 'r-', 'LineWidth', 1.5); hold on;
plot(t_sim/60, dec_error, 'b-', 'LineWidth', 1.5);
yline(0, 'k--', 'LineWidth', 1);
ylabel('Error (degrees)'); title('Component Errors');
legend('RA error', 'Dec error', 'Location', 'best');
grid on;

subplot(2, 1, 2);
semilogy(t_sim/60, total_error, 'k-', 'LineWidth', 2);
ylabel('Total Error (degrees)'); xlabel('Time (minutes)');
title('Total Observation Error (log scale)');
grid on;

% % ===================================================================
% % Plot 4: 3D Animated Trajectory (Controlled)
% % ===================================================================
% 
% fprintf('Generating 3D trajectory animation...\n');
% 
% figure('Name', '3D Satellite and Volcano Trajectory (Optimal Control)', 'Position', [100 100 1200 800]);
% 
% % Earth sphere
% [xs, ys, zs] = sphere(50);
% xs = xs * params.R_e;
% ys = ys * params.R_e;
% zs = zs * params.R_e;
% 
% % Animation with frame skipping for smooth rendering
% n_frames = length(t_sim);
% frame_skip = max(1, floor(n_frames / 25));
% 
% for i = 1:frame_skip:n_frames
%     clf;
% 
%     % Draw Earth
%     surf(xs, ys, zs, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.3 0.5 0.8]);
%     hold on;
%     axis equal;
%     grid on;
%     xlabel('X (ECI) [m]'); 
%     ylabel('Y (ECI) [m]'); 
%     zlabel('Z (ECI) [m]');
% 
%     t = t_sim(i);
%     title(sprintf('Volcano Observation with Optimal Control - Time: %.2f / %.2f min', t/60, t_f_optimal/60));
% 
%     % Satellite position
%     r_sat = X_opt_traj(i, 1:3)';
%     plot3(r_sat(1), r_sat(2), r_sat(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% 
%     % Satellite trajectory (past)
%     if i > 1
%         plot3(X_opt_traj(1:i, 1), X_opt_traj(1:i, 2), X_opt_traj(1:i, 3), ...
%             'r-', 'LineWidth', 1.5);
%     end
% 
%     % Volcano position at time t
%     Rotz_t = [cos(params.omega_earth * t), -sin(params.omega_earth * t), 0;
%               sin(params.omega_earth * t), cos(params.omega_earth * t), 0;
%               0, 0, 1];
%     r_volcano_t = Rotz_t * r_volcano_0;
% 
%     plot3(r_volcano_t(1), r_volcano_t(2), r_volcano_t(3), ...
%         'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
% 
%     % Visibility check and line
%     visibility_check = dot(r_volcano_t, r_sat - r_volcano_t);
%     if visibility_check > 0
%         plot3([r_volcano_t(1), r_sat(1)], [r_volcano_t(2), r_sat(2)], ...
%             [r_volcano_t(3), r_sat(3)], 'g-', 'LineWidth', 2.5);
%         line_label = 'Line (Visible)';
%     else
%         plot3([r_volcano_t(1), r_sat(1)], [r_volcano_t(2), r_sat(2)], ...
%             [r_volcano_t(3), r_sat(3)], 'r--', 'LineWidth', 1.5);
%         line_label = 'Line (Not Visible)';
%     end
% 
%     legend('Earth', 'Satellite', 'Trajectory', 'Volcano', line_label, 'Location', 'best');
%     view(45, 30);
%     drawnow;
% end
% 
% % ===================================================================
% % Plot 5: Zenith Distance vs Time
% % ===================================================================
% 
% fprintf('Computing zenith distance...\n');
% 
% zenith_distance_opt = zeros(n_sim_steps, 1);
% 
% for i = 1:n_sim_steps
%     t = t_sim(i);
%     r_sat = X_opt_traj(i, 1:3)';
% 
%     % Volcano position at time t
%     Rotz_t = [cos(params.omega_earth * t), -sin(params.omega_earth * t), 0;
%               sin(params.omega_earth * t), cos(params.omega_earth * t), 0;
%               0, 0, 1];
%     r_volcano_t = Rotz_t * r_volcano_0;
% 
%     % Zenith distance calculation
%     zenith_at_volcano = r_volcano_t / norm(r_volcano_t);
%     volcano_to_sat = r_sat - r_volcano_t;
%     volcano_to_sat_unit = volcano_to_sat / norm(volcano_to_sat);
%     cos_zenith = dot(zenith_at_volcano, volcano_to_sat_unit);
%     zenith_distance_opt(i) = rad2deg(acos(max(-1, min(1, cos_zenith))));
% end
% 
% figure('Position', [100 100 1400 500]);
% 
% plot(t_sim/60, zenith_distance_opt, 'b-', 'LineWidth', 2.5); hold on;
% yline(90, 'r--', 'LineWidth', 2, 'Label', 'Visibility Cutoff (90°)');
% xlabel('Time (minutes)');
% ylabel('Zenith Distance (degrees)');
% title('Zenith Distance of Satellite from Volcano vs Time (Optimal Control)');
% ylim([0 180]);
% grid on;
% legend('Zenith Distance', 'Location', 'best');
% 
% % Mark final target time
% xline(t_target/60, 'k--', 'LineWidth', 1.5, 'Label', sprintf('Uncontrolled Target (%.2f min)', t_target/60));
% xline(t_f_optimal/60, 'g-', 'LineWidth', 2, 'Label', sprintf('Optimal Final Time (%.2f min)', t_f_optimal/60));

fprintf('\nAll plots generated successfully!\n\n');

fprintf('=== OPTIMIZATION COMPLETE ===\n');
fprintf('See figures for trajectory and control analysis\n');


% ===================================================================
% HELPER FUNCTION: Terminal State Constraint
% ===================================================================

function [c, ceq] = constraint_terminal_state(x, X0, r_volcano_0, params, ...
    ra_target, dec_target, ra_dot_target, dec_dot_target, n_intervals)

    t_f = x(1);
    tau_all = reshape(x(2:end), n_intervals, 3);
    
    % Propagate with piecewise constant control
    t_intervals = linspace(0, t_f, n_intervals+1);
    X_current = X0;
    
    for i = 1:n_intervals
        t_start = t_intervals(i);
        t_end = t_intervals(i+1);
        tau_i = tau_all(i, :)';
        
        % Propagate interval i
        [~, X_seg] = ode45(@(t,X) Sat_template(t, X, tau_i, params), ...
            [t_start, t_end], X_current, ...
            odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
        
        X_current = X_seg(end, :)';
    end
    
    % Get final observation
    obs_final = state_to_observation(X_current, params);
    ra_final = deg2rad(obs_final.ra);
    dec_final = deg2rad(obs_final.dec);
    
    % Compute final angular rates (numerical differentiation)
    dt_check = 0.5;
    [~, X_seg_plus] = ode45(@(t,X) Sat_template(t, X, zeros(3,1), params), ...
        [t_f, t_f+dt_check], X_current, ...
        odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    X_plus = X_seg_plus(end, :)';
    
    obs_plus = state_to_observation(X_plus, params);
    ra_dot_final = (deg2rad(obs_plus.ra) - ra_final) / dt_check;
    dec_dot_final = (deg2rad(obs_plus.dec) - dec_final) / dt_check;
    
    % Equality constraints (terminal state matching)
    ceq = [
        ra_final - ra_target;
        dec_final - dec_target;
        ra_dot_final - ra_dot_target;
        dec_dot_final - dec_dot_target
    ];
    
    % No inequality constraints
    c = [];
end