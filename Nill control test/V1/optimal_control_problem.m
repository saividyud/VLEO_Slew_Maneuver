% Minimum Time Trajectory with Bounded Acceleration
% Given: x0=[0;0;0], v0, x1, v1, amax
% Find: minimum time tmin

clear; clc; close all;

% Problem parameters
x0 = [0; 0; 0];           % Initial position (origin)
v0 = [1; 0.5; 0.2];       % Initial velocity (random)
x1 = [5; 3; 2];           % Final position
v1 = [-0.5; 0.3; 0.1];    % Final velocity (random)
amax = 2.0;               % Maximum acceleration magnitude

% Scaling factors for better numerical conditioning
L_scale = norm(x1 - x0);  % Length scale
if L_scale < 1e-6
    L_scale = 1.0;
end
V_scale = max(norm(v0), norm(v1));  % Velocity scale
if V_scale < 1e-6
    V_scale = 1.0;
end
T_scale = L_scale / V_scale;  % Time scale

% Scale the problem
x0_s = x0 / L_scale;
x1_s = x1 / L_scale;
v0_s = v0 / V_scale;
v1_s = v1 / V_scale;
amax_s = amax * T_scale / V_scale;

% Estimate initial guess for final time (in scaled units)
dx = x1_s - x0_s;
dist = norm(dx);
avg_vel = (norm(v0_s) + norm(v1_s)) / 2;
tf_guess = dist / max(avg_vel, 0.1) * 2;

% Improved initial conditions with better scaling
n_mesh = 50;
solinit = bvpinit(linspace(0, 1, n_mesh), ...
    @(tau) guess_init(tau, x0_s, v0_s, x1_s, v1_s, tf_guess));

% Improved bvp4c options for better convergence
options = bvpset('RelTol', 1e-3, 'AbsTol', 1e-5, 'NMax', 10000);

try
    sol = bvp4c(@(tau,y) ode_sys(tau, y, amax_s), ...
                @(ya,yb) bc_conditions(ya, yb, x0_s, v0_s, x1_s, v1_s), ...
                solinit, options);
    
    % Extract solution
    tau = linspace(0, 1, 500);
    y = deval(sol, tau);
    tf_s = y(13, 1);  % Final time (scaled)
    
    % Unscale the solution
    tf = tf_s * T_scale;
    t = tau * tf;
    x_traj = y(1:3, :) * L_scale;
    v_traj = y(4:6, :) * V_scale;
    lambda_x = y(7:9, :);
    lambda_v = y(10:12, :);
    
    % Compute optimal acceleration
    a_traj = zeros(3, length(t));
    for i = 1:length(t)
        lambda_norm = norm(lambda_x(:,i));
        if lambda_norm > 1e-10
            a_traj(:,i) = -amax * lambda_x(:,i) / lambda_norm;
        end
    end
    
    % Compute acceleration magnitude
    a_mag = vecnorm(a_traj);
    
    fprintf('Minimum time: %.4f seconds\n', tf);
    fprintf('Initial position: [%.2f, %.2f, %.2f]\n', x0);
    fprintf('Final position: [%.2f, %.2f, %.2f]\n', x_traj(:,end));
    fprintf('Position error: %.6f\n', norm(x_traj(:,end) - x1));
    fprintf('Velocity error: %.6f\n', norm(v_traj(:,end) - v1));
    fprintf('Acceleration magnitude (max): %.4f m/s^2\n', max(a_mag));
    fprintf('Acceleration magnitude (mean): %.4f m/s^2\n', mean(a_mag));
    
    % Plot acceleration magnitude and lambda_r
    plot_analysis(t, a_traj, a_mag, lambda_x, lambda_v, amax);
    
    % Animate trajectory with initial velocities passed
    animate_trajectory(t, x_traj, v_traj, a_traj, x0, x1, v0, amax);
    
catch ME
    fprintf('BVP solver failed: %s\n', ME.message);
    fprintf('Try adjusting initial guess or problem parameters.\n');
end

%% Helper Functions

function dydt = ode_sys(tau, y, amax_s)
    % State: y = [x(1:3); v(4:6); lambda_r(7:9); lambda_v(10:12); tf(13)]
    v = y(4:6);
    lambda_r = y(7:9);  % Changed name for clarity
    lambda_v = y(10:12);
    tf = y(13);
    
    % Optimal control (bang-bang)
    % Note: We use lambda_r (not lambda_v!) for control direction
    lambda_r_norm = norm(lambda_r);
    lambda_v_norm = norm(lambda_v);
    if lambda_r_norm > 1e-10
        a = amax_s * lambda_v / lambda_v_norm;
    else
        a = [0; 0; 0];
    end
    
    % State equations: dx/dtau = tf * dx/dt
    dxdt = v;
    dvdt = a;
    
    % Costate equations (CORRECTED!)
    dlambda_rdt = zeros(3,1);  % lambda_r is constant!
    dlambda_vdt = -lambda_r;   % lambda_v evolves linearly
    
    % Final time is constant
    dtfdt = 0;
    
    dydt = tf * [dxdt; dvdt; dlambda_rdt; dlambda_vdt; dtfdt];
end

function res = bc_conditions(ya, yb, x0_s, v0_s, x1_s, v1_s)
    % Initial conditions
    res_initial = [
        ya(1:3) - x0_s;
        ya(4:6) - v0_s
    ];
    
    % Final conditions
    res_final = [
        yb(1:3) - x1_s;
        yb(4:6) - v1_s
    ];
    
    % Transversality condition: H(tf) = 0
    lambda_r_f = yb(7:9);
    lambda_v_f = yb(10:12);
    v_f = yb(4:6);
    
    % At optimal control: a = amax * lambda_v / |lambda_v|
    % H = lambda_r^T v + lambda_v^T a + 1
    %   = lambda_r^T v + |lambda_v| * amax + 1
    H_final = dot(lambda_r_f, v_f) + norm(lambda_v_f) + 1;
    
    res = [res_initial; res_final; H_final];
end



function y0 = guess_init(tau, x0_s, v0_s, x1_s, v1_s, tf_guess)
    % Improved initial guess using polynomial interpolation
    
    % Position: cubic polynomial to match positions and velocities
    x_guess = x0_s + v0_s * tau * tf_guess + ...
              (3 * (x1_s - x0_s) / tf_guess - 2 * v0_s - v1_s) * tau^2 * tf_guess + ...
              (2 * (x0_s - x1_s) / tf_guess + v0_s + v1_s) * tau^3 * tf_guess;
    
    % Velocity: derivative of position
    v_guess = v0_s + ...
              2 * (3 * (x1_s - x0_s) / tf_guess - 2 * v0_s - v1_s) * tau + ...
              3 * (2 * (x0_s - x1_s) / tf_guess + v0_s + v1_s) * tau^2;
    
    % Costates: rough estimate based on needed acceleration
    a_guess = (v1_s - v0_s) / tf_guess;
    lambda_x_guess = -a_guess / max(norm(a_guess), 1e-6);
    lambda_v_guess = -(x1_s - x0_s) / (tf_guess^2);
    
    y0 = [x_guess; v_guess; lambda_x_guess; lambda_v_guess; tf_guess];
end

function plot_analysis(t, a_traj, a_mag, lambda_x, lambda_v, amax)
    figure('Position', [100 100 1400 800]);
    
    % Acceleration magnitude
    subplot(3,2,1);
    plot(t, a_mag, 'b-', 'LineWidth', 2); hold on;
    plot(t, amax*ones(size(t)), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Acceleration Magnitude [m/s^2]');
    title('Acceleration Magnitude vs Time');
    legend('||a(t)||', 'a_{max}', 'Location', 'best');
    
    % Acceleration components
    subplot(3,2,2);
    plot(t, a_traj(1,:), 'b-', 'LineWidth', 2); hold on;
    plot(t, a_traj(2,:), 'r-', 'LineWidth', 2);
    plot(t, a_traj(3,:), 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('Acceleration [m/s^2]');
    title('Acceleration Components vs Time');
    legend('a_x', 'a_y', 'a_z', 'Location', 'best');
    
    % Lambda_r (costate for position) magnitude
    subplot(3,2,3);
    lambda_r_mag = vecnorm(lambda_x);
    plot(t, lambda_r_mag, 'm-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('||\lambda_r|| (scaled)');
    title('Position Costate Magnitude vs Time');
    
    % Lambda_r components
    subplot(3,2,4);
    plot(t, lambda_x(1,:), 'b-', 'LineWidth', 2); hold on;
    plot(t, lambda_x(2,:), 'r-', 'LineWidth', 2);
    plot(t, lambda_x(3,:), 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('\lambda_r (scaled)');
    title('Position Costate Components vs Time');
    legend('\lambda_{r,x}', '\lambda_{r,y}', '\lambda_{r,z}', 'Location', 'best');
    
    % Lambda_v (costate for velocity) magnitude
    subplot(3,2,5);
    lambda_v_mag = vecnorm(lambda_v);
    plot(t, lambda_v_mag, 'c-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('||\lambda_v|| (scaled)');
    title('Velocity Costate Magnitude vs Time');
    
    % Lambda_v components
    subplot(3,2,6);
    plot(t, lambda_v(1,:), 'b-', 'LineWidth', 2); hold on;
    plot(t, lambda_v(2,:), 'r-', 'LineWidth', 2);
    plot(t, lambda_v(3,:), 'g-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('\lambda_v (scaled)');
    title('Velocity Costate Components vs Time');
    legend('\lambda_{v,x}', '\lambda_{v,y}', '\lambda_{v,z}', 'Location', 'best');
    
    sgtitle('Optimal Control Analysis: Acceleration and Costates');
end

function animate_trajectory(t, x_traj, v_traj, a_traj, x0, x1, v0_init, amax)
    figure('Position', [100 100 1400 600]);
    
    % 3D trajectory plot
    subplot(1,2,1);
    plot3(x_traj(1,:), x_traj(2,:), x_traj(3,:), 'b-', 'LineWidth', 2);
    hold on;
    plot3(x0(1), x0(2), x0(3), 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
    plot3(x1(1), x1(2), x1(3), 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
    
    % Animated point
    h_point = plot3(x0(1), x0(2), x0(3), 'ko', 'MarkerSize', 10, ...
                    'MarkerFaceColor', 'k');
    
    % Create quiver objects for velocity and acceleration
    v_scale = 0.5;
    a_scale = 0.3;
    h_vel = quiver3(x0(1), x0(2), x0(3), ...
                    v_scale*v0_init(1), v_scale*v0_init(2), v_scale*v0_init(3), ...
                    'r', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    h_acc = quiver3(x0(1), x0(2), x0(3), ...
                    a_scale*a_traj(1,1), a_scale*a_traj(2,1), a_scale*a_traj(3,1), ...
                    'm', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
    
    grid on; axis equal;
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    title('3D Trajectory Animation');
    legend('Trajectory', 'Start', 'End', 'Current Position', ...
           'Velocity', 'Acceleration', 'Location', 'best');
    view(45, 30);
    
    % State vs time plots
    subplot(1,2,2);
    plot(t, x_traj(1,:), 'b-', 'LineWidth', 2); hold on;
    plot(t, x_traj(2,:), 'r-', 'LineWidth', 2);
    plot(t, x_traj(3,:), 'g-', 'LineWidth', 2);
    plot(t, v_traj(1,:), 'b--', 'LineWidth', 1.5);
    plot(t, v_traj(2,:), 'r--', 'LineWidth', 1.5);
    plot(t, v_traj(3,:), 'g--', 'LineWidth', 1.5);
    ylims = ylim;
    h_time = plot([0 0], ylims, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('Time [s]');
    ylabel('State Variables');
    legend('x', 'y', 'z', 'v_x', 'v_y', 'v_z', 'Current Time', ...
           'Location', 'best');
    title('Position and Velocity vs Time');
    
    % Animation loop
    dt_anim = 0.05;  % Animation time step
    skip = max(1, floor(length(t) / 100));  % Skip points for smoother animation
    
    for i = 1:skip:length(t)
        % Update 3D plot
        subplot(1,2,1);
        set(h_point, 'XData', x_traj(1,i), 'YData', x_traj(2,i), ...
                    'ZData', x_traj(3,i));
        
        % Update velocity quiver (using proper Quiver object properties)
        set(h_vel, 'XData', x_traj(1,i), ...
                   'YData', x_traj(2,i), ...
                   'ZData', x_traj(3,i), ...
                   'UData', v_scale*v_traj(1,i), ...
                   'VData', v_scale*v_traj(2,i), ...
                   'WData', v_scale*v_traj(3,i));
        
        % Update acceleration quiver
        set(h_acc, 'XData', x_traj(1,i), ...
                   'YData', x_traj(2,i), ...
                   'ZData', x_traj(3,i), ...
                   'UData', a_scale*a_traj(1,i), ...
                   'VData', a_scale*a_traj(2,i), ...
                   'WData', a_scale*a_traj(3,i));
        
        % Update time indicator
        subplot(1,2,2);
        set(h_time, 'XData', [t(i) t(i)], 'YData', ylims);
        
        drawnow;
        pause(dt_anim);
    end
    
    fprintf('\nAnimation complete!\n');
end
