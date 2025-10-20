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
    
    % Compute optimal acceleration
    a_traj = zeros(3, length(t));
    for i = 1:length(t)
        lambda_norm = norm(lambda_x(:,i));
        if lambda_norm > 1e-10
            a_traj(:,i) = -amax * lambda_x(:,i) / lambda_norm;
        end
    end
    
    fprintf('Minimum time: %.4f seconds\n', tf);
    fprintf('Initial position: [%.2f, %.2f, %.2f]\n', x0);
    fprintf('Final position: [%.2f, %.2f, %.2f]\n', x_traj(:,end));
    fprintf('Position error: %.6f\n', norm(x_traj(:,end) - x1));
    fprintf('Velocity error: %.6f\n', norm(v_traj(:,end) - v1));
    fprintf('Acceleration magnitude (max): %.4f\n', max(vecnorm(a_traj)));
    
    % Animate trajectory with initial velocities passed
    animate_trajectory(t, x_traj, v_traj, a_traj, x0, x1, v0, amax);
    
catch ME
    fprintf('BVP solver failed: %s\n', ME.message);
    fprintf('Try adjusting initial guess or problem parameters.\n');
end

%% Helper Functions

function dydt = ode_sys(tau, y, amax_s)
    % State: y = [x(1:3); v(4:6); lambda_x(7:9); lambda_v(10:12); tf(13)]
    v = y(4:6);
    lambda_x = y(7:9);
    lambda_v = y(10:12);
    tf = y(13);
    
    % Optimal control (bang-bang)
    lambda_norm = norm(lambda_x);
    if lambda_norm > 1e-10
        a = -amax_s * lambda_x / lambda_norm;
    else
        a = [0; 0; 0];
    end
    
    % State equations: dx/dtau = tf * dx/dt
    dxdt = v;
    dvdt = a;
    
    % Costate equations: dlambda/dtau = -tf * dH/dx
    dlambda_xdt = -lambda_v;  % -dH/dx
    dlambda_vdt = -lambda_x;  % -dH/dv
    
    % Final time is constant
    dtfdt = 0;
    
    dydt = tf * [dxdt; dvdt; dlambda_xdt; dlambda_vdt; dtfdt];
end

function res = bc_conditions(ya, yb, x0_s, v0_s, x1_s, v1_s)
    % Boundary conditions at tau=0 and tau=1
    
    % Initial conditions
    res_initial = [
        ya(1:3) - x0_s;  % x(0) = x0
        ya(4:6) - v0_s   % v(0) = v0
    ];
    
    % Final conditions
    res_final = [
        yb(1:3) - x1_s;  % x(tf) = x1
        yb(4:6) - v1_s   % v(tf) = v1
    ];
    
    % Transversality condition for free final time
    lambda_x_f = yb(7:9);
    lambda_v_f = yb(10:12);
    v_f = yb(4:6);
    
    % For bang-bang: H = 1 + lambda_v^T v - |lambda_x|
    H_final = 1 + dot(lambda_v_f, v_f) - norm(lambda_x_f);
    
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
