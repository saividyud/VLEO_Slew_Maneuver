projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

clc;
clear;
close all;

c = vleo.util.constants();
[t, X] = simulate_attitude_history(c);

torques = zeros(numel(t), 3);
for idx = 1:numel(t)
    torques(idx, :) = vleo.control.attitude_pd_controller(t(idx), X(idx, :)')';
end

results = struct( ...
    't', t, ...
    'rs', X(:, 1:3), ...
    'vs', X(:, 4:6), ...
    'betas', X(:, 7:10), ...
    'omegas', X(:, 11:13), ...
    'torques', torques);
vleo.viz.animate_results(struct('results', results));

function [t, X] = simulate_attitude_history(c)
    duration = 30 * 60;
    timeStep = 1;
    ts = 0:timeStep:duration;

    [rEci, vEci] = vleo.analysis.keplerian_to_eci_safe( ...
        c.R_earth + 250e3, 0, 0, 0, 0, 0, ...
        'GravitationalParameter', c.mu_earth, 'Action', 'None');

    b1 = rvec(vEci / norm(vEci));
    b3 = rvec(-rEci / norm(rEci));
    b2 = cross(b3, b1);
    b2 = b2 / norm(b2);
    rBodyFromInertial = [b1'; b2'; b3'];
    q0 = dcm2quat(rBodyFromInertial);
    if q0(1) < 0
        q0 = -q0;
    end

    X0 = [rvec(rEci); rvec(vEci); q0(:); 0; 0; 0];
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [t, X] = ode45(@vleo.dynamics.sat_dynamics_nonlinear, ts, X0, opts);
end

function v = rvec(v)
    v = reshape(v, [], 1);
end
