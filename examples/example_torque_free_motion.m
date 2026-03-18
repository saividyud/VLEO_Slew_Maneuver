%% example_torque_free_motion.m
%  Propagate a rigid body with no external torque and plot Euler angles.
%  Demonstrates torque-free precession / nutation for a 6U CubeSat with
%  asymmetric inertia and an initial angular rate off the principal axes.
%  Also shows sensitivity to +/-1% perturbation in initial body rates,
%  and compares the nonlinear error against a linearized error model.
%
%  Linearized error dynamics (torque-free, about nominal trajectory):
%    state: dx = [dtheta(3x1); domega(3x1)]
%    dtheta_dot = domega
%    domega_dot = J^{-1} * ( skew(J*w_nom)*domega - skew(w_nom)*J*domega )
%  where w_nom(t) is the nominal body rate at time t, obtained by
%  interpolating the nominal nonlinear trajectory.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

clc;
close all;

%% Parameters
c = vleo.util.constants();
vehicle = vleo.dynamics.default_vehicle_config();
I_body = vehicle.inertiaBody;
Jinv = inv(I_body);
mu = c.mu_earth;

fprintf('=== Torque-Free Rigid Body Motion ===\n');
fprintf('Inertia [kg m^2]: diag([%.4f, %.4f, %.4f])\n', ...
    I_body(1,1), I_body(2,2), I_body(3,3));

%% Initial orbit (250 km circular, equatorial)
altitude_km = 250;
[rEci, vEci] = vleo.analysis.keplerian_to_eci_safe( ...
    c.R_earth + altitude_km * 1e3, 0, 0, 0, 0, 0);

%% Initial attitude — identity quaternion (body aligned with ECI)
q0 = [1; 0; 0; 0];  % scalar-first: [q0; q1; q2; q3]

%% Initial angular rate — off-axis spin to trigger precession
omega0_deg_per_sec = [2; 5; 1];  % [deg/s]
omega0 = deg2rad(omega0_deg_per_sec);

fprintf('Initial body rates [deg/s]: [%.1f, %.1f, %.1f]\n', ...
    omega0_deg_per_sec(1), omega0_deg_per_sec(2), omega0_deg_per_sec(3));

%% Propagation settings
tEnd_sec = 120;
dt = 0.1;
tSpan = 0:dt:tEnd_sec;
odeOpts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
odeFun = @(t, X) vleo.dynamics.sat_dynamics_base(t, X, mu, [], [0;0;0], I_body, false);

%% Propagate three nonlinear cases: nominal, +1%, -1%
perturbFrac = 0.01;
cases = struct( ...
    'name',  {'Nominal', '+1%', '-1%'}, ...
    'omega0', {omega0, omega0 * (1 + perturbFrac), omega0 * (1 - perturbFrac)}, ...
    'color', {'b', 'r', [0 0.6 0]}, ...
    'style', {'-', '--', ':'});

nCases = numel(cases);
for ic = 1:nCases
    X0 = [rEci(:); vEci(:); q0; cases(ic).omega0];
    fprintf('Propagating %s case (nonlinear) ...\n', cases(ic).name);
    [tOut, XHist] = ode45(odeFun, tSpan, X0, odeOpts);

    qHist = vleo.util.normalize_quaternion_rows(XHist(:, 7:10));
    qHist = vleo.util.align_quaternion_signs(qHist);

    cases(ic).tOut = tOut;
    cases(ic).qHist = qHist;
    cases(ic).eulerDeg = vleo.util.quat_history_to_euler_deg(qHist);
    cases(ic).omegaRad = XHist(:, 11:13);
    cases(ic).omegaDeg = rad2deg(XHist(:, 11:13));
end

%% Compute nonlinear error (true difference)
%  Attitude error via error quaternion: q_err = inv(q_nom) * q_pert
%  For small angles, the vector part ≈ 0.5 * dtheta
%  Rate error is simply domega = omega_pert - omega_nom
for ic = 2:nCases
    nPts = numel(cases(1).tOut);
    cases(ic).dthetaNL_deg = zeros(nPts, 3);
    cases(ic).domegaNL_deg = zeros(nPts, 3);
    for k = 1:nPts
        qNom = cases(1).qHist(k, :)';   % scalar-first [q0;q1;q2;q3]
        qPert = cases(ic).qHist(k, :)';
        qErr = quat_multiply(quat_conj(qNom), qPert);
        % Ensure scalar part positive for consistent small-angle extraction
        if qErr(1) < 0
            qErr = -qErr;
        end
        % dtheta ≈ 2 * q_vec for small angles
        cases(ic).dthetaNL_deg(k, :) = rad2deg(2 * qErr(2:4)');
        cases(ic).domegaNL_deg(k, :) = cases(ic).omegaDeg(k, :) - cases(1).omegaDeg(k, :);
    end
end

%% Propagate linearized error dynamics
%  dx = [dtheta(3); domega(3)]
%  dtheta_dot = domega
%  domega_dot = J^{-1} * ( skew(J*w_nom)*domega - skew(w_nom)*J*domega )
%
%  This is the linearization of Euler's equation:
%    J*omega_dot = -omega x (J*omega)
%  about omega_nom(t), keeping only first-order terms in domega.

% Build interpolant for nominal body rates
tNom = cases(1).tOut;
omegaNomRad = cases(1).omegaRad;  % [N x 3]

odeOptsLin = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

for ic = 2:nCases
    domega0 = cases(ic).omega0 - omega0;  % initial rate perturbation
    dtheta0 = [0; 0; 0];                  % no initial attitude error
    dx0 = [dtheta0; domega0];

    fprintf('Propagating %s case (linearized) ...\n', cases(ic).name);
    linOdeFun = @(t, dx) linearized_error_dynamics( ...
        t, dx, I_body, Jinv, tNom, omegaNomRad);
    [~, dxHist] = ode45(linOdeFun, tSpan, dx0, odeOptsLin);

    cases(ic).dthetaLin_deg = rad2deg(dxHist(:, 1:3));
    cases(ic).domegaLin_deg = rad2deg(dxHist(:, 4:6));
end

%% Verify conserved quantities (nominal case)
H_norm = zeros(size(cases(1).tOut));
T_rot = zeros(size(cases(1).tOut));
for k = 1:numel(cases(1).tOut)
    w = cases(1).omegaRad(k, :)';
    H_norm(k) = norm(I_body * w);
    T_rot(k) = 0.5 * w' * I_body * w;
end
fprintf('Nominal angular momentum drift: %.2e (relative)\n', ...
    max(abs(H_norm - H_norm(1))) / H_norm(1));
fprintf('Nominal rotational energy drift: %.2e (relative)\n', ...
    max(abs(T_rot - T_rot(1))) / T_rot(1));

%% Plot 1: Euler angles — all three nonlinear cases overlaid
axLabels = {'Roll [deg]', 'Pitch [deg]', 'Yaw [deg]'};
figure('Name', 'Torque-Free Euler Angles', 'Position', [100 100 800 600]);
for ax = 1:3
    subplot(3, 1, ax);
    hold on;
    for ic = 1:nCases
        plot(cases(ic).tOut, cases(ic).eulerDeg(:, ax), ...
            cases(ic).style, 'Color', cases(ic).color, 'LineWidth', 1.3);
    end
    ylabel(axLabels{ax});
    grid on;
    if ax == 1
        title('Torque-Free Euler Angles (ZYX) — Nominal vs \pm1% Rate Perturbation');
        legend({cases.name}, 'Location', 'best');
    end
end
xlabel('Time [s]');

%% Plot 2: Body rates — all three nonlinear cases overlaid
rateLabels = {'\omega_x [deg/s]', '\omega_y [deg/s]', '\omega_z [deg/s]'};
figure('Name', 'Torque-Free Body Rates', 'Position', [100 100 800 600]);
for ax = 1:3
    subplot(3, 1, ax);
    hold on;
    for ic = 1:nCases
        plot(cases(ic).tOut, cases(ic).omegaDeg(:, ax), ...
            cases(ic).style, 'Color', cases(ic).color, 'LineWidth', 1.3);
    end
    ylabel(rateLabels{ax});
    grid on;
    if ax == 1
        title('Torque-Free Body Rates — Nominal vs \pm1% Rate Perturbation');
        legend({cases.name}, 'Location', 'best');
    end
end
xlabel('Time [s]');

%% Plot 3: Attitude error — nonlinear vs linearized
errLabels = {'\delta\theta_x [deg]', '\delta\theta_y [deg]', '\delta\theta_z [deg]'};
for ic = 2:nCases
    figure('Name', ['Attitude Error — ' cases(ic).name], ...
        'Position', [100 100 800 600]);
    for ax = 1:3
        subplot(3, 1, ax);
        hold on;
        plot(cases(1).tOut, cases(ic).dthetaNL_deg(:, ax), ...
            '-', 'Color', cases(ic).color, 'LineWidth', 1.5);
        plot(cases(1).tOut, cases(ic).dthetaLin_deg(:, ax), ...
            '--k', 'LineWidth', 1.3);
        ylabel(errLabels{ax});
        grid on;
        if ax == 1
            title(['Attitude Error (' cases(ic).name ...
                ') — Nonlinear (solid) vs Linearized (dashed)']);
            legend('Nonlinear', 'Linearized', 'Location', 'best');
        end
    end
    xlabel('Time [s]');
end

%% Plot 4: Rate error — nonlinear vs linearized
dRateLabels = {'\delta\omega_x [deg/s]', '\delta\omega_y [deg/s]', '\delta\omega_z [deg/s]'};
for ic = 2:nCases
    figure('Name', ['Rate Error — ' cases(ic).name], ...
        'Position', [100 100 800 600]);
    for ax = 1:3
        subplot(3, 1, ax);
        hold on;
        plot(cases(1).tOut, cases(ic).domegaNL_deg(:, ax), ...
            '-', 'Color', cases(ic).color, 'LineWidth', 1.5);
        plot(cases(1).tOut, cases(ic).domegaLin_deg(:, ax), ...
            '--k', 'LineWidth', 1.3);
        ylabel(dRateLabels{ax});
        grid on;
        if ax == 1
            title(['Rate Error (' cases(ic).name ...
                ') — Nonlinear (solid) vs Linearized (dashed)']);
            legend('Nonlinear', 'Linearized', 'Location', 'best');
        end
    end
    xlabel('Time [s]');
end

%% Plot 5: Linearization quality — difference between nonlinear and linearized error
for ic = 2:nCases
    figure('Name', ['Linearization Residual — ' cases(ic).name], ...
        'Position', [100 100 800 600]);
    for ax = 1:3
        subplot(3, 1, ax);
        hold on;
        plot(cases(1).tOut, ...
            cases(ic).dthetaNL_deg(:, ax) - cases(ic).dthetaLin_deg(:, ax), ...
            '-', 'Color', cases(ic).color, 'LineWidth', 1.3);
        ylabel(['\Delta' errLabels{ax}]);
        grid on;
        if ax == 1
            title(['Linearization Residual (' cases(ic).name ...
                ') — Nonlinear minus Linearized']);
        end
    end
    xlabel('Time [s]');
end

%% Plot 6: Conserved quantities (nominal)
figure('Name', 'Conserved Quantities', 'Position', [100 100 800 400]);

subplot(2, 1, 1);
plot(cases(1).tOut, H_norm, 'b-', 'LineWidth', 1.3);
ylabel('|H| [N m s]');
title('Angular Momentum and Rotational Energy (nominal, should be constant)');
grid on;

subplot(2, 1, 2);
plot(cases(1).tOut, T_rot * 1e3, 'r-', 'LineWidth', 1.3);
ylabel('T_{rot} [mJ]');
xlabel('Time [s]');
grid on;

fprintf('Done.\n');

%% ---- Local helper functions ----

function dxd = linearized_error_dynamics(t, dx, J, Jinv, tNom, omegaNomRad)
    % Linearized torque-free attitude error dynamics.
    %   dx = [dtheta(3x1); domega(3x1)]
    %
    % Kinematics (Markley & Crassidis):
    %   dtheta_dot = domega - skew(w_nom) * dtheta
    % The cross term accounts for the fact that dtheta is expressed in the
    % rotating nominal body frame.
    %
    % Dynamics (linearized Euler equation, torque-free):
    %   domega_dot = J^{-1} * [ skew(J*w_nom)*domega - skew(w_nom)*J*domega ]

    dtheta = dx(1:3);
    domega = dx(4:6);

    % Interpolate nominal angular rate at current time
    wNom = interp1(tNom, omegaNomRad, t, 'pchip')';

    Jw = J * wNom;

    dtheta_dot = domega - skew(wNom) * dtheta;
    domega_dot = Jinv * (skew(Jw) * domega - skew(wNom) * J * domega);

    dxd = [dtheta_dot; domega_dot];
end

function S = skew(v)
    S = [  0    -v(3)  v(2);
          v(3)   0    -v(1);
         -v(2)  v(1)   0  ];
end

function qc = quat_conj(q)
    % Conjugate of scalar-first quaternion [q0; q1; q2; q3]
    qc = [q(1); -q(2); -q(3); -q(4)];
end

function qp = quat_multiply(a, b)
    % Hamilton product, scalar-first convention
    qp = [a(1)*b(1) - a(2)*b(2) - a(3)*b(3) - a(4)*b(4);
          a(1)*b(2) + a(2)*b(1) + a(3)*b(4) - a(4)*b(3);
          a(1)*b(3) - a(2)*b(4) + a(3)*b(1) + a(4)*b(2);
          a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + a(4)*b(1)];
end
