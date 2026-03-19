% Force-Limited Minimum Slew Time vs Angular Separation
% Computes the minimum rest-to-rest maneuver duration for angular
% separations from 0 to 180 degrees, subject to a peak actuator force
% of 145 mN with moment arms [0.3, 0.3, 0.2] m.
%
% For a rest-to-rest minimax rotation of angle theta about principal
% axis k: T_min = 2 * sqrt(theta * I_k / (L_k * F_max))

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc; close all;

vehicle = vleo.dynamics.default_vehicle_config();
momentArms = [0.3; 0.3; 0.2];
forceLimit_N = 0.145;
principalInertia = diag(vehicle.inertiaBody);
ratioIL = principalInertia ./ momentArms;

fprintf('=== FORCE-LIMITED MINIMUM SLEW TIME ===\n');
fprintf('Vehicle: 6U CubeSat (%.0f kg)\n', vehicle.mass);
fprintf('Principal inertias [kg m^2]: [%.4f, %.4f, %.4f]\n', principalInertia);
fprintf('Moment arms [m]: [%.2f, %.2f, %.2f]\n', momentArms);
fprintf('Force limit: %.1f mN\n', forceLimit_N * 1e3);
fprintf('I/L ratios [kg m]: [%.4f, %.4f, %.4f]\n\n', ratioIL);

%% Analytic minimum time for rest-to-rest principal-axis rotations
angles_deg = linspace(0.5, 180, 360);
angles_rad = deg2rad(angles_deg);

T_min = zeros(numel(angles_deg), 3);
for k = 1:3
    T_min(:, k) = 2 * sqrt(angles_rad(:) * ratioIL(k) / forceLimit_N);
end

%% Optimal rotation axis (equalizes force across all axes)
% n_k = gamma * L_k / I_k, gamma = 1 / sqrt(sum((L_k/I_k)^2))
sumLI2 = sum((momentArms ./ principalInertia).^2);
gammaOpt = 1 / sqrt(sumLI2);
nOpt = gammaOpt * (momentArms ./ principalInertia);
T_opt = 2 * sqrt(angles_rad(:) * gammaOpt / forceLimit_N);

fprintf('Optimal rotation axis (equalizes force): [%.3f, %.3f, %.3f]\n', nOpt);
fprintf('Optimal I/L ratio: %.4f kg m  (vs worst %.4f, best %.4f)\n\n', gammaOpt, max(ratioIL), min(ratioIL));

%% Verify a few points with the solver
verifyAngles_deg = [30, 90, 180];
T_ref = 30;
qIdentity = [1; 0; 0; 0];
omegaZero = zeros(3, 1);

fprintf('Verification (analytic vs solver):\n');
for vIdx = 1:numel(verifyAngles_deg)
    theta_rad = deg2rad(verifyAngles_deg(vIdx));
    for k = 1:3
        axisVec = zeros(3, 1); axisVec(k) = 1;
        qTarget = rotation_vector_to_quaternion(theta_rad * axisVec);
        result = vleo.control.solve_minimax_attitude_slew_fast(qIdentity, omegaZero, qTarget, ...
            omegaZero, vehicle.inertiaBody, T_ref, 20, 'momentArms', momentArms);
        T_solver = T_ref * sqrt(result.gamma / forceLimit_N);
        T_analytic = 2 * sqrt(theta_rad * ratioIL(k) / forceLimit_N);
        fprintf('  %3d deg, axis %d: T=%.3f s, solver=%.3f s, diff=%.2e s\n', ...
            verifyAngles_deg(vIdx), k, T_analytic, T_solver, abs(T_analytic - T_solver));
    end

    qTarget = rotation_vector_to_quaternion(theta_rad * nOpt);
    result = vleo.control.solve_minimax_attitude_slew_fast(qIdentity, omegaZero, qTarget, ...
        omegaZero, vehicle.inertiaBody, T_ref, 20, 'momentArms', momentArms);
    T_solver = T_ref * sqrt(result.gamma / forceLimit_N);
    T_analytic = 2 * sqrt(theta_rad * gammaOpt / forceLimit_N);
    fprintf('  %3d deg, optimal: T=%.3f s, solver=%.3f s, diff=%.2e s\n', ...
        verifyAngles_deg(vIdx), T_analytic, T_solver, abs(T_analytic - T_solver));
end

%% Reference points
fprintf('\nReference slew times:\n');
for angleDeg = [45, 90, 135, 180]
    theta_rad = deg2rad(angleDeg);
    T_worst = 2 * sqrt(theta_rad * max(ratioIL) / forceLimit_N);
    T_best = 2 * sqrt(theta_rad * gammaOpt / forceLimit_N);
    fprintf('  %3d deg: worst-case %.2f s, optimal-axis %.2f s\n', angleDeg, T_worst, T_best);
end

%% Plot
figure('Name', 'Min Slew Time vs Angular Separation', ...
    'Color', 'w', 'Position', [100, 100, 1000, 600]);

colors = [0.85, 0.15, 0.15; 0.15, 0.45, 0.85; 0.15, 0.65, 0.25];
axisNames = {'x', 'y', 'z'};

hold on;
for k = 1:3
    plot(angles_deg, T_min(:, k), '-', 'LineWidth', 2, 'Color', colors(k, :), ...
        'DisplayName', sprintf('%s-axis (I/L = %.3f)', axisNames{k}, ratioIL(k)));
end
plot(angles_deg, T_opt, 'k--', 'LineWidth', 2.2, ...
    'DisplayName', sprintf('optimal axis (I/L_{eff} = %.3f)', gammaOpt));
hold off;

grid on;
xlabel('Angular Separation (deg)', 'FontSize', 16);
ylabel('Minimum Slew Time (s)', 'FontSize', 16);
title(sprintf('Rest-to-Rest Minimax Slew Time (F_{max} = %g mN)', forceLimit_N * 1e3), 'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 13);
xlim([0, 180]);
ylim([0, inf]);
set(gca, 'FontSize', 13, 'LineWidth', 1.0);

hold on;
referenceAnglesDeg = [45, 90, 135, 180];
referenceAnglesRad = deg2rad(referenceAnglesDeg);
worstAxisIdx = find(ratioIL == max(ratioIL), 1, 'first');
worstAxisTimes = 2 * sqrt(referenceAnglesRad(:) * ratioIL(worstAxisIdx) / forceLimit_N);
optimalAxisTimes = 2 * sqrt(referenceAnglesRad(:) * gammaOpt / forceLimit_N);
plot(referenceAnglesDeg, worstAxisTimes, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('reference points (%s-axis)', axisNames{worstAxisIdx}));
plot(referenceAnglesDeg, optimalAxisTimes, 'kd', 'MarkerSize', 6, 'MarkerFaceColor', [0.3, 0.3, 0.3], ...
    'DisplayName', 'reference points (optimal axis)');
hold off;

annotationText = sprintf([ ...
    'Worst-axis 180 deg: %.2f s\n', ...
    'Optimal-axis 180 deg: %.2f s\n', ...
    '145 mN force limit, L = [0.3, 0.3, 0.2] m'], ...
    worstAxisTimes(end), optimalAxisTimes(end));
text(0.58, 0.23, annotationText, 'Units', 'normalized', 'FontSize', 12, ...
    'BackgroundColor', 'w', 'Margin', 8);

figureExport = vleo.util.export_paper_figure(gcf, 'force_limited_slew_sweep');
fprintf('Paper figure saved to: %s\n', figureExport.pngPath);
fprintf('Paper figure saved to: %s\n', figureExport.pdfPath);

fprintf('\nDone.\n');

function q = rotation_vector_to_quaternion(rotationVector)
    r = rotation_vector_to_dcm(rotationVector);
    q = dcm2quat(r)';
    q = q / norm(q);
    if q(1) < 0
        q = -q;
    end
end

function r = rotation_vector_to_dcm(rotationVector)
    angle = norm(rotationVector);
    kx = skew_matrix(rotationVector);
    if angle < 1e-12
        r = eye(3) + kx + 0.5 * (kx * kx);
        return;
    end
    axis = rotationVector / angle;
    k = skew_matrix(axis);
    r = eye(3) + sin(angle) * k + (1 - cos(angle)) * (k * k);
end

function s = skew_matrix(v)
    s = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end
