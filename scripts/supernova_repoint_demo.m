projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
 

if ~exist('demoSeed', 'var')
    demoSeed = [];
end

if ~exist('generateVisuals', 'var')
    generateVisuals = usejava('desktop');
end

clearvars -except demoSeed generateVisuals
clc;

if isempty(demoSeed)
    rng('shuffle', 'twister');
    demoSeedText = 'shuffle';
elseif (ischar(demoSeed) || (isstring(demoSeed) && isscalar(demoSeed))) && strcmpi(strtrim(char(demoSeed)), 'shuffle')
    rng('shuffle', 'twister');
    demoSeedText = 'shuffle';
elseif isscalar(demoSeed) && isnumeric(demoSeed) && isfinite(demoSeed)
    rng(double(demoSeed), 'twister');
    demoSeedText = num2str(double(demoSeed));
else
    error('supernova_repoint_demo:InvalidDemoSeed', ...
        'demoSeed must be empty, ''shuffle'', or a finite numeric scalar.');
end

vehicle = vleo.dynamics.default_vehicle_config();
cameraBody = [0; 0; 1];

initialBrightness = 1.0;
brightnessThreshold = 0.25;
brightnessDecayTimeSec = 42;
safetyMarginSec = 6;
captureThresholdDeg = 0.3;

latestUsableTimeSec = brightnessDecayTimeSec * log(initialBrightness / brightnessThreshold);
maneuverDurationSec = max(12, latestUsableTimeSec - safetyMarginSec);
omegaInitialBody = zeros(3, 1);
omegaTargetBody = zeros(3, 1);

odeOpts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

initialTargetEci = random_unit_vector();
supernovaTargetEci = random_unit_vector();
qInitial = camera_pointing_quaternion(initialTargetEci);
qTarget = camera_pointing_quaternion(supernovaTargetEci);

result = vleo.control.solve_minimax_attitude_slew_fast(qInitial, omegaInitialBody, qTarget, ...
    omegaTargetBody, vehicle.inertiaBody, maneuverDurationSec, 40, ...
    'refineNonlinear', true, 'nRefinementIntervals', [12, 16, 20]);
[tReplay, xReplay] = simulate_replay(result, qInitial, omegaInitialBody, vehicle.inertiaBody, odeOpts);
actualPointingEci = pointing_history_from_quaternions(xReplay(:, 1:4), cameraBody);
actualPointingErrorDeg = pointing_error_history_deg(actualPointingEci, supernovaTargetEci);
captureIdx = find(actualPointingErrorDeg <= captureThresholdDeg, 1, 'first');
if isempty(captureIdx)
    captureTimeSec = NaN;
else
    captureTimeSec = tReplay(captureIdx);
end
finalAttitudeErrorDeg = vleo.util.rotation_error_angle_deg(xReplay(end, 1:4)', qTarget);
finalRateErrorDegPerSec = norm(rad2deg(xReplay(end, 5:7)' - omegaTargetBody));

plannedPointingEci = pointing_history_from_quaternions(result.qVerify, cameraBody);
plannedPointingErrorDeg = pointing_error_history_deg(plannedPointingEci, supernovaTargetEci);
brightnessHistory = initialBrightness * exp(-tReplay / brightnessDecayTimeSec);

[raInitialDeg, decInitialDeg] = direction_to_ra_dec_deg(initialTargetEci);
[raTargetDeg, decTargetDeg] = direction_to_ra_dec_deg(supernovaTargetEci);
separationDeg = rad2deg(acos(vleo.util.clamp_scalar(dot(initialTargetEci, supernovaTargetEci), -1, 1)));

fprintf('=== SUPERNOVA REPOINT DEMO ===\n');
fprintf('Demo seed: %s\n', demoSeedText);
fprintf('Initial and event directions are sampled independently and uniformly on the celestial sphere.\n');
fprintf('Initial target RA/Dec: [%.2f, %.2f] deg\n', raInitialDeg, decInitialDeg);
fprintf('Supernova target RA/Dec: [%.2f, %.2f] deg\n', raTargetDeg, decTargetDeg);
fprintf('Angular separation: %.2f deg\n', separationDeg);
fprintf('Brightness threshold time: %.2f s\n', latestUsableTimeSec);
fprintf('Planned maneuver duration: %.2f s\n', maneuverDurationSec);
fprintf('Solver runtime: %.3f ms\n', 1e3 * result.solveTimeSec);
fprintf('Peak component torque gamma: %.3e N m\n', result.gamma);
fprintf('Planned terminal attitude error: %.3e deg\n', result.terminalAngleErrorDeg);
fprintf('Nonlinear replay terminal attitude error: %.3e deg\n', finalAttitudeErrorDeg);
fprintf('Nonlinear replay terminal rate error: %.3e deg/s\n', finalRateErrorDegPerSec);

if isnan(captureTimeSec)
    fprintf('Pointing threshold of %.2f deg was not reached before the event faded.\n', captureThresholdDeg);
else
    fprintf('Pointing threshold of %.2f deg reached at %.2f s.\n', captureThresholdDeg, captureTimeSec);
    fprintf('Brightness at capture: %.3f\n', initialBrightness * exp(-captureTimeSec / brightnessDecayTimeSec));
end

if ~isnan(captureTimeSec) && captureTimeSec <= latestUsableTimeSec
    fprintf('Result: repoint completes before the supernova becomes too dim.\n');
else
    fprintf('Result: repoint misses the usable brightness window.\n');
end

if generateVisuals
    figure('Name', 'Supernova Repoint Demo', 'Color', 'w');

    subplot(2, 2, 1);
    plot(tReplay, brightnessHistory, 'LineWidth', 1.6, 'Color', [0.85, 0.35, 0.15]);
    hold on;
    yline(brightnessThreshold, '--k', 'Threshold', 'LineWidth', 1.0);
    xline(latestUsableTimeSec, '--', 'Latest usable time', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.0);
    if ~isnan(captureTimeSec)
        xline(captureTimeSec, '-', 'Capture', 'Color', [0.1, 0.55, 0.2], 'LineWidth', 1.0);
    end
    grid on;
    xlabel('Time [s]');
    ylabel('Relative brightness');
    title('Supernova brightness window');

    subplot(2, 2, 2);
    plot(result.timeVerify, plannedPointingErrorDeg, '--', 'LineWidth', 1.2, 'Color', [0.2, 0.4, 0.75], ...
        'DisplayName', 'planned');
    hold on;
    plot(tReplay, actualPointingErrorDeg, 'LineWidth', 1.6, 'Color', [0.05, 0.05, 0.1], ...
        'DisplayName', 'nonlinear replay');
    yline(captureThresholdDeg, ':', 'Capture threshold', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    ylabel('Pointing error [deg]');
    title('Camera repoint performance');
    legend('Location', 'best');

    subplot(2, 2, 3);
    plot(result.timeVerify, result.tauVerify(:, 1), 'LineWidth', 1.3, 'DisplayName', 'tau_x');
    hold on;
    plot(result.timeVerify, result.tauVerify(:, 2), 'LineWidth', 1.3, 'DisplayName', 'tau_y');
    plot(result.timeVerify, result.tauVerify(:, 3), 'LineWidth', 1.3, 'DisplayName', 'tau_z');
    yline(result.gamma, '--k', 'gamma', 'LineWidth', 1.0);
    yline(-result.gamma, '--k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    grid on;
    xlabel('Time [s]');
    ylabel('Torque [N m]');
    title('Certified minimax torque command');
    legend('Location', 'best');

    subplot(2, 2, 4);
    plot(result.timeVerify, rad2deg(result.omegaVerifyBody(:, 1)), '--', 'LineWidth', 1.1, 'Color', [0.2, 0.4, 0.75], ...
        'HandleVisibility', 'off');
    hold on;
    plot(result.timeVerify, rad2deg(result.omegaVerifyBody(:, 2)), '--', 'LineWidth', 1.1, 'Color', [0.25, 0.6, 0.35], ...
        'HandleVisibility', 'off');
    plot(result.timeVerify, rad2deg(result.omegaVerifyBody(:, 3)), '--', 'LineWidth', 1.1, 'Color', [0.85, 0.35, 0.2], ...
        'HandleVisibility', 'off');
    plot(tReplay, rad2deg(xReplay(:, 5)), 'LineWidth', 1.3, 'Color', [0.2, 0.4, 0.75], 'DisplayName', 'omega_x');
    plot(tReplay, rad2deg(xReplay(:, 6)), 'LineWidth', 1.3, 'Color', [0.25, 0.6, 0.35], 'DisplayName', 'omega_y');
    plot(tReplay, rad2deg(xReplay(:, 7)), 'LineWidth', 1.3, 'Color', [0.85, 0.35, 0.2], 'DisplayName', 'omega_z');
    grid on;
    xlabel('Time [s]');
    ylabel('Body rate [deg/s]');
    title('Nonlinear body-rate replay');
    legend('Location', 'best');
end

function xDot = rigid_body_attitude_dynamics(t, x, torqueInterp, inertiaBody)
    q = x(1:4);
    q = q / norm(q);
    omegaBody = x(5:7);
    tauBody = torqueInterp(t);

    xDot = zeros(7, 1);
    xDot(1:4) = quaternion_kinematics(q, omegaBody);
    xDot(5:7) = inertiaBody \ (tauBody - cross(omegaBody, inertiaBody * omegaBody));
end

function [tReplay, xReplay] = simulate_replay(result, qInitial, omegaInitialBody, inertiaBody, odeOpts)
    torqueInterp = @(t) interp1(result.timeNodes, result.tauNodes, ...
        min(max(t, result.timeNodes(1)), result.timeNodes(end)), 'previous')';
    [tReplay, xReplay] = ode45(@(t, x) rigid_body_attitude_dynamics(t, x, torqueInterp, inertiaBody), ...
        result.timeVerify, [qInitial; omegaInitialBody], odeOpts);
    xReplay(:, 1:4) = vleo.util.normalize_quaternion_rows(xReplay(:, 1:4));
    xReplay(:, 1:4) = vleo.util.align_quaternion_signs(xReplay(:, 1:4));
end

function qDot = quaternion_kinematics(qScalarFirst, omegaBody)
    omegaMatrix = [0, -omegaBody(1), -omegaBody(2), -omegaBody(3); ...
        omegaBody(1), 0, omegaBody(3), -omegaBody(2); ...
        omegaBody(2), -omegaBody(3), 0, omegaBody(1); ...
        omegaBody(3), omegaBody(2), -omegaBody(1), 0];
    qDot = 0.5 * omegaMatrix * qScalarFirst;
end

function q = camera_pointing_quaternion(pointingEci)
    pointingEci = pointingEci / norm(pointingEci);
    referenceEci = [0; 0; 1];
    if abs(dot(referenceEci, pointingEci)) > 0.95
        referenceEci = [1; 0; 0];
    end

    xBodyInEci = referenceEci - dot(referenceEci, pointingEci) * pointingEci;
    xBodyInEci = xBodyInEci / norm(xBodyInEci);
    yBodyInEci = cross(pointingEci, xBodyInEci);
    yBodyInEci = yBodyInEci / norm(yBodyInEci);
    zBodyInEci = pointingEci;

    rBodyToEci = [xBodyInEci, yBodyInEci, zBodyInEci];
    q = dcm2quat(rBodyToEci')';
    q = q / norm(q);
    if q(1) < 0
        q = -q;
    end
end

function pointingHistory = pointing_history_from_quaternions(qHistory, cameraBody)
    nRows = size(qHistory, 1);
    pointingHistory = zeros(nRows, 3);
    for rowIdx = 1:nRows
        rBodyToEci = quat2dcm(qHistory(rowIdx, :))';
        pointingVector = rBodyToEci * cameraBody;
        pointingHistory(rowIdx, :) = (pointingVector / norm(pointingVector))';
    end
end

function errorDeg = pointing_error_history_deg(pointingHistory, targetDirectionEci)
    targetDirectionEci = targetDirectionEci / norm(targetDirectionEci);
    cosAngle = pointingHistory * targetDirectionEci;
    cosAngle = max(min(cosAngle, 1), -1);
    errorDeg = rad2deg(acos(cosAngle));
end

function [raDeg, decDeg] = direction_to_ra_dec_deg(directionEci)
    directionEci = directionEci / norm(directionEci);
    raRad = atan2(directionEci(2), directionEci(1));
    if raRad < 0
        raRad = raRad + 2 * pi;
    end
    decRad = asin(directionEci(3));
    raDeg = rad2deg(raRad);
    decDeg = rad2deg(decRad);
end

function unitVector = random_unit_vector()
    phi = 2 * pi * rand();
    z = 2 * rand() - 1;
    radial = sqrt(max(0, 1 - z^2));
    unitVector = [radial * cos(phi); radial * sin(phi); z];
end
