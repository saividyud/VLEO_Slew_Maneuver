projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('benchmarkSeed', 'var')
    benchmarkSeed = 17;
end

if ~exist('nBenchmarkCases', 'var')
    nBenchmarkCases = 200;
end

if ~exist('generateVisuals', 'var')
    generateVisuals = usejava('desktop');
end

clearvars -except benchmarkSeed nBenchmarkCases generateVisuals
clc;

rng(benchmarkSeed, 'twister');
vehicle = vleo.dynamics.default_vehicle_config();
principalInertia = diag(vehicle.inertiaBody);
nVerifyCases = min(12, nBenchmarkCases);

solveTimesSec = zeros(nBenchmarkCases, 1);
rotationAnglesDeg = zeros(nBenchmarkCases, 1);
maneuverDurationsSec = zeros(nBenchmarkCases, 1);
peakTorqueAbs = zeros(nBenchmarkCases, 1);
terminalAngleErrorsDeg = zeros(nBenchmarkCases, 1);
terminalRateErrorsDegPerSec = zeros(nBenchmarkCases, 1);
denseVerifierGamma = nan(nVerifyCases, 1);
denseVerifierCaseIdx = (1:nVerifyCases)';

fprintf('=== FAST MINIMAX ATTITUDE SLEW BENCHMARK ===\n');
fprintf('Benchmark seed: %d\n', benchmarkSeed);
fprintf('Cases: %d\n', nBenchmarkCases);
fprintf('Dense verification subset: %d\n', nVerifyCases);
fprintf('Principal inertias [kg m^2]: [%.4f, %.4f, %.4f]\n\n', principalInertia);

for caseIdx = 1:nBenchmarkCases
    [qInitial, qTarget, omegaInitialBody, omegaTargetBody, maneuverDuration, rotationAngleDeg] = random_case();

    result = vleo.control.solve_minimax_attitude_slew_fast(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, vehicle.inertiaBody, maneuverDuration, 24);

    solveTimesSec(caseIdx) = result.solveTimeSec;
    rotationAnglesDeg(caseIdx) = rotationAngleDeg;
    maneuverDurationsSec(caseIdx) = maneuverDuration;
    peakTorqueAbs(caseIdx) = result.gamma;
    terminalAngleErrorsDeg(caseIdx) = result.terminalAngleErrorDeg;
    terminalRateErrorsDegPerSec(caseIdx) = result.terminalRateErrorDegPerSec;

    if caseIdx <= nVerifyCases
        denseGamma = 0;
        for axisIdx = 1:3
            denseGamma = max(denseGamma, dense_axis_peak_torque( ...
                result.relativeRotationVectorBody(axisIdx), omegaInitialBody(axisIdx), ...
                omegaTargetBody(axisIdx), principalInertia(axisIdx), maneuverDuration));
        end
        denseVerifierGamma(caseIdx) = denseGamma;
    end
end

medianSolveMs = 1e3 * median(solveTimesSec);
percentile95SolveMs = 1e3 * simple_percentile(solveTimesSec, 95);
maxSolveMs = 1e3 * max(solveTimesSec);
maxTerminalAngleErrorDeg = max(terminalAngleErrorsDeg);
maxTerminalRateErrorDegPerSec = max(terminalRateErrorsDegPerSec);
maxDenseMismatch = max(abs(peakTorqueAbs(denseVerifierCaseIdx) - denseVerifierGamma));

fprintf('Median solve time: %.3f ms\n', medianSolveMs);
fprintf('95th percentile solve time: %.3f ms\n', percentile95SolveMs);
fprintf('Worst solve time: %.3f ms\n', maxSolveMs);
fprintf('Max terminal attitude error in solver model: %.3e deg\n', maxTerminalAngleErrorDeg);
fprintf('Max terminal rate error in solver model: %.3e deg/s\n', maxTerminalRateErrorDegPerSec);
fprintf('Max dense-verifier gamma mismatch: %.3e N m\n', maxDenseMismatch);
fprintf('Peak commanded torque range: [%.3e, %.3e] N m\n', min(peakTorqueAbs), max(peakTorqueAbs));
fprintf('Rotation-angle range: [%.2f, %.2f] deg\n', min(rotationAnglesDeg), max(rotationAnglesDeg));
fprintf('Maneuver-duration range: [%.2f, %.2f] s\n', min(maneuverDurationsSec), max(maneuverDurationsSec));

if generateVisuals
    figure('Name', 'Fast Minimax Solver Benchmark', 'Color', 'w');

    subplot(2, 2, 1);
    scatter(rotationAnglesDeg, 1e3 * solveTimesSec, 24, maneuverDurationsSec, 'filled');
    grid on;
    xlabel('Relative rotation angle [deg]');
    ylabel('Solve time [ms]');
    title('Solve time vs maneuver angle');
    colorbar;

    subplot(2, 2, 2);
    plot(sort(1e3 * solveTimesSec), 'LineWidth', 1.5);
    grid on;
    xlabel('Sorted case index');
    ylabel('Solve time [ms]');
    title('Solve-time distribution');

    subplot(2, 2, 3);
    plot(denseVerifierCaseIdx, peakTorqueAbs(denseVerifierCaseIdx), 'o-', 'LineWidth', 1.2, ...
        'DisplayName', 'analytic gamma');
    hold on;
    plot(denseVerifierCaseIdx, denseVerifierGamma, 'x--', 'LineWidth', 1.2, ...
        'DisplayName', 'dense verifier');
    grid on;
    xlabel('Verification case');
    ylabel('Peak torque [N m]');
    title('Analytic gamma vs dense verifier');
    legend('Location', 'best');

    subplot(2, 2, 4);
    scatter(maneuverDurationsSec, peakTorqueAbs, 24, rotationAnglesDeg, 'filled');
    grid on;
    xlabel('Maneuver duration [s]');
    ylabel('Peak component torque [N m]');
    title('Torque demand across benchmark set');
    colorbar;
end

function [qInitial, qTarget, omegaInitialBody, omegaTargetBody, maneuverDuration, rotationAngleDeg] = random_case()
    qInitial = random_quaternion();

    rotationAngleDeg = 4 + 28 * rand();
    relativeRotationVector = deg2rad(rotationAngleDeg) * random_unit_vector();
    qTarget = apply_relative_rotation(qInitial, relativeRotationVector);

    omegaInitialBody = deg2rad(0.8 * (2 * rand(3, 1) - 1));
    omegaTargetBody = deg2rad(0.8 * (2 * rand(3, 1) - 1));
    maneuverDuration = 18 + 72 * rand();
end

function q = random_quaternion()
    angle = 2 * pi * rand();
    q = rotation_vector_to_quaternion(angle * random_unit_vector());
end

function unitVector = random_unit_vector()
    unitVector = randn(3, 1);
    unitVector = unitVector / norm(unitVector);
end

function qTarget = apply_relative_rotation(qInitial, relativeRotationVector)
    rInitial = quat2dcm(qInitial');
    rRelative = rotation_vector_to_dcm(relativeRotationVector);
    qTarget = dcm2quat(rRelative * rInitial)';
    qTarget = qTarget / norm(qTarget);
    if qTarget(1) < 0
        qTarget = -qTarget;
    end
end

function peakTorque = dense_axis_peak_torque(thetaTarget, omegaInitial, omegaTarget, inertiaAxis, maneuverDuration)
    deltaRate = omegaTarget - omegaInitial;
    epsilon = max(1e-9 * maneuverDuration, 1e-12);
    timeGrid = linspace(epsilon, maneuverDuration - epsilon, 8001);
    peakTorque = inf;

    for initialSign = [-1, 1]
        residualGrid = nan(size(timeGrid));
        validGrid = false(size(timeGrid));

        for gridIdx = 1:numel(timeGrid)
            switchTime = timeGrid(gridIdx);
            denominator = initialSign * (2 * switchTime - maneuverDuration);
            if abs(denominator) <= 1e-12
                continue;
            end

            accelMag = deltaRate / denominator;
            if accelMag <= 0 || ~isfinite(accelMag)
                continue;
            end

            residualGrid(gridIdx) = axis_position_residual(switchTime, thetaTarget, ...
                omegaInitial, omegaTarget, maneuverDuration, initialSign);
            validGrid(gridIdx) = true;
        end

        exactIdx = find(validGrid & abs(residualGrid) <= 1e-10);
        for idx = exactIdx(:)'
            switchTime = timeGrid(idx);
            accelMag = acceleration_from_switch_time(switchTime, omegaInitial, omegaTarget, ...
                maneuverDuration, initialSign);
            peakTorque = min(peakTorque, inertiaAxis * accelMag);
        end

        signChangeIdx = find(validGrid(1:end - 1) & validGrid(2:end) & ...
            residualGrid(1:end - 1) .* residualGrid(2:end) < 0);
        for idx = signChangeIdx(:)'
            bracket = [timeGrid(idx), timeGrid(idx + 1)];
            root = fzero(@(switchTime) axis_position_residual(switchTime, thetaTarget, ...
                omegaInitial, omegaTarget, maneuverDuration, initialSign), bracket);
            accelMag = acceleration_from_switch_time(root, omegaInitial, omegaTarget, ...
                maneuverDuration, initialSign);
            if accelMag > 0 && isfinite(accelMag)
                peakTorque = min(peakTorque, inertiaAxis * accelMag);
            end
        end
    end

    if ~isfinite(peakTorque)
        error('benchmark_minimax_solver:DenseVerifierFailed', ...
            'Dense single-switch verifier failed to find a feasible root.');
    end
end

function residual = axis_position_residual(switchTime, thetaTarget, omegaInitial, omegaTarget, maneuverDuration, initialSign)
    accelMag = acceleration_from_switch_time(switchTime, omegaInitial, omegaTarget, maneuverDuration, initialSign);
    shapeTerm = initialSign * accelMag * (2 * switchTime * maneuverDuration - switchTime^2 - 0.5 * maneuverDuration^2);
    thetaFinal = omegaInitial * maneuverDuration + shapeTerm;
    residual = thetaFinal - thetaTarget;
end

function accelMag = acceleration_from_switch_time(switchTime, omegaInitial, omegaTarget, maneuverDuration, initialSign)
    denominator = initialSign * (2 * switchTime - maneuverDuration);
    accelMag = (omegaTarget - omegaInitial) / denominator;
end

function percentileValue = simple_percentile(values, percentile)
    values = sort(values(:));
    if isempty(values)
        percentileValue = NaN;
        return;
    end

    index = 1 + (numel(values) - 1) * percentile / 100;
    lowerIdx = floor(index);
    upperIdx = ceil(index);
    if lowerIdx == upperIdx
        percentileValue = values(lowerIdx);
        return;
    end

    weight = index - lowerIdx;
    percentileValue = (1 - weight) * values(lowerIdx) + weight * values(upperIdx);
end

function q = rotation_vector_to_quaternion(rotationVector)
    q = dcm2quat(rotation_vector_to_dcm(rotationVector))';
    q = q / norm(q);
    if q(1) < 0
        q = -q;
    end
end

function r = rotation_vector_to_dcm(rotationVector)
    angle = norm(rotationVector);
    if angle < 1e-12
        kx = skew_matrix(rotationVector);
        r = eye(3) + kx + 0.5 * (kx * kx);
        return;
    end

    axis = rotationVector / angle;
    k = skew_matrix(axis);
    r = eye(3) + sin(angle) * k + (1 - cos(angle)) * (k * k);
end

function s = skew_matrix(v)
    s = [0, -v(3), v(2); ...
         v(3), 0, -v(1); ...
         -v(2), v(1), 0];
end
