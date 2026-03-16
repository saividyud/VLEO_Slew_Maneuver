function results = test_solve_minimax_attitude_slew_fast()
% TEST_SOLVE_MINIMAX_ATTITUDE_SLEW_FAST Regression tests for the fast minimax slew solver.

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    results = struct('passed', 0, 'failed', 0, 'tests', {{}});
    testFuncs = {@test_zero_maneuver_returns_zero_torque, ...
        @test_solver_hits_requested_terminal_boundary, ...
        @test_certificate_matches_dense_single_switch_verifier};

    for idx = 1:numel(testFuncs)
        results = run_test(results, testFuncs{idx});
    end

    print_summary(results);
end

function passed = test_zero_maneuver_returns_zero_torque()
    vehicle = vleo.dynamics.default_vehicle_config();
    qIdentity = [1; 0; 0; 0];
    omegaZero = zeros(3, 1);
    result = vleo.control.solve_minimax_attitude_slew_fast( ...
        qIdentity, omegaZero, qIdentity, omegaZero, vehicle.inertiaBody, 30, 12);

    checks = { ...
        result.gamma, 0, 1e-12, 'Zero-maneuver gamma'; ...
        result.peakAxisTorqueAbs, zeros(3, 1), 1e-12, 'Zero-maneuver peak torque'; ...
        result.omegaVerifyBody(end, :)', omegaZero, 1e-12, 'Zero-maneuver terminal rate'; ...
        result.tauVerify, zeros(size(result.tauVerify)), 1e-12, 'Zero-maneuver torque history' ...
    };

    passed = run_checks(checks);
    passed = verify_true(result.isCertifiedMinimax, 'Zero maneuver should be certified minimax.') & passed;
    passed = verify_true(result.terminalAngleErrorDeg <= 1e-10, ...
        sprintf('Zero maneuver should finish at zero angle error, got %.3e deg.', result.terminalAngleErrorDeg)) & passed;
end

function passed = test_solver_hits_requested_terminal_boundary()
    vehicle = vleo.dynamics.default_vehicle_config();
    qInitial = [1; 0; 0; 0];
    qTarget = rotation_vector_to_quaternion(deg2rad([9; -6; 7]));
    omegaInitialBody = deg2rad([0.4; -0.3; 0.2]);
    omegaTargetBody = deg2rad([-0.2; 0.15; -0.1]);

    result = vleo.control.solve_minimax_attitude_slew_fast(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, vehicle.inertiaBody, 42, 18);

    checks = { ...
        result.omegaVerifyBody(end, :)', omegaTargetBody, 1e-12, 'Terminal body rate'; ...
        max(abs(result.terminalRotationVectorError)), 0, 1e-12, 'Terminal rotation-vector error'; ...
        max(abs(result.peakAxisTorqueAbs - max(abs(result.tauVerify), [], 1)')), 0, 1e-12, 'Peak torque extraction' ...
    };

    passed = run_checks(checks);
    passed = verify_true(result.terminalAngleErrorDeg <= 1e-10, ...
        sprintf('Terminal quaternion error should be near zero, got %.3e deg.', result.terminalAngleErrorDeg)) & passed;
    passed = verify_true(all(result.switchTimes >= 0 & result.switchTimes <= 42), ...
        'All switch times should stay inside the maneuver horizon.') & passed;
end

function passed = test_certificate_matches_dense_single_switch_verifier()
    vehicle = vleo.dynamics.default_vehicle_config();
    qInitial = rotation_vector_to_quaternion(deg2rad([2; -1; 1.5]));
    qTarget = rotation_vector_to_quaternion(deg2rad([10; -7; 5]));
    omegaInitialBody = deg2rad([0.55; -0.25; 0.3]);
    omegaTargetBody = deg2rad([-0.15; 0.35; -0.2]);
    maneuverDuration = 36;

    result = vleo.control.solve_minimax_attitude_slew_fast(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, vehicle.inertiaBody, maneuverDuration, 16);

    relativeRotationVector = result.relativeRotationVectorBody;
    denseGamma = 0;
    principalInertia = diag(vehicle.inertiaBody);
    for axisIdx = 1:3
        denseGamma = max(denseGamma, dense_axis_peak_torque( ...
            relativeRotationVector(axisIdx), omegaInitialBody(axisIdx), omegaTargetBody(axisIdx), ...
            principalInertia(axisIdx), maneuverDuration));
    end

    passed = verify_equal(result.gamma, denseGamma, 1e-10, 'Certified gamma vs dense verifier');
    passed = verify_true(result.isCertifiedMinimax, 'Fast solver should return a minimax certificate.') & passed;
    passed = verify_true(abs(result.optimalityCertificate.gap) <= 1e-12, ...
        'The optimality certificate should report zero gap.') & passed;
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
        error('test_solve_minimax_attitude_slew_fast:DenseVerifierFailed', ...
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

function passed = run_checks(checks)
    passed = true;
    for idx = 1:size(checks, 1)
        passed = verify_equal(checks{idx, 1}, checks{idx, 2}, checks{idx, 3}, checks{idx, 4}) & passed;
    end
end

function passed = verify_equal(actual, expected, tolerance, label)
    passed = numel(actual) == numel(expected) && all(abs(actual(:) - expected(:)) <= tolerance);
    if passed
        return;
    end

    if isscalar(actual) && isscalar(expected)
        fprintf('    %s mismatch. Expected %.16g, got %.16g\n', label, expected, actual);
        return;
    end

    fprintf('    %s mismatch.\n', label);
    fprintf('    Expected: [%s]\n', num2str(expected(:).', ' %.16g'));
    fprintf('    Actual:   [%s]\n', num2str(actual(:).', ' %.16g'));
end

function passed = verify_true(condition, message)
    passed = logical(condition);
    if ~passed
        fprintf('    %s\n', message);
    end
end

function results = run_test(results, testFunc)
    testName = func2str(testFunc);

    try
        passed = testFunc();
        if passed
            results.passed = results.passed + 1;
            status = 'PASS';
        else
            results.failed = results.failed + 1;
            status = 'FAIL';
        end
    catch ME
        results.failed = results.failed + 1;
        status = 'ERROR';
        fprintf('    Exception: %s\n', ME.message);
    end

    results.tests{end + 1} = struct('name', testName, 'status', status);
    fprintf('[%s] %s\n', status, testName);
end

function print_summary(results)
    total = results.passed + results.failed;
    fprintf('\n=== Test Summary ===\n');
    fprintf('Total:  %d\n', total);
    fprintf('Passed: %d\n', results.passed);
    fprintf('Failed: %d\n', results.failed);

    if results.failed == 0
        fprintf('\nAll tests passed!\n');
    else
        fprintf('\nSome tests failed. Review output above.\n');
    end
end
