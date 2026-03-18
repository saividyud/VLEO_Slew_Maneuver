function results = test_too_minimax_slew_exported_cases()
% TEST_TOO_MINIMAX_SLEW_EXPORTED_CASES Data-driven regression using exported TOO minimax slew cases.

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    results = struct('passed', 0, 'failed', 0, 'tests', {{}});
    results = run_test(results, @test_exported_too_minimax_slew_cases);
    print_summary(results);
end

function passed = test_exported_too_minimax_slew_cases()
    caseBundlePath = fullfile(vleo.util.generated_data_dir(), 'too_minimax_slew_solver_cases.mat');
    if ~exist(caseBundlePath, 'file')
        error('test_too_minimax_slew_exported_cases:MissingCaseBundle', ...
            'Missing exported case bundle. Run scripts/export_too_minimax_slew_cases.m first.');
    end

    data = load(caseBundlePath, 'cases', 'metadata');
    if ~isfield(data, 'cases') || isempty(data.cases)
        error('test_too_minimax_slew_exported_cases:EmptyCaseBundle', ...
            'The exported case bundle does not contain any solver cases.');
    end

    vehicle = vleo.dynamics.default_vehicle_config();
    angleToleranceDeg = 5e-2;
    rateToleranceDegPerSec = 5e-3;
    passMask = true(numel(data.cases), 1);

    fprintf('Loaded %d exported TOO minimax slew cases from %s\n', numel(data.cases), caseBundlePath);
    if isfield(data, 'metadata') && isfield(data.metadata, 'createdAt')
        fprintf('Case bundle created at: %s\n', data.metadata.createdAt);
    end

    for caseIdx = 1:numel(data.cases)
        solverCase = data.cases(caseIdx);
        result = vleo.control.solve_minimax_attitude_slew_fast( ...
            solverCase.qInitial, solverCase.omegaInitialBody, solverCase.qTarget, ...
            solverCase.omegaTargetBody, vehicle.inertiaBody, solverCase.maneuverDurationSec, 20, ...
            'refineNonlinear', true, 'nRefinementIntervals', [12, 16, 20]);

        replay = replay_piecewise_torque(result, solverCase.qInitial, solverCase.omegaInitialBody, vehicle.inertiaBody);
        replayAngleErrorDeg = vleo.util.rotation_error_angle_deg(replay.qEnd, solverCase.qTarget);
        replayRateErrorDegPerSec = norm(rad2deg(replay.omegaEnd - solverCase.omegaTargetBody));

        passMask(caseIdx) = replayAngleErrorDeg <= angleToleranceDeg && ...
            replayRateErrorDegPerSec <= rateToleranceDegPerSec;

        if passMask(caseIdx)
            fprintf('[PASS] %s  angle=%.3e deg  rate=%.3e deg/s  gamma=%.3e N m\n', ...
                solverCase.caseId, replayAngleErrorDeg, replayRateErrorDegPerSec, result.gamma);
        else
            fprintf('[FAIL] %s  exit=%d  angle=%.3e deg  rate=%.3e deg/s  gamma=%.3e N m\n', ...
                solverCase.caseId, result.exitflag, replayAngleErrorDeg, replayRateErrorDegPerSec, result.gamma);
        end
    end

    passed = all(passMask);
end

function replay = replay_piecewise_torque(result, qInitial, omegaInitialBody, inertiaBody)
    odeOpts = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
    torqueInterp = @(t) interp1(result.timeNodes, result.tauNodes, ...
        min(max(t, result.timeNodes(1)), result.timeNodes(end)), 'previous')';
    [~, xReplay] = ode45(@(t, x) rigid_body_attitude_dynamics(t, x, torqueInterp, inertiaBody), ...
        result.timeNodes, [qInitial; omegaInitialBody], odeOpts);
    xReplay(:, 1:4) = vleo.util.normalize_quaternion_rows(xReplay(:, 1:4));
    xReplay(:, 1:4) = vleo.util.align_quaternion_signs(xReplay(:, 1:4));

    replay = struct();
    replay.qEnd = xReplay(end, 1:4)';
    replay.omegaEnd = xReplay(end, 5:7)';
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

function qDot = quaternion_kinematics(qScalarFirst, omegaBody)
    omegaMatrix = [0, -omegaBody(1), -omegaBody(2), -omegaBody(3); ...
        omegaBody(1), 0, omegaBody(3), -omegaBody(2); ...
        omegaBody(2), -omegaBody(3), 0, omegaBody(1); ...
        omegaBody(3), omegaBody(2), -omegaBody(1), 0];
    qDot = 0.5 * omegaMatrix * qScalarFirst;
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
