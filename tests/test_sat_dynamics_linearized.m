function results = test_sat_dynamics_linearized()
% TEST_SAT_DYNAMICS_LINEARIZED Smoke tests for the small-angle linearized attitude model.

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    results = struct('passed', 0, 'failed', 0, 'tests', {{}});
    testFuncs = {@test_trim_state_has_zero_rotational_dynamics, ...
        @test_small_attitude_error_generates_restoring_torque, ...
        @test_control_disable_zeros_linearized_torque};

    for idx = 1:numel(testFuncs)
        results = run_test(results, testFuncs{idx});
    end

    print_summary(results);
end

function passed = test_trim_state_has_zero_rotational_dynamics()
    simParams = vleo.gui.default_simulation_params();
    [X0, Xr] = vleo.dynamics.initial_state_from_sim_params(simParams);
    Xtrim = X0;
    Xtrim(7:10) = Xr(1:4);
    Xtrim(11:13) = Xr(5:7);

    Xd = vleo.dynamics.sat_dynamics_linearized(0, Xtrim, Xr, simParams);
    passed = verify_vector(Xd(7:13), zeros(7, 1), 1e-12, 'Trim rotational derivative');
end

function passed = test_small_attitude_error_generates_restoring_torque()
    simParams = vleo.gui.default_simulation_params();
    [X0, Xr] = vleo.dynamics.initial_state_from_sim_params(simParams);
    Xperturbed = X0;
    Xperturbed(7:10) = normalize_quaternion(Xr(1:4) + [0; 1e-2; 0; 0]);
    tauBody = vleo.dynamics.linearized_control_torque(Xperturbed, Xr, simParams);
    Xd = vleo.dynamics.sat_dynamics_linearized(0, Xperturbed, Xr, simParams);

    passed = verify_true(tauBody(1) < 0, 'Positive x-axis attitude error should produce negative restoring torque.') & ...
        verify_true(Xd(11) < 0, 'Positive x-axis attitude error should produce negative x angular acceleration.');
end

function passed = test_control_disable_zeros_linearized_torque()
    simParams = vleo.gui.default_simulation_params();
    simParams.initParams.Modes.enableControl = false;
    [X0, Xr] = vleo.dynamics.initial_state_from_sim_params(simParams);
    X0(7:10) = normalize_quaternion(Xr(1:4) + [0; 5e-3; -3e-3; 2e-3]);
    X0(11:13) = [0.01; -0.02; 0.03];

    tauBody = vleo.dynamics.linearized_control_torque(X0, Xr, simParams);
    passed = verify_vector(tauBody, zeros(3, 1), 1e-12, 'Control-disabled linearized torque');
end

function q = normalize_quaternion(q)
    q = reshape(q, [], 1);
    q = q / norm(q);
    if q(1) < 0
        q = -q;
    end
end

function passed = verify_vector(actual, expected, tolerance, label)
    passed = numel(actual) == numel(expected) && all(abs(actual(:) - expected(:)) <= tolerance);
    if ~passed
        fprintf('    %s mismatch.\n', label);
        fprintf('    Expected: [%s]\n', num2str(expected(:).', ' %.16g'));
        fprintf('    Actual:   [%s]\n', num2str(actual(:).', ' %.16g'));
    end
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
