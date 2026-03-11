%% Unit Test Template
% Copy this template when creating new tests.
% Save as: test_<function_name>.m
%
% Author: [Your Name]
% Date: YYYY-MM-DD

function results = test_template()
%TEST_TEMPLATE Unit tests for [function_name]
%   results = test_template() runs all tests and returns results structure
%
%   Test naming convention:
%       test_<function_name>_<test_description>

    % Initialize test framework
    results = struct('passed', 0, 'failed', 0, 'tests', {{}});

    % Add project paths
    project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(project_root);
    setup_project();

    % Run tests
    results = run_test(results, @test_basic_functionality);
    results = run_test(results, @test_edge_case_1);
    results = run_test(results, @test_edge_case_2);
    results = run_test(results, @test_against_reference);

    % Print summary
    print_summary(results);
end

%% Test Functions
% Each test function should:
%   - Take no inputs
%   - Return true if test passes, false if fails
%   - Have descriptive name explaining what is tested

function passed = test_basic_functionality()
    % Test basic/nominal case
    % Description: [What this test verifies]

    % Setup
    input_value = 1;  % Replace with actual test input
    expected_output = 1;  % Replace with expected output
    tolerance = 1e-10;

    % Execute
    actual_output = input_value;  % Replace with actual function call

    % Verify
    passed = abs(actual_output - expected_output) < tolerance;

    if ~passed
        fprintf('    Expected: %g, Got: %g\n', expected_output, actual_output);
    end
end

function passed = test_edge_case_1()
    % Test edge case: [describe edge case]
    % Example: Zero input, singular conditions, limits

    % Setup
    input_value = 0;
    expected_output = 0;
    tolerance = 1e-10;

    % Execute
    actual_output = input_value;  % Replace with actual function call

    % Verify
    passed = abs(actual_output - expected_output) < tolerance;
end

function passed = test_edge_case_2()
    % Test edge case: [describe another edge case]

    passed = true;  % Replace with actual test
end

function passed = test_against_reference()
    % Test against known reference values (e.g., from textbook, paper)
    % Source: [Citation]

    % Reference values from [Source]
    % Example: Vallado, Table 2-4, p. 113

    passed = true;  % Replace with actual test
end

%% Helper Functions

function results = run_test(results, test_func)
    % Run a single test and record results
    test_name = func2str(test_func);

    try
        passed = test_func();
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

    % Record test
    results.tests{end+1} = struct('name', test_name, 'status', status);
    fprintf('[%s] %s\n', status, test_name);
end

function print_summary(results)
    % Print test summary
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
