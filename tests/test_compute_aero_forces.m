function results = test_compute_aero_forces()
% TEST_COMPUTEAEROFORCES Regression and interface tests for compute_aero_forces.

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    results = struct('passed', 0, 'failed', 0, 'tests', {{}});
    testFuncs = {@test_cook_reference_case, @test_sentman_reference_case, ...
        @test_position_eci_matches_geodetic_inputs, @test_unsupported_gsi_model_errors};

    for idx = 1:numel(testFuncs)
        results = run_test(results, testFuncs{idx});
    end

    print_summary(results);
end

function passed = test_cook_reference_case()
    actual = compute_results(nominal_inputs('cook'));
    checks = { ...
        actual.AreaRef, 0.110000005513, 1e-12, 'AreaRef'; ...
        actual.AreaProj, 0.060000003278, 1e-12, 'AreaProj'; ...
        actual.LenRef, 0.100000001490, 1e-12, 'LenRef'; ...
        actual.Cf_w, [-1.146797966946; 0; 0], 1e-10, 'Cf_w'; ...
        actual.Cm_B, [0; 1.72019699314; -1.146797966946], 1e-10, 'Cm_B'; ...
        actual.Cf_s, [-0.645454548369; 0; 0], 1e-10, 'Cf_s'; ...
        actual.rho, 1.726247250065e-10, 1e-20, 'rho'; ...
        actual.vinf, 7784.261748565626, 1e-9, 'vinf' ...
    };
    passed = run_checks(checks);
end

function passed = test_sentman_reference_case()
    sentmanResults = compute_results(nominal_inputs('sentman'));
    cookResults = compute_results(nominal_inputs('cook'));
    checks = { ...
        sentmanResults.Cf_w, [-1.203529093977; 0; 0], 1e-10, 'Sentman Cf_w'; ...
        sentmanResults.Cm_B, [0; 1.8052936858; -1.203529093977], 1e-10, 'Sentman Cm_B'; ...
        sentmanResults.Cf_s, cookResults.Cf_s, 1e-12, 'Sentman Cf_s' ...
    };
    passed = run_checks(checks);
    passed = verify_true(sentmanResults.Cf_w(1) < cookResults.Cf_w(1), ...
        'Sentman should produce a larger drag magnitude than Cook for the reference case.') & passed;
end

function passed = test_position_eci_matches_geodetic_inputs()
    inputs = nominal_inputs('cook');
    positionEci = geodetic_to_eci(inputs.location, inputs.time);
    directResults = compute_results(inputs);
    inputs.location = struct('positionECI', positionEci);
    eciResults = compute_results(inputs);
    checks = { ...
        eciResults.Cf_w, directResults.Cf_w, 1e-12, 'Cf_w from positionECI'; ...
        eciResults.Cm_B, directResults.Cm_B, 1e-12, 'Cm_B from positionECI'; ...
        eciResults.Cf_s, directResults.Cf_s, 1e-12, 'Cf_s from positionECI'; ...
        eciResults.location.latitude, directResults.location.latitude, 1e-10, 'Resolved latitude'; ...
        eciResults.location.longitude, directResults.location.longitude, 1e-10, 'Resolved longitude'; ...
        eciResults.location.altitude, directResults.location.altitude, 1e-6, 'Resolved altitude' ...
    };
    passed = run_checks(checks);
end

function passed = test_unsupported_gsi_model_errors()
    try
        compute_results(nominal_inputs('cll'));
        passed = verify_true(false, 'Expected compute_aero_forces to reject unsupported gsi_model ''cll''.');
    catch ME
        passed = verify_true(strcmp(ME.identifier, 'computeAeroForces:UnsupportedGsiModel'), ...
            sprintf('Unexpected error identifier: %s', ME.identifier));
    end
end

function inputs = nominal_inputs(gsiModel)
    inputs = struct();
    inputs.objFile = vleo.util.geometry_asset_path('6U CubeSat.obj');
    inputs.location = struct('altitude', 200e3, 'latitude', 25.8, 'longitude', 0);
    inputs.attitude = struct('aoa', 0, 'aos', 0);
    inputs.time = struct('year', 2002, 'dayOfYear', 106, 'UTseconds', 0);
    inputs.options = struct( ...
        'f107Average', 65, ...
        'f107Daily', 65, ...
        'magneticIndex', ones(1, 7) * 3, ...
        'anomalousOxygen', 0, ...
        'gsi_model', gsiModel, ...
        'alpha', 1, ...
        'Tw', 300, ...
        'enableShadow', 1, ...
        'enableSolar', 1, ...
        'sol_cR', 0.15, ...
        'sol_cD', 0.25 ...
    );
end

function results = compute_results(inputs)
    results = vleo.aero.compute_aero_forces(inputs.objFile, inputs.location, inputs.attitude, inputs.time, inputs.options);
end

function positionEci = geodetic_to_eci(location, time)
    c = vleo.util.constants();
    latitude = deg2rad(location.latitude);
    longitude = deg2rad(location.longitude);
    primeVertical = c.R_earth / sqrt(1 - c.eccentricity_squared_earth * sin(latitude)^2);
    positionEcef = [ ...
        (primeVertical + location.altitude) * cos(latitude) * cos(longitude); ...
        (primeVertical + location.altitude) * cos(latitude) * sin(longitude); ...
        (primeVertical * (1 - c.eccentricity_squared_earth) + location.altitude) * sin(latitude) ...
    ];
    thetaGmst = greenwich_sidereal_time(julian_date_from_year_and_day(time.year, time.dayOfYear, time.UTseconds));
    rotation = [cos(thetaGmst), -sin(thetaGmst), 0; ...
                sin(thetaGmst),  cos(thetaGmst), 0; ...
                0,               0,              1];
    positionEci = rotation * positionEcef;
end

function jd = julian_date_from_year_and_day(year, dayOfYear, utSeconds)
    [month, day] = day_of_year_to_month_day(year, dayOfYear);
    hour = floor(utSeconds / 3600);
    minute = floor((utSeconds - 3600 * hour) / 60);
    second = utSeconds - 3600 * hour - 60 * minute;
    if month <= 2
        year = year - 1;
        month = month + 12;
    end
    century = floor(year / 100);
    correction = 2 - century + floor(century / 4);
    jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + ...
        day + correction - 1524.5 + (hour + (minute + second / 60) / 60) / 24;
end

function theta = greenwich_sidereal_time(julianDateUt1)
    tut1 = (julianDateUt1 - 2451545.0) / 36525.0;
    theta = -6.2e-6 * tut1^3 + 0.093104 * tut1^2 + ...
        (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
    theta = rem(theta * (pi / 180) / 240.0, 2 * pi);
    if theta < 0
        theta = theta + 2 * pi;
    end
end

function [month, day] = day_of_year_to_month_day(year, dayOfYear)
    monthLengths = [31, 28 + (mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)), ...
        31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    month = 1;
    day = dayOfYear;
    while day > monthLengths(month)
        day = day - monthLengths(month);
        month = month + 1;
    end
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
    if isscalar(expected) && isscalar(actual)
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
