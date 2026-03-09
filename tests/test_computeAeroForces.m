function results = test_computeAeroForces()
% TEST_COMPUTEAEROFORCES Regression and interface tests for computeAeroForces

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    results = struct('passed', 0, 'failed', 0, 'tests', {{}});

    results = runTest(results, @test_cook_reference_case);
    results = runTest(results, @test_sentman_reference_case);
    results = runTest(results, @test_position_eci_matches_geodetic_inputs);
    results = runTest(results, @test_unsupported_gsi_model_errors);

    printSummary(results);
end

function passed = test_cook_reference_case()
    [objFile, location, attitude, time, options] = nominalInputs('cook');
    results = computeAeroForces(objFile, location, attitude, time, options);

    passed = true;
    passed = passed && verifyClose(results.AreaRef, 0.110000005513, 1e-12, 'AreaRef');
    passed = passed && verifyClose(results.AreaProj, 0.060000003278, 1e-12, 'AreaProj');
    passed = passed && verifyClose(results.LenRef, 0.100000001490, 1e-12, 'LenRef');
    passed = passed && verifyVector(results.Cf_w, [-1.146797964779; 0; 0], 1e-10, 'Cf_w');
    passed = passed && verifyVector(results.Cm_B, [0; 1.720196989891; -1.146797964779], 1e-10, 'Cm_B');
    passed = passed && verifyVector(results.Cf_s, [-0.645454548369; 0; 0], 1e-10, 'Cf_s');
    passed = passed && verifyClose(results.rho, 1.726247250065e-10, 1e-20, 'rho');
    passed = passed && verifyClose(results.vinf, 7784.262050320798, 1e-9, 'vinf');
end

function passed = test_sentman_reference_case()
    [objFile, location, attitude, time, options] = nominalInputs('sentman');
    sentmanResults = computeAeroForces(objFile, location, attitude, time, options);

    options.gsi_model = 'cook';
    cookResults = computeAeroForces(objFile, location, attitude, time, options);

    passed = true;
    passed = passed && verifyVector(sentmanResults.Cf_w, [-1.203529089431; 0; 0], 1e-10, 'Sentman Cf_w');
    passed = passed && verifyVector(sentmanResults.Cm_B, [0; 1.805293678982; -1.203529089431], 1e-10, 'Sentman Cm_B');
    passed = passed && verifyVector(sentmanResults.Cf_s, cookResults.Cf_s, 1e-12, 'Sentman Cf_s');
    passed = passed && sentmanResults.Cf_w(1) < cookResults.Cf_w(1);

    if ~(sentmanResults.Cf_w(1) < cookResults.Cf_w(1))
        fprintf('    Sentman should produce a larger drag magnitude than Cook for the reference case.\n');
    end
end

function passed = test_position_eci_matches_geodetic_inputs()
    [objFile, location, attitude, time, options] = nominalInputs('cook');
    [positionEci] = geodeticToEci(location.latitude, location.longitude, location.altitude, ...
        time.year, time.dayOfYear, time.UTseconds);

    directResults = computeAeroForces(objFile, location, attitude, time, options);
    eciResults = computeAeroForces(objFile, struct('positionECI', positionEci), attitude, time, options);

    passed = true;
    passed = passed && verifyVector(eciResults.Cf_w, directResults.Cf_w, 1e-12, 'Cf_w from positionECI');
    passed = passed && verifyVector(eciResults.Cm_B, directResults.Cm_B, 1e-12, 'Cm_B from positionECI');
    passed = passed && verifyVector(eciResults.Cf_s, directResults.Cf_s, 1e-12, 'Cf_s from positionECI');
    passed = passed && verifyClose(eciResults.location.latitude, location.latitude, 1e-10, 'Resolved latitude');
    passed = passed && verifyClose(eciResults.location.longitude, location.longitude, 1e-10, 'Resolved longitude');
    passed = passed && verifyClose(eciResults.location.altitude, location.altitude, 1e-6, 'Resolved altitude');
end

function passed = test_unsupported_gsi_model_errors()
    [objFile, location, attitude, time, options] = nominalInputs('cll');

    try
        computeAeroForces(objFile, location, attitude, time, options);
        passed = false;
        fprintf('    Expected computeAeroForces to reject unsupported gsi_model ''cll''.\n');
    catch ME
        passed = strcmp(ME.identifier, 'computeAeroForces:UnsupportedGsiModel');
        if ~passed
            fprintf('    Unexpected error identifier: %s\n', ME.identifier);
        end
    end
end

function [objFile, location, attitude, time, options] = nominalInputs(gsiModel)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    objFile = fullfile(projectRoot, 'src', 'dynamics', '6U CubeSat.obj');

    location = struct('altitude', 200e3, 'latitude', 25.8, 'longitude', 0);
    attitude = struct('aoa', 0, 'aos', 0);
    time = struct('year', 2002, 'dayOfYear', 106, 'UTseconds', 0);
    options = struct( ...
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

function [positionEci] = geodeticToEci(latitudeDeg, longitudeDeg, altitude, year, dayOfYear, utSeconds)
    equatorialRadius = 6378137.0;
    flattening = 1 / 298.257223563;
    eccentricitySquared = flattening * (2 - flattening);

    latitude = deg2rad(latitudeDeg);
    longitude = deg2rad(longitudeDeg);
    primeVertical = equatorialRadius / sqrt(1 - eccentricitySquared * sin(latitude)^2);

    positionEcef = [ ...
        (primeVertical + altitude) * cos(latitude) * cos(longitude); ...
        (primeVertical + altitude) * cos(latitude) * sin(longitude); ...
        (primeVertical * (1 - eccentricitySquared) + altitude) * sin(latitude) ...
    ];

    thetaGmst = greenwichSiderealTime(julianDateFromYearAndDay(year, dayOfYear, utSeconds));
    rotation = [cos(thetaGmst), -sin(thetaGmst), 0; ...
                sin(thetaGmst),  cos(thetaGmst), 0; ...
                0,               0,              1];
    positionEci = rotation * positionEcef;
end

function jd = julianDateFromYearAndDay(year, dayOfYear, utSeconds)
    [month, day] = dayOfYearToMonthDay(year, dayOfYear);
    hour = floor(utSeconds / 3600);
    minute = floor((utSeconds - 3600 * hour) / 60);
    second = utSeconds - 3600 * hour - 60 * minute;

    if month <= 2
        year = year - 1;
        month = month + 12;
    end

    a = floor(year / 100);
    b = 2 - a + floor(a / 4);
    jd = floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) + ...
        day + b - 1524.5 + (hour + (minute + second / 60) / 60) / 24;
end

function theta = greenwichSiderealTime(julianDateUt1)
    tut1 = (julianDateUt1 - 2451545.0) / 36525.0;
    theta = -6.2e-6 * tut1^3 + 0.093104 * tut1^2 + ...
        (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
    theta = rem(theta * (pi / 180) / 240.0, 2 * pi);
    if theta < 0
        theta = theta + 2 * pi;
    end
end

function [month, day] = dayOfYearToMonthDay(year, dayOfYear)
    monthLengths = [31, 28 + isLeapYear(year), 31, 30, 31, 30, ...
                    31, 31, 30, 31, 30, 31];
    month = 1;
    day = dayOfYear;

    while day > monthLengths(month)
        day = day - monthLengths(month);
        month = month + 1;
    end
end

function tf = isLeapYear(year)
    tf = mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0);
end

function passed = verifyClose(actual, expected, tolerance, label)
    passed = abs(actual - expected) <= tolerance;
    if ~passed
        fprintf('    %s mismatch. Expected %.16g, got %.16g\n', label, expected, actual);
    end
end

function passed = verifyVector(actual, expected, tolerance, label)
    passed = all(abs(actual(:) - expected(:)) <= tolerance);
    if ~passed
        fprintf('    %s mismatch.\n', label);
        fprintf('    Expected: [%s]\n', num2str(expected(:).', ' %.16g'));
        fprintf('    Actual:   [%s]\n', num2str(actual(:).', ' %.16g'));
    end
end

function results = runTest(results, testFunc)
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

function printSummary(results)
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
