function results = test_fixed_axis_nadir_drag_case()
% TEST_FIXED_AXIS_NADIR_DRAG_CASE Compare drag for three fixed nadir-pointing axis assignments.

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(projectRoot);
    setup_project();

    clc;
    close all;

    [tSec, orbitRadius_m, meanMotion_rad_sec, time0, objFile, options, caseSpecs] = configure_test_case();
    caseResults(numel(caseSpecs), 1) = case_result_template();

    print_test_header(tSec, meanMotion_rad_sec, orbitRadius_m, options);

    for caseIdx = 1:numel(caseSpecs)
        caseResults(caseIdx) = evaluate_case(caseSpecs(caseIdx), tSec, orbitRadius_m, meanMotion_rad_sec, ...
            time0, objFile, options);
        plot_case(caseResults(caseIdx));
        print_case_summary(caseResults(caseIdx));
    end

    plot_combined_drag(caseResults);

    hypothesis = evaluate_hypothesis(caseResults);
    print_overall_summary(caseResults, hypothesis);

    results = struct( ...
        'passed', 1, ...
        'failed', 0, ...
        'tests', {{struct('name', 'test_fixed_axis_nadir_drag_case_runs', 'status', 'PASS')}}, ...
        'caseResults', caseResults, ...
        'hypothesis', hypothesis ...
    );
end

function [tSec, orbitRadius_m, meanMotion_rad_sec, time0, objFile, options, caseSpecs] = configure_test_case()
    c = vleo.util.constants();
    altitude_m = 200e3;
    orbitRadius_m = c.R_earth + altitude_m;
    meanMotion_rad_sec = sqrt(c.mu_earth / orbitRadius_m^3);
    orbitPeriod_sec = 2 * pi / meanMotion_rad_sec;

    tSec = linspace(0, orbitPeriod_sec, 121).';
    time0 = struct('year', 2002, 'dayOfYear', 106, 'UTseconds', 43200);
    objFile = vleo.util.geometry_asset_path('6U CubeSat.obj');
    options = struct( ...
        'f107Average', 65, ...
        'f107Daily', 65, ...
        'magneticIndex', ones(1, 7) * 3, ...
        'anomalousOxygen', 0, ...
        'gsi_model', 'cook', ...
        'alpha', 1, ...
        'Tw', 300, ...
        'enableShadow', 1, ...
        'enableSolar', 0, ...
        'sol_cR', 0.15, ...
        'sol_cD', 0.25 ...
    );
    caseSpecs = struct( ...
        'name', {'Case 1 - x y z', 'Case 2 - y z x', 'Case 3 - z x y'}, ...
        'northAxis', {'x', 'y', 'z'}, ...
        'velocityAxis', {'y', 'z', 'x'}, ...
        'nadirAxis', {'z', 'x', 'y'} ...
    );
end

function template = case_result_template()
    template = struct( ...
        'name', '', ...
        'northAxis', '', ...
        'velocityAxis', '', ...
        'nadirAxis', '', ...
        'attitude', struct('aoa', NaN, 'aos', NaN), ...
        'tSec', zeros(0, 1), ...
        'forceBody_N', zeros(0, 3), ...
        'coeffBody', zeros(0, 3), ...
        'dragForce_N', zeros(0, 1), ...
        'cdReference', zeros(0, 1), ...
        'cdProjected', zeros(0, 1), ...
        'density_kg_m3', zeros(0, 1), ...
        'projectedArea_m2', zeros(0, 1), ...
        'meanDragForce_N', NaN, ...
        'meanCdReference', NaN, ...
        'meanCdProjected', NaN, ...
        'meanForceBody_N', zeros(1, 3), ...
        'meanProjectedArea_m2', NaN ...
    );
end

function print_test_header(tSec, meanMotion_rad_sec, orbitRadius_m, options)
    c = vleo.util.constants();
    orbitPeriod_sec = tSec(end);
    altitude_km = (orbitRadius_m - c.R_earth) / 1000;

    fprintf('=== FIXED-AXIS NADIR DRAG TEST CASE ===\n');
    fprintf('Orbit: circular equatorial, altitude = %.0f km\n', altitude_km);
    fprintf('Orbital period: %.2f minutes\n', orbitPeriod_sec / 60);
    fprintf('Mean motion: %.6e rad/s\n', meanMotion_rad_sec);
    fprintf('Samples: %d over one orbit\n', numel(tSec));
    fprintf('Aero model: %s, solar = off\n\n', options.gsi_model);
end

function caseResult = evaluate_case(caseSpec, tSec, orbitRadius_m, meanMotion_rad_sec, time0, objFile, options)
    sampleCount = numel(tSec);
    attitude = aero_angles_from_velocity_body(axis_unit_vector(caseSpec.velocityAxis));

    forceBody_N = zeros(sampleCount, 3);
    coeffBody = zeros(sampleCount, 3);
    dragForce_N = zeros(sampleCount, 1);
    cdReference = zeros(sampleCount, 1);
    cdProjected = zeros(sampleCount, 1);
    density_kg_m3 = zeros(sampleCount, 1);
    projectedArea_m2 = zeros(sampleCount, 1);

    for idx = 1:sampleCount
        trueLongitude_rad = meanMotion_rad_sec * tSec(idx);
        positionEci_m = orbitRadius_m * [cos(trueLongitude_rad); sin(trueLongitude_rad); 0];

        time = time0;
        time.UTseconds = time0.UTseconds + tSec(idx);

        aeroResults = vleo.aero.compute_aero_forces(objFile, ...
            struct('positionECI', positionEci_m), attitude, time, options);

        dynamicPressure_Pa = 0.5 * aeroResults.rho * aeroResults.vinf^2;
        forceBody_N(idx, :) = (dynamicPressure_Pa * aeroResults.AreaRef * aeroResults.Cf_b).';
        coeffBody(idx, :) = aeroResults.Cf_b.';
        dragForce_N(idx) = -aeroResults.F_aero(1);
        cdReference(idx) = -aeroResults.Cf_w(1);
        cdProjected(idx) = dragForce_N(idx) / (dynamicPressure_Pa * aeroResults.AreaProj);
        density_kg_m3(idx) = aeroResults.rho;
        projectedArea_m2(idx) = aeroResults.AreaProj;
    end

    caseResult = case_result_template();
    caseResult.name = caseSpec.name;
    caseResult.northAxis = caseSpec.northAxis;
    caseResult.velocityAxis = caseSpec.velocityAxis;
    caseResult.nadirAxis = caseSpec.nadirAxis;
    caseResult.attitude = attitude;
    caseResult.tSec = tSec;
    caseResult.forceBody_N = forceBody_N;
    caseResult.coeffBody = coeffBody;
    caseResult.dragForce_N = dragForce_N;
    caseResult.cdReference = cdReference;
    caseResult.cdProjected = cdProjected;
    caseResult.density_kg_m3 = density_kg_m3;
    caseResult.projectedArea_m2 = projectedArea_m2;
    caseResult.meanDragForce_N = mean(dragForce_N);
    caseResult.meanCdReference = mean(cdReference);
    caseResult.meanCdProjected = mean(cdProjected);
    caseResult.meanForceBody_N = mean(forceBody_N, 1);
    caseResult.meanProjectedArea_m2 = mean(projectedArea_m2);
end

function vector = axis_unit_vector(axisName)
    vector = zeros(3, 1);
    vector(axis_index(axisName)) = 1;
end

function idx = axis_index(axisName)
    idx = find('xyz' == lower(axisName), 1, 'first');
    if isempty(idx) || numel(axisName) ~= 1
        error('test_fixed_axis_nadir_drag_case:InvalidAxis', ...
            'Axis name must be x, y, or z.');
    end
end

function attitude = aero_angles_from_velocity_body(velocityBody)
    velocityBody = velocityBody / norm(velocityBody);
    attitude = struct( ...
        'aoa', atan2d(velocityBody(3), velocityBody(1)), ...
        'aos', asind(velocityBody(2)) ...
    );
end

function plot_case(caseResult)
    [axisLabels, axisColors] = axis_plot_metadata();
    timeMinutes = caseResult.tSec / 60;

    figure('Color', 'w', 'Name', caseResult.name, 'Visible', default_figure_visibility());
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for axisIdx = 1:3
        nexttile(axisIdx);
        plot(timeMinutes, caseResult.forceBody_N(:, axisIdx), 'LineWidth', 1.6, 'Color', axisColors(axisIdx, :));
        grid on;
        xlabel('Time [min]');
        ylabel(sprintf('F_{b,%s} [N]', axisLabels{axisIdx}));
        title(sprintf('Body %s force', upper(axisLabels{axisIdx})));

        nexttile(axisIdx + 3);
        plot(timeMinutes, caseResult.coeffBody(:, axisIdx), 'LineWidth', 1.6, 'Color', axisColors(axisIdx, :));
        grid on;
        xlabel('Time [min]');
        ylabel(sprintf('C_{F,b,%s} [-]', axisLabels{axisIdx}));
        title(sprintf('Body %s coefficient', upper(axisLabels{axisIdx})));
    end

    sgtitle(sprintf(['%s | %s=NCP, %s=prograde, %s=nadir\n', ...
        'Mean drag = %.3e N | Mean C_D(A_{ref}) = %.3f | Mean C_D(A_{proj}) = %.3f'], ...
        caseResult.name, upper(caseResult.northAxis), upper(caseResult.velocityAxis), upper(caseResult.nadirAxis), ...
        caseResult.meanDragForce_N, caseResult.meanCdReference, caseResult.meanCdProjected));
end

function plot_combined_drag(caseResults)
    [~, axisColors] = axis_plot_metadata();
    progradeAxes = {'x', 'y', 'z'};

    figure('Color', 'w', 'Name', 'Combined prograde-axis drag', 'Visible', default_figure_visibility());
    ax = axes();
    hold(ax, 'on');
    grid(ax, 'on');

    for axisIdx = 1:numel(progradeAxes)
        caseIdx = find(strcmp({caseResults.velocityAxis}, progradeAxes{axisIdx}), 1, 'first');
        dragAlongAxis_N = abs(caseResults(caseIdx).forceBody_N(:, axis_index(progradeAxes{axisIdx})));
        dragAlongAxis_N = dragAlongAxis_N / max(dragAlongAxis_N);

        plot(ax, caseResults(caseIdx).tSec / 60, dragAlongAxis_N, ...
            'LineWidth', 1.8, 'Color', axisColors(axisIdx, :), ...
            'DisplayName', sprintf('%s prograde / max', upper(progradeAxes{axisIdx})));
    end

    xlabel(ax, 'Time [min]');
    ylabel(ax, 'Normalized drag [-]');
    title(ax, 'Normalized drag vs time with X, Y, and Z as the prograde axis');
    legend(ax, 'Location', 'best');
    ylim(ax, [0, 1.05]);
    hold(ax, 'off');
end

function [axisLabels, axisColors] = axis_plot_metadata()
    axisLabels = {'x', 'y', 'z'};
    axisColors = [0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880; 0, 0.4470, 0.7410];
end

function figureVisibility = default_figure_visibility()
    figureVisibility = 'on';
    if ~usejava('desktop')
        figureVisibility = 'off';
    end
end

function hypothesis = evaluate_hypothesis(caseResults)
    meanDragForce_N = [caseResults.meanDragForce_N];
    meanCdProjected = [caseResults.meanCdProjected];
    meanCdReference = [caseResults.meanCdReference];
    [~, dragOrderIdx] = sort(meanDragForce_N, 'descend');

    hypothesis = struct();
    hypothesis.dragOrderByNadirAxis = {caseResults(dragOrderIdx).nadirAxis};
    hypothesis.userExpectedDragOrder = {'z', 'y', 'x'};
    hypothesis.dragOrderMatchesUserExpectation = isequal(hypothesis.dragOrderByNadirAxis, hypothesis.userExpectedDragOrder);
    hypothesis.meanCdProjected = meanCdProjected;
    hypothesis.meanCdReference = meanCdReference;
    hypothesis.projectedCdRelativeSpread = (max(meanCdProjected) - min(meanCdProjected)) / mean(meanCdProjected);
    hypothesis.referenceCdRelativeSpread = (max(meanCdReference) - min(meanCdReference)) / mean(meanCdReference);
end

function print_case_summary(caseResult)
    fprintf('%s\n', caseResult.name);
    fprintf('  Axis mapping: %s=NCP, %s=prograde, %s=nadir\n', ...
        upper(caseResult.northAxis), upper(caseResult.velocityAxis), upper(caseResult.nadirAxis));
    fprintf('  Attitude inputs: aoa = %.1f deg, aos = %.1f deg\n', caseResult.attitude.aoa, caseResult.attitude.aos);
    fprintf('  Mean drag force:      %.3e N\n', caseResult.meanDragForce_N);
    fprintf('  Mean C_D (A_ref):     %.3f\n', caseResult.meanCdReference);
    fprintf('  Mean C_D (A_proj):    %.3f\n', caseResult.meanCdProjected);
    fprintf('  Mean body force [N]:  [%.3e, %.3e, %.3e]\n\n', caseResult.meanForceBody_N);
end

function print_overall_summary(caseResults, hypothesis)
    fprintf('Overall ranking by mean drag force: %s-down > %s-down > %s-down\n', ...
        upper(hypothesis.dragOrderByNadirAxis{1}), upper(hypothesis.dragOrderByNadirAxis{2}), upper(hypothesis.dragOrderByNadirAxis{3}));
    fprintf('Mean C_D (A_ref):    [%.3f, %.3f, %.3f]\n', hypothesis.meanCdReference);
    fprintf('Mean C_D (A_proj):   [%.3f, %.3f, %.3f]\n', hypothesis.meanCdProjected);
    fprintf('Projected-area C_D relative spread: %.2f %%\n', 100 * hypothesis.projectedCdRelativeSpread);
    fprintf('Fixed-reference C_D relative spread: %.2f %%\n', 100 * hypothesis.referenceCdRelativeSpread);

    if hypothesis.dragOrderMatchesUserExpectation
        fprintf('User drag-order hypothesis matches this case.\n');
    else
        fprintf('User drag-order hypothesis does not match this case.\n');
    end

    fprintf(['Note: C_D(A_ref) changes with attitude because the code uses a fixed reference area. ', ...
        'C_D(A_proj) is the better apples-to-apples comparison across orientations.\n']);

    for caseIdx = 1:numel(caseResults)
        fprintf('  %s projected area mean: %.4f m^2\n', caseResults(caseIdx).name, caseResults(caseIdx).meanProjectedArea_m2);
    end
end
