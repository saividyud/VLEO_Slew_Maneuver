% nadir_pointing_aero_survey.m
% Normal-case survey: nadir-pointing satellite in equatorial orbit.
%   z-body (camera)  -> nadir
%   y-body           -> North Celestial Pole (NCP)
%   x-body           -> completes right-hand frame (anti-velocity)

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc; close all;

%% ---- Configuration ----------------------------------------------------
altitude_km = 200;
nOrbits = 3;
samplesPerOrbit = 120;
gsi_model = 'cook';
thrusterLayoutMode = 'standard';
Isp_s = 760;
usableFuelMass_kg = 1.5;

%% ---- Run survey -------------------------------------------------------
fprintf('=== NADIR-POINTING EQUATORIAL AERO SURVEY ===\n');
fprintf('Altitude:       %d km\n', altitude_km);
fprintf('Duration:       %d orbits\n', nOrbits);
fprintf('Samples/orbit:  %d\n', samplesPerOrbit);
fprintf('Pointing:       z-body = nadir, y-body = NCP\n');
fprintf('Body aoa/aos:   180 deg / 0 deg\n');
fprintf('GSI model:      %s\n', gsi_model);
fprintf('Thrusters:      %s layout\n', thrusterLayoutMode);
fprintf('Isp / Fuel:     %.0f s / %.2f kg\n\n', Isp_s, usableFuelMass_kg);

computeTimer = tic;
survey = vleo.analysis.compute_nadir_pointing_aero_survey(altitude_km, ...
    'nOrbits', nOrbits, ...
    'samplesPerOrbit', samplesPerOrbit, ...
    'gsi_model', gsi_model, ...
    'thrusterLayoutMode', thrusterLayoutMode, ...
    'Isp_s', Isp_s, ...
    'usableFuelMass_kg', usableFuelMass_kg, ...
    'showProgress', true);
fprintf('Computed %d points in %.1f s\n\n', numel(survey.tSec), toc(computeTimer));

%% ---- Summary statistics -----------------------------------------------
axNames = {'x', 'y', 'z'};
axColors = [0.8500 0.3250 0.0980; ...
            0.4660 0.6740 0.1880; ...
            0      0.4470 0.7410];

fprintf('Orbit period: %.1f min\n', survey.orbitalPeriod_s / 60);
fprintf('Density range: %.3e to %.3e kg/m^3\n', ...
    min(survey.density_kg_m3), max(survey.density_kg_m3));

fprintf('\nBody-frame force [mN]:\n');
for ax = 1:3
    fprintf('  %s: min = %+.4f   max = %+.4f   mean = %+.4f\n', ...
        axNames{ax}, min(survey.forceBody_N(:, ax)) * 1e3, ...
        max(survey.forceBody_N(:, ax)) * 1e3, mean(survey.forceBody_N(:, ax)) * 1e3);
end

fprintf('\nBody-frame torque [uNm]:\n');
for ax = 1:3
    fprintf('  %s: min = %+.4f   max = %+.4f   mean = %+.4f\n', ...
        axNames{ax}, min(survey.torqueBody_Nm(:, ax)) * 1e6, ...
        max(survey.torqueBody_Nm(:, ax)) * 1e6, mean(survey.torqueBody_Nm(:, ax)) * 1e6);
end

fprintf('\nThruster layout unit direction: [%.4f, %.4f, %.4f]\n', ...
    survey.layout.exhaust_direction_unit);
fprintf('Average total scalar thrust: %.4f mN\n', survey.averageTotalThrusterForce_N * 1e3);
fprintf('Peak individual thruster:    %.4f mN\n', max(survey.thrusterForces_N(:)) * 1e3);
fprintf('Mass flow rate:              %.4e kg/s  (%.2f g/day)\n', ...
    survey.massFlowRate_kg_s, survey.massFlowRate_g_day);
fprintf('Mission life:                %.2f days  (%.3f years)\n', ...
    survey.missionLife_days, survey.missionLife_years);

%% ---- Figure 1: aero force and torque ----------------------------------
figure('Color', 'w', 'Name', 'Nadir-Pointing Equatorial Aero Survey', ...
    'Position', [80, 80, 1400, 700]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for ax = 1:3
    nexttile(ax);
    plot(survey.timeMin, survey.forceBody_N(:, ax) * 1e3, ...
        'LineWidth', 1.4, 'Color', axColors(ax, :));
    grid on;
    xlabel('Time [min]');
    ylabel(sprintf('F_{%s} [mN]', axNames{ax}));
    title(sprintf('Aero force %s-body', upper(axNames{ax})));
    add_orbit_dividers(survey.orbitalPeriod_s, nOrbits);
end

for ax = 1:3
    nexttile(ax + 3);
    plot(survey.timeMin, survey.torqueBody_Nm(:, ax) * 1e6, ...
        'LineWidth', 1.4, 'Color', axColors(ax, :));
    grid on;
    xlabel('Time [min]');
    ylabel(sprintf('\\tau_{%s} [\\muNm]', axNames{ax}));
    title(sprintf('Aero torque %s-body', upper(axNames{ax})));
    add_orbit_dividers(survey.orbitalPeriod_s, nOrbits);
end

title(tl, sprintf(['Nadir-Pointing Equatorial Orbit | h = %d km | %s model\n', ...
    'z_{body} = nadir, y_{body} = NCP, x_{body} = anti-velocity'], ...
    altitude_km, gsi_model), 'FontSize', 13);

%% ---- Figure 2: required thruster forces --------------------------------
figure('Color', 'w', 'Name', 'Thruster Allocation', ...
    'Position', [140, 140, 1100, 850]);
tl2 = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for thrusterIdx = 1:8
    nexttile;
    plot(survey.timeMin, survey.thrusterForces_N(:, thrusterIdx) * 1e3, ...
        'LineWidth', 1.4, 'Color', [0.85 0.15 0.15]);
    grid on;
    xlabel('Time [min]');
    ylabel(sprintf('T%d [mN]', thrusterIdx));
    title(sprintf('Thruster %d', thrusterIdx));
    add_orbit_dividers(survey.orbitalPeriod_s, nOrbits);
end

title(tl2, sprintf(['Thruster Forces To Counter Aerodynamics | %s layout\n', ...
    'Avg total = %.3f mN, Mission life = %.2f days'], ...
    survey.layout.mode, survey.averageTotalThrusterForce_N * 1e3, survey.missionLife_days), ...
    'FontSize', 13);

fprintf('\nDone.\n');

%% ---- Local helpers ----------------------------------------------------
function add_orbit_dividers(orbitalPeriod_s, nOrbits)
    for orbitIdx = 1:nOrbits - 1
        xline(orbitIdx * orbitalPeriod_s / 60, 'k--', ...
            'Alpha', 0.3, 'HandleVisibility', 'off');
    end
end
