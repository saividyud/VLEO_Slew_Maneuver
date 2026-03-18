% nadir_pointing_lifetime_vs_altitude.m
% Sweeps altitude and estimates propellant-limited mission life for
% continuous nadir-pointing aerodynamic compensation.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc; close all;

%% ---- Configuration ----------------------------------------------------
altitudes_km = 150:10:250;
nOrbits = 3;
samplesPerOrbit = 120;
gsi_model = 'cook';
thrusterLayoutMode = 'standard';
Isp_values_s = [760, 3000, 19300];
usableFuelMass_kg = 1.5;

nIsp = numel(Isp_values_s);
missionLife_days = zeros(nIsp, numel(altitudes_km));
averageTotalThrust_mN = zeros(size(altitudes_km));
peakThrusterForce_mN = zeros(size(altitudes_km));

fprintf('=== NADIR-POINTING LIFETIME VS ALTITUDE ===\n');
fprintf('Altitude range: %d:%d:%d km\n', altitudes_km(1), altitudes_km(2) - altitudes_km(1), altitudes_km(end));
fprintf('Duration per case: %d orbits\n', nOrbits);
fprintf('Thrusters: %s layout\n', thrusterLayoutMode);
fprintf('Isp values: %.0f s, %.0f s, %.0f s\n', Isp_values_s);
fprintf('Usable fuel: %.2f kg\n\n', usableFuelMass_kg);

g0 = 9.80665;

for idx = 1:numel(altitudes_km)
    survey = vleo.analysis.compute_nadir_pointing_aero_survey(altitudes_km(idx), ...
        'nOrbits', nOrbits, ...
        'samplesPerOrbit', samplesPerOrbit, ...
        'gsi_model', gsi_model, ...
        'thrusterLayoutMode', thrusterLayoutMode, ...
        'Isp_s', Isp_values_s(1), ...
        'usableFuelMass_kg', usableFuelMass_kg, ...
        'showProgress', false);

    averageTotalThrust_mN(idx) = survey.averageTotalThrusterForce_N * 1e3;
    peakThrusterForce_mN(idx) = max(survey.thrusterForces_N(:)) * 1e3;

    massFlowRate_kg_s = survey.averageTotalThrusterForce_N ./ (Isp_values_s(:) * g0);
    missionLife_days(:, idx) = usableFuelMass_kg ./ massFlowRate_kg_s / 86400;

    fprintf('[%02d/%02d] %3d km | avg total thrust = %7.3f mN | peak thruster = %6.3f mN | life = [%8.2f, %8.2f, %8.2f] days\n', ...
        idx, numel(altitudes_km), altitudes_km(idx), averageTotalThrust_mN(idx), ...
        peakThrusterForce_mN(idx), missionLife_days(1, idx), missionLife_days(2, idx), missionLife_days(3, idx));
end

figure('Color', 'w', 'Name', 'Nadir-Pointing Lifetime vs Altitude', ...
    'Position', [120, 120, 1100, 520]);

lineColors = [0.10 0.35 0.75; 0.85 0.20 0.15; 0.10 0.60 0.25];
lineStyles = {'o-', 's-', 'd-'};

hold on;
for ispIdx = 1:nIsp
    plot(altitudes_km, missionLife_days(ispIdx, :), lineStyles{ispIdx}, ...
        'LineWidth', 1.8, 'MarkerSize', 6, ...
        'Color', lineColors(ispIdx, :), ...
        'MarkerFaceColor', lineColors(ispIdx, :), ...
        'DisplayName', sprintf('Isp = %.0f s', Isp_values_s(ispIdx)));
end
hold off;
grid on;
xlabel('Altitude [km]');
ylabel('Mission life [days]');
title(sprintf(['Fuel-limited mission life | %s layout | %s model\n', ...
    '3-orbit averaging, Isp sweep = [%.0f, %.0f, %.0f] s, usable fuel = %.2f kg'], ...
    thrusterLayoutMode, gsi_model, Isp_values_s(1), Isp_values_s(2), Isp_values_s(3), usableFuelMass_kg));
legend('Location', 'northwest');

fprintf('\nDone.\n');
