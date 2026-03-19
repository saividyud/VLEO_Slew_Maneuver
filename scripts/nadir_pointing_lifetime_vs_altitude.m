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
lowestDemonstratedOrbit_km = 167.4;
actuatorForceBudget_mN = 145;

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
    'Position', [120, 120, 1100, 760]);

lineColors = [0.10 0.35 0.75; 0.85 0.20 0.15; 0.10 0.60 0.25];
lineStyles = {'o-', 's-', 'd-'};

layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axThrust = nexttile(layout, 1);
plot(axThrust, altitudes_km, averageTotalThrust_mN, 'o-', ...
    'LineWidth', 1.9, 'MarkerSize', 6, ...
    'Color', [0.05 0.32 0.68], 'MarkerFaceColor', [0.05 0.32 0.68], ...
    'DisplayName', 'Average total thrust');
hold(axThrust, 'on');
plot(axThrust, altitudes_km, peakThrusterForce_mN, 's--', ...
    'LineWidth', 1.8, 'MarkerSize', 6, ...
    'Color', [0.85 0.25 0.12], 'MarkerFaceColor', [0.85 0.25 0.12], ...
    'DisplayName', 'Peak single-thruster force');
referenceLine = xline(axThrust, lowestDemonstratedOrbit_km, 'k-.', 'LineWidth', 1.2, ...
    'Label', 'SLATS record: 167.4 km', 'LabelOrientation', 'aligned', 'FontSize', 11);
vleo.util.hide_from_legend(referenceLine);
budgetLine = yline(axThrust, actuatorForceBudget_mN, ':', 'LineWidth', 1.4, 'Color', [0.35 0.35 0.35], ...
    'Label', '145 mN actuator budget', 'FontSize', 11);
vleo.util.hide_from_legend(budgetLine);
hold(axThrust, 'off');
grid(axThrust, 'on');
xlabel(axThrust, 'Altitude [km]');
ylabel(axThrust, 'Thrust [mN]');
title(axThrust, sprintf('%s nadir-pointing aerodynamic compensation', upper(thrusterLayoutMode)));
legend(axThrust, 'Location', 'northeast');
set(axThrust, 'FontSize', 13, 'LineWidth', 1.0);

axLife = nexttile(layout, 2);
hold(axLife, 'on');
lifeLabels = { ...
    'Electrospray-like, I_{sp} = 760 s', ...
    'Conventional ion, I_{sp} = 3000 s', ...
    'Advanced DS4G-like, I_{sp} = 19300 s'};
for ispIdx = 1:nIsp
    semilogy(axLife, altitudes_km, missionLife_days(ispIdx, :), lineStyles{ispIdx}, ...
        'LineWidth', 1.9, 'MarkerSize', 6, ...
        'Color', lineColors(ispIdx, :), ...
        'MarkerFaceColor', lineColors(ispIdx, :), ...
        'DisplayName', lifeLabels{ispIdx});
end
referenceLine = xline(axLife, lowestDemonstratedOrbit_km, 'k-.', 'LineWidth', 1.2, ...
    'Label', 'SLATS record altitude', 'LabelOrientation', 'aligned', 'FontSize', 11);
vleo.util.hide_from_legend(referenceLine);
line30 = yline(axLife, 30, ':', 'LineWidth', 1.2, 'Color', [0.35 0.35 0.35], 'Label', '30 days', 'FontSize', 11);
vleo.util.hide_from_legend(line30);
line365 = yline(axLife, 365, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5], 'Label', '1 year', 'FontSize', 11);
vleo.util.hide_from_legend(line365);
hold(axLife, 'off');
grid(axLife, 'on');
xlabel(axLife, 'Altitude [km]');
ylabel(axLife, 'Mission life [days]');
title(axLife, sprintf('Fuel-limited life for %.2f kg usable propellant (%s GSI model)', usableFuelMass_kg, gsi_model));
legend(axLife, 'Location', 'northwest');
set(axLife, 'FontSize', 13, 'LineWidth', 1.0, 'YMinorGrid', 'on', 'YScale', 'log');
ylim(axLife, [1, 3e3]);
yticks(axLife, [1, 3, 10, 30, 100, 300, 1000, 3000]);

title(layout, sprintf(['Nadir-pointing VLEO endurance vs altitude\n', ...
    '%d-orbit averaging, %d samples/orbit, %s layout'], ...
    nOrbits, samplesPerOrbit, thrusterLayoutMode), 'FontSize', 15);

figureExport = vleo.util.export_paper_figure(gcf, 'nadir_pointing_lifetime_vs_altitude');
fprintf('Paper figure saved to: %s\n', figureExport.pngPath);
fprintf('Paper figure saved to: %s\n', figureExport.pdfPath);

fprintf('\nDone.\n');
