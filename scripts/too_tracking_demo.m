% Control Test Earth - Target of Opportunity (TOO) Observation Scenario
% Simulates a 6U CubeSat observing a TOO event with Earth rotation.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('includeAerodynamics', 'var')
    includeAerodynamics = true;
end

if ~exist('scenarioSeed', 'var')
    scenarioSeed = [];
end

if ~exist('generateVisuals', 'var')
    generateVisuals = usejava('desktop');
end

clearvars -except includeAerodynamics scenarioSeed generateVisuals
clc;

if ~(isscalar(includeAerodynamics) && (islogical(includeAerodynamics) || isnumeric(includeAerodynamics)))
    error('control_test3:InvalidAeroToggle', ...
        'includeAerodynamics must be a numeric or logical scalar.');
end
includeAerodynamics = logical(includeAerodynamics);

if isempty(scenarioSeed)
    rng('shuffle');
    scenarioSeedText = 'shuffle';
elseif isnumeric(scenarioSeed) && isscalar(scenarioSeed)
    scenarioSeedValue = double(scenarioSeed);
    if ~isfinite(scenarioSeedValue) || scenarioSeedValue < 0 || scenarioSeedValue ~= floor(scenarioSeedValue)
        error('control_test3:InvalidScenarioSeed', ...
            'scenarioSeed must be empty or a nonnegative integer scalar.');
    end
    rng(scenarioSeedValue, 'twister');
    scenarioSeedText = num2str(scenarioSeedValue);
else
    error('control_test3:InvalidScenarioSeed', ...
        'scenarioSeed must be empty or a nonnegative integer scalar.');
end

if ~(isscalar(generateVisuals) && (islogical(generateVisuals) || isnumeric(generateVisuals)))
    error('control_test3:InvalidGenerateVisuals', ...
        'generateVisuals must be a numeric or logical scalar.');
end
generateVisuals = logical(generateVisuals);

params = vleo.dynamics.default_control_test_params(includeAerodynamics);
scenario = vleo.analysis.generate_too_scenario(params);

fprintf('=== TOO OBSERVATION SCENARIO ===\n');
fprintf('Spacecraft: 6U CubeSat (10 x 20 x 30 cm, 12 kg, uniform density)\n');
fprintf('Aerodynamic disturbance torque: %s\n', vleo.util.on_off_text(params.includeAerodynamics));
fprintf('Scenario seed: %s\n', scenarioSeedText);
fprintf('Principal inertias [kg m^2]: [%.4f, %.4f, %.4f]\n\n', diag(params.I_CB));

fprintf('TOO event at t = 0\n');
fprintf('TOO Location (at t=0): Lat=%.2f deg, Lon=%.2f deg\n', scenario.too_lat, scenario.too_lon);
fprintf('Orbital Period: %.2f minutes\n', scenario.orbital_period / 60);
fprintf('Orbital Altitude: %.2f km\n', scenario.altitude / 1000);
fprintf('Max horizon angle (theta_max): %.2f degrees\n', rad2deg(scenario.theta_max));
fprintf('Selected theta: %.2f degrees, phi: %.2f degrees\n\n', rad2deg(scenario.theta), rad2deg(scenario.phi));
fprintf('Visibility check at t=0: %.2e (should be > 0)\n', scenario.visibility_at_eruption);
fprintf('Zenith distance at t=0: %.2f degrees (should be < 90)\n\n', scenario.zenith_distance_0);
fprintf('Simulation from t = %.2f to %.2f seconds\n', scenario.t_start, scenario.t_end);
fprintf('Total simulation time: %.2f minutes\n\n', (scenario.t_end - scenario.t_start) / 60);

fprintf('Propagating satellite trajectory...\n');
odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
odeFun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, scenario.z_body_eci, params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false);
[tspan, X_history] = vleo.dynamics.propagate_bidirectional(odeFun, scenario.X0, ...
    scenario.t_start, scenario.t_end, scenario.dt, odeOpts);

trackingHistory = vleo.analysis.compute_observation_history(tspan, X_history, scenario.r_too_0, params);
tau_total_history = trackingHistory.tau_track_total_history;

if params.includeAerodynamics
    fprintf('Evaluating aerodynamic torque history every %.0f seconds...\n', params.aeroSampleStep);
    tau_aero_history = vleo.dynamics.estimate_aero_torque_history(tspan, X_history, trackingHistory.R_B2E_hist, params);
else
    fprintf('Aerodynamic torque disabled; using zero disturbance torque.\n');
    tau_aero_history = zeros(size(tau_total_history));
end
tau_control_history = tau_total_history - tau_aero_history;

fprintf('Max |tau_total|   = %.3e N m\n', max(vecnorm(tau_total_history, 2, 2)));
fprintf('Max |tau_control| = %.3e N m\n', max(vecnorm(tau_control_history, 2, 2)));
fprintf('Max |tau_aero|    = %.3e N m\n\n', max(vecnorm(tau_aero_history, 2, 2)));

visibilityInfo = vleo.analysis.analyze_visibility(tspan, trackingHistory.visibility_history);
if ~isempty(visibilityInfo.idx_start_visible)
    t_start_visible = visibilityInfo.t_start_visible;
    fprintf('Satellite starts seeing TOO at t = %.2f s (%.2f min)\n', ...
        t_start_visible, t_start_visible / 60);
else
    t_start_visible = NaN;
    if trackingHistory.visibility_history(1)
        fprintf('TOO visible at simulation start\n');
    else
        fprintf('TOO not visible at simulation start\n');
    end
end

fprintf('TOO event starts at t = %.2f s (%.2f min)\n', scenario.t_eruption, scenario.t_eruption / 60);

if ~isempty(visibilityInfo.idx_end_visible)
    t_end_visible = visibilityInfo.t_end_visible;
    fprintf('TOO goes out of sight at t = %.2f s (%.2f min)\n\n', ...
        t_end_visible, t_end_visible / 60);
else
    t_end_visible = NaN;
    if trackingHistory.visibility_history(end)
        fprintf('TOO still visible at simulation end\n\n');
    else
        fprintf('TOO not visible at simulation end\n\n');
    end
end

simulations_dir = vleo.util.simulations_dir();
anim_data_file = fullfile(simulations_dir, sprintf('control_test3_anim_data_aero_%s_seed_%s.mat', ...
    vleo.util.on_off_text(params.includeAerodynamics), scenarioSeedText));
visibility_history = trackingHistory.visibility_history;
pointing_eci_history = trackingHistory.pointing_eci_history;
r_too_0 = scenario.r_too_0;
save(anim_data_file, 'tspan', 'X_history', 'visibility_history', ...
    'pointing_eci_history', 'r_too_0', 'params', 'tau_total_history', ...
    'tau_control_history', 'tau_aero_history');
fprintf('Animation data saved to: %s\n', anim_data_file);
fprintf('Simulation outputs directory: %s\n', simulations_dir);
fprintf('Run tools/render_gif.py to generate the GIF from the latest animation data.\n');

torque_csv_path = vleo.util.save_torque_history_csv('control_test3', tspan, visibility_history, ...
    tau_total_history, tau_control_history, tau_aero_history, params.includeAerodynamics, scenarioSeedText);
fprintf('Torque history CSV saved to: %s\n', torque_csv_path);

verification = vleo.analysis.run_tracking_verification(tspan, X_history, trackingHistory, visibilityInfo, ...
    scenario.r_too_0, tau_control_history, tau_aero_history, params, odeOpts);
plotTime = verification.plotTime;
plotDesiredEulerDeg = verification.plotDesiredEulerDeg;
plotActualEulerDeg = verification.plotActualEulerDeg;
ra_error = verification.ra_error;
dec_error = verification.dec_error;
tspan_verif = verification.tspan_verif;
verification_pointing_error_deg = verification.pointing_error_deg;

if ~isempty(ra_error)
    fprintf('Running Verification Simulation...\n');
    fprintf('Max RA Error: %.2e deg\n', max(abs(ra_error)));
    fprintf('Max Dec Error: %.2e deg\n', max(abs(dec_error)));
end

plotDesiredEulerDeg = vleo.util.unwrap_angle_history_deg(plotDesiredEulerDeg);
plotActualEulerDeg = vleo.util.unwrap_angle_history_deg(plotActualEulerDeg);

euler_csv_path = vleo.util.save_euler_history_csv('control_test3', plotTime, plotDesiredEulerDeg, ...
    plotActualEulerDeg, params.includeAerodynamics, scenarioSeedText);
fprintf('Euler history CSV saved to: %s\n', euler_csv_path);

animationIndices = 1:numel(tspan);
if ~isnan(t_start_visible) && ~isnan(t_end_visible)
    animationIndices = visibilityInfo.idx_start_visible:visibilityInfo.idx_end_visible;
end

animationTime = tspan(animationIndices);
animationStateHistory = X_history(animationIndices, :);
animationControlHistory = tau_control_history(animationIndices, :);
animationAeroHistory = tau_aero_history(animationIndices, :);

if ~isempty(tspan_verif)
    animationTime = tspan_verif;
    animationStateHistory = verification.X_verif;
end

guiResults = struct();
guiResults.t = animationTime(:);
guiResults.rs = animationStateHistory(:, 1:3);
guiResults.vs = animationStateHistory(:, 4:6);
guiResults.betas = animationStateHistory(:, 7:10);
guiResults.omegas = animationStateHistory(:, 11:13);
guiResults.torques = animationControlHistory;
guiResults.aeroTorques = animationAeroHistory;
guiResults.totalTorques = animationControlHistory + animationAeroHistory;
if isempty(tspan_verif)
    guiResults.eulerDeg = plotActualEulerDeg(animationIndices, :);
else
    guiResults.eulerDeg = plotActualEulerDeg;
end

if generateVisuals
    if usejava('desktop')
        fprintf('Generating GUI-compatible animated results plot...\n');
        vleo.viz.animate_results(struct('results', guiResults));
    else
        fprintf('Skipping GUI animation because MATLAB desktop is unavailable.\n');
    end

    visibleMask = trackingHistory.visibility_history;
    visibleTimesMin = tspan(visibleMask) / 60;
    actuatorForceHistory_mN = 1e3 * abs(tau_control_history) ./ params.momentArms';
    peakAxisForceHistory_mN = max(actuatorForceHistory_mN, [], 2);
    visibleActuatorForceHistory_mN = actuatorForceHistory_mN(visibleMask, :);
    visiblePeakAxisForceHistory_mN = peakAxisForceHistory_mN(visibleMask);
    peakForceVisible_mN = max(visiblePeakAxisForceHistory_mN);
    forceMarginFraction = peakForceVisible_mN / (145);

    fprintf('Generating paper-ready TOO tracking figure...\n');
    figure('Name', 'TOO Tracking Force and Error', 'Color', 'w', 'Position', [120 120 1180 760]);
    paperLayout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    axForce = nexttile(paperLayout, 1);
    plot(axForce, visibleTimesMin, visibleActuatorForceHistory_mN(:, 1), '-', 'LineWidth', 1.7, ...
        'Color', [0.84 0.22 0.16], 'DisplayName', 'Body-x');
    hold(axForce, 'on');
    plot(axForce, visibleTimesMin, visibleActuatorForceHistory_mN(:, 2), '-', 'LineWidth', 1.7, ...
        'Color', [0.15 0.42 0.78], 'DisplayName', 'Body-y');
    plot(axForce, visibleTimesMin, visibleActuatorForceHistory_mN(:, 3), '-', 'LineWidth', 1.7, ...
        'Color', [0.12 0.62 0.34], 'DisplayName', 'Body-z');
    plot(axForce, visibleTimesMin, visiblePeakAxisForceHistory_mN, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Peak axis demand');
    grid(axForce, 'on');
    xlabel(axForce, 'Time (minutes)', 'FontSize', 15);
    ylabel(axForce, 'Required actuator force (mN)', 'FontSize', 15);
    title(axForce, 'Representative TOO tracking force demand during the visible pass', 'FontSize', 16);
    legend(axForce, 'Location', 'northwest', 'FontSize', 12);
    mark_visibility_window(axForce, t_start_visible, scenario.t_eruption, t_end_visible, 11);
    vleo.util.set_visible_domain(axForce, t_start_visible, t_end_visible);
    text(axForce, 0.58, 0.92, sprintf('Peak = %.3f mN (%.3f%% of 145 mN budget)', ...
        peakForceVisible_mN, 100 * forceMarginFraction), ...
        'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 6);
    set(axForce, 'FontSize', 12, 'LineWidth', 1.0);

    axError = nexttile(paperLayout, 2);
    if ~isempty(verification_pointing_error_deg)
        plot(axError, tspan_verif / 60, verification_pointing_error_deg, '-', 'LineWidth', 1.9, ...
            'Color', [0.10 0.10 0.10], 'DisplayName', 'Angular pointing error');
        grid(axError, 'on');
        xlabel(axError, 'Time (minutes)', 'FontSize', 15);
        ylabel(axError, 'Pointing error (deg)', 'FontSize', 15);
        title(axError, 'Open-loop verification error using the commanded tracking torques', 'FontSize', 16);
        legend(axError, 'Location', 'northwest', 'FontSize', 12);
        xlim(axError, [tspan_verif(1) / 60, tspan_verif(end) / 60]);
        text(axError, 0.58, 0.90, sprintf('Max angular error = %.3f deg', ...
            max(verification_pointing_error_deg)), ...
            'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 6);
        set(axError, 'FontSize', 12, 'LineWidth', 1.0);
    else
        axis(axError, 'off');
        text(0.5, 0.5, 'Verification errors unavailable for this run.', ...
            'Parent', axError, 'HorizontalAlignment', 'center', 'FontSize', 13);
    end

    title(paperLayout, sprintf('TOO tracking at %.0f km, seed %s, aerodynamics %s', ...
        scenario.altitude / 1000, scenarioSeedText, upper(vleo.util.on_off_text(params.includeAerodynamics))), ...
        'FontSize', 17);

    figureStem = sprintf('too_tracking_force_and_error_seed_%s_aero_%s', ...
        scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
    figureExport = vleo.util.export_paper_figure(gcf, figureStem);
    fprintf('Paper figure saved to: %s\n', figureExport.pngPath);
    fprintf('Paper figure saved to: %s\n', figureExport.pdfPath);
else
    fprintf('Skipping GUI animation and figure generation (generateVisuals=false).\n');
end

fprintf('\nSimulation complete!\n');

function mark_visibility_window(ax, tStartVisible, tEruption, tEndVisible, fontSize)
    hold(ax, 'on');
    if ~isnan(tStartVisible)
        marker = xline(ax, tStartVisible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', fontSize);
        vleo.util.hide_from_legend(marker);
    end
    marker = xline(ax, tEruption / 60, 'k--', 'LineWidth', 2, 'Label', 'TOO Event (t=0)', 'FontSize', fontSize);
    vleo.util.hide_from_legend(marker);
    if ~isnan(tEndVisible)
        marker = xline(ax, tEndVisible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', fontSize);
        vleo.util.hide_from_legend(marker);
    end
end
