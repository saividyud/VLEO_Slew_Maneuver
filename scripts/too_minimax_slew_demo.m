% Control Test 4 - Nadir Hold, Minimax Slew, and Target of Opportunity (TOO) Tracking
% Holds nadir pointing before event, slews with a minimax torque profile,
% then matches the TOO-tracking attitude/rate so tracking can continue.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('includeAerodynamics', 'var')
    includeAerodynamics = true;
end

if ~exist('scenarioSeed', 'var')
    scenarioSeed = [];
end

if ~exist('targetMatchTime', 'var')
    targetMatchTime = [];
end

if ~exist('generateVisuals', 'var')
    generateVisuals = usejava('desktop');
end

clearvars -except includeAerodynamics scenarioSeed targetMatchTime generateVisuals
clc;

if ~(isscalar(includeAerodynamics) && (islogical(includeAerodynamics) || isnumeric(includeAerodynamics)))
    error('control_test4:InvalidAeroToggle', ...
        'includeAerodynamics must be a numeric or logical scalar.');
end
includeAerodynamics = logical(includeAerodynamics);

if isempty(scenarioSeed)
    rng('shuffle');
    scenarioSeedText = 'shuffle';
elseif isscalar(scenarioSeed) && isnumeric(scenarioSeed) && isfinite(scenarioSeed)
    scenarioSeedValue = double(scenarioSeed);
    if scenarioSeedValue >= 0 && scenarioSeedValue == floor(scenarioSeedValue)
        rng(scenarioSeedValue, 'twister');
        scenarioSeedText = num2str(scenarioSeedValue);
    else
        error('control_test4:InvalidScenarioSeed', ...
            'scenarioSeed must be empty or a nonnegative integer scalar.');
    end
else
    error('control_test4:InvalidScenarioSeed', ...
        'scenarioSeed must be empty or a nonnegative integer scalar.');
end

if ~(isempty(targetMatchTime) || (isscalar(targetMatchTime) && isnumeric(targetMatchTime) && isfinite(targetMatchTime)))
    error('control_test4:InvalidTargetMatchTime', ...
        'targetMatchTime must be empty or a finite numeric scalar.');
end

if ~(isscalar(generateVisuals) && (islogical(generateVisuals) || isnumeric(generateVisuals)))
    error('control_test4:InvalidGenerateVisuals', ...
        'generateVisuals must be a numeric or logical scalar.');
end
generateVisuals = logical(generateVisuals);

params = vleo.dynamics.default_control_test_params(includeAerodynamics);
scenario = vleo.analysis.generate_too_scenario(params);

fprintf('=== TOO OBSERVATION WITH MINIMAX ATTITUDE SLEW ===\n');
fprintf('Spacecraft: 6U CubeSat (10 x 20 x 30 cm, 12 kg, uniform density)\n');
fprintf('Aerodynamic disturbance torque estimate: %s\n', vleo.util.on_off_text(params.includeAerodynamics));
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

fprintf('Propagating satellite trajectory with nadir hold before eruption...\n');
odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
odeFun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, scenario.z_body_eci, params, ...
    'useJ2', false, 'useAtmDrag', false, 'useControl', false);
[tspan, X_orbit_history] = vleo.dynamics.propagate_bidirectional(odeFun, scenario.X0, ...
    scenario.t_start, scenario.t_end, scenario.dt, odeOpts);
n_steps = size(X_orbit_history, 1);
idx_eruption = find(abs(tspan - scenario.t_eruption) < 1e-12, 1, 'first');

trackingHistory = vleo.analysis.compute_observation_history(tspan, X_orbit_history, scenario.r_too_0, params);
ra_TOO_history = trackingHistory.ra_TOO_history;
dec_TOO_history = trackingHistory.dec_TOO_history;
visibility_history = trackingHistory.visibility_history;

visibilityInfo = vleo.analysis.analyze_visibility(tspan, visibility_history);
idx_start_visible = visibilityInfo.idx_start_visible;
idx_end_visible = visibilityInfo.idx_end_visible;

if ~isempty(idx_start_visible)
    t_start_visible = visibilityInfo.t_start_visible;
    fprintf('Satellite starts seeing TOO at t = %.2f s (%.2f min)\n', ...
        t_start_visible, t_start_visible / 60);
else
    t_start_visible = NaN;
    if visibility_history(1)
        fprintf('TOO visible at simulation start\n');
    else
        fprintf('TOO not visible at simulation start\n');
    end
end

fprintf('TOO event starts at t = %.2f s (%.2f min)\n', scenario.t_eruption, scenario.t_eruption / 60);

if ~isempty(idx_end_visible)
    t_end_visible = visibilityInfo.t_end_visible;
    fprintf('TOO goes out of sight at t = %.2f s (%.2f min)\n', ...
        t_end_visible, t_end_visible / 60);
else
    if visibility_history(end)
        error('control_test4:NoVisibilityEnd', ...
            ['TOO remains visible at simulation end, so the default match-time ', ...
             'rule is undefined in this window. Increase the simulation window.']);
    end
    error('control_test4:NoVisibilityEnd', ...
        'TOO is not visible after event in the current window.');
end

if t_end_visible <= scenario.t_eruption
    error('control_test4:InvalidVisibilityWindow', ...
        'The TOO is not visible long enough after event to define a post-event match time.');
end

t_match_default = scenario.t_eruption + 0.75 * (t_end_visible - scenario.t_eruption);
X_eruption = X_orbit_history(idx_eruption, :)';
q_initial = X_eruption(7:10) / norm(X_eruption(7:10));
omega_initial_body = X_eruption(11:13);
usedFallbackMatchScan = false;

if isempty(targetMatchTime)
    t_match_request = t_match_default;
    [~, idx_match_default] = min(abs(tspan - t_match_request));
    t_match_default_sampled = tspan(idx_match_default);
    if idx_match_default <= idx_eruption || idx_match_default >= idx_end_visible
        error('control_test4:SnappedMatchTimeOutOfRange', ...
            'The default match time snaps outside the valid post-eruption visible interval.');
    end

    q_target_default = trackingHistory.q_track_history(idx_match_default, :)';
    q_target_default = q_target_default / norm(q_target_default);
    if dot(q_initial, q_target_default) < 0
        q_target_default = -q_target_default;
    end
    omega_target_body_default = trackingHistory.omega_track_body_history(idx_match_default, :)';

    n_opt_intervals = 20;
    fprintf('Solving fast rotational slew with %d output intervals...\n', n_opt_intervals);
    fprintf('Using an analytic warm start with nonlinear direct-shooting refinement.\n');
    opt_result_default = vleo.control.solve_minimax_attitude_slew_fast_legacy(q_initial, omega_initial_body, q_target_default, ...
        omega_target_body_default, params.I_CB, t_match_default_sampled - scenario.t_eruption, n_opt_intervals, ...
        'refineNonlinear', true, 'nRefinementIntervals', [12, 16, 20], ...
        'momentArms', params.momentArms);

    if is_match_result_feasible(opt_result_default)
        idx_match = idx_match_default;
        t_match = t_match_default_sampled;
        q_target = q_target_default;
        omega_target_body = omega_target_body_default;
        opt_result = opt_result_default;
        matchSourceText = 'default 3/4 visibility rule';

        matchSelection = struct();
        matchSelection.candidateTimes = t_match;
        matchSelection.gammaHistory_N = opt_result.gammaForceEquivalent;
        matchSelection.angleErrorDegHistory = opt_result.terminalAngleErrorDeg;
        matchSelection.rateErrorDegPerSecHistory = opt_result.terminalRateErrorDegPerSec;
        matchSelection.isFeasibleHistory = true;
        matchSelection.selectedTime = t_match;
        matchSelection.defaultTime = t_match_request;
        matchSelection.defaultSampledTime = t_match;
        matchSelection.selectionText = matchSourceText;
    else
        usedFallbackMatchScan = true;
        fprintf(['Default handoff at %.2f s was not solver-feasible ', ...
            '(attitude error %.3e deg, rate error %.3e deg/s).\n'], ...
            t_match_default_sampled, opt_result_default.terminalAngleErrorDeg, ...
            opt_result_default.terminalRateErrorDegPerSec);
        fprintf('Scanning nearby candidate handoff times to recover a robust solution...\n');

        matchSelection = select_robust_match_time(tspan, idx_eruption, idx_end_visible, ...
            trackingHistory, q_initial, omega_initial_body, params.I_CB, params.momentArms, t_match_default);
        idx_match = matchSelection.idxSelected;
        t_match = matchSelection.selectedTime;
        q_target = matchSelection.qTarget;
        omega_target_body = matchSelection.omegaTargetBody;
        opt_result = matchSelection.resultSelected;
        matchSourceText = matchSelection.selectionText;
    end
else
    t_match_request = double(targetMatchTime);
    matchSourceText = 'user override';

    if t_match_request <= scenario.t_eruption || t_match_request >= t_end_visible
        error('control_test4:MatchTimeOutOfRange', ...
            'targetMatchTime must be strictly between eruption time and end-visible time.');
    end

    [~, idx_match] = min(abs(tspan - t_match_request));
    t_match = tspan(idx_match);
    if idx_match <= idx_eruption || idx_match >= idx_end_visible
        error('control_test4:SnappedMatchTimeOutOfRange', ...
            'The requested match time snaps outside the valid post-eruption visible interval.');
    end

    q_target = trackingHistory.q_track_history(idx_match, :)';
    q_target = q_target / norm(q_target);
    if dot(q_initial, q_target) < 0
        q_target = -q_target;
    end
    omega_target_body = trackingHistory.omega_track_body_history(idx_match, :)';

    n_opt_intervals = 20;
    fprintf('Solving fast rotational slew with %d output intervals...\n', n_opt_intervals);
    fprintf('Using an analytic warm start with nonlinear direct-shooting refinement.\n');
    opt_result = vleo.control.solve_minimax_attitude_slew_fast_legacy(q_initial, omega_initial_body, q_target, ...
        omega_target_body, params.I_CB, t_match - scenario.t_eruption, n_opt_intervals, ...
        'refineNonlinear', true, 'nRefinementIntervals', [12, 16, 20], ...
        'momentArms', params.momentArms);

    matchSelection = struct();
    matchSelection.candidateTimes = t_match;
    matchSelection.gammaHistory_N = opt_result.gammaForceEquivalent;
    matchSelection.angleErrorDegHistory = opt_result.terminalAngleErrorDeg;
    matchSelection.rateErrorDegPerSecHistory = opt_result.terminalRateErrorDegPerSec;
    matchSelection.isFeasibleHistory = is_match_result_feasible(opt_result);
    matchSelection.selectedTime = t_match;
    matchSelection.defaultTime = t_match_request;
    matchSelection.defaultSampledTime = t_match;
    matchSelection.selectionText = matchSourceText;
end

if usedFallbackMatchScan
    fprintf('Requested match time (default 3/4 visibility rule): %.2f s\n', t_match_request);
    fprintf('Auto-selection mode: %s\n', matchSourceText);
else
    fprintf('Requested match time (%s): %.2f s\n', matchSourceText, t_match_request);
end
fprintf('Using sampled match time: %.2f s (%.2f min)\n', t_match, t_match / 60);
fprintf('Maneuver duration after eruption: %.2f s\n\n', t_match - scenario.t_eruption);

matchScanFigureData = matchSelection;
if isempty(targetMatchTime) && isscalar(matchSelection.candidateTimes)
    fprintf('Evaluating nearby match-time scan for plotting...\n');
    matchScanFigureData = select_robust_match_time(tspan, idx_eruption, idx_end_visible, ...
        trackingHistory, q_initial, omega_initial_body, params.I_CB, params.momentArms, t_match_default);
end

fprintf('Target tracking state at t_match:\n');
fprintf('  Target camera RA/Dec: [%.3f, %.3f] deg\n', ...
    ra_TOO_history(idx_match), dec_TOO_history(idx_match));
fprintf('  Target body-rate magnitude: %.3e rad/s\n\n', norm(omega_target_body));

fprintf('Slew solver: %s\n', opt_result.solverUsed);
fprintf('Slew objective mode: %s\n', opt_result.objectiveMode);
fprintf('Optimization exit flag: %d\n', opt_result.exitflag);
fprintf('Optimization message: %s\n', opt_result.message);
fprintf('Peak component torque gamma: %.3e N m\n', opt_result.gammaTorqueEquivalent);
fprintf('Peak equivalent force gamma: %.3e N\n', opt_result.gammaForceEquivalent);
fprintf('Peak axis forces [N]: [%.3e, %.3e, %.3e]\n', opt_result.peakAxisForceAbs);
fprintf('Moment arms [m]: [%.2f, %.2f, %.2f]\n', params.momentArms);
fprintf('Grid max |tau_i|: %.3e N m\n', opt_result.maxTorqueNodeComponentAbs);
fprintf('Grid max ||tau||_2: %.3e N m\n', opt_result.maxTorqueNodeNorm);
fprintf('Max equality residual: %.3e\n', opt_result.maxEqualityResidual);
fprintf('Max inequality violation: %.3e\n', opt_result.maxInequalityViolation);
fprintf('Continuous terminal attitude error: %.3e deg\n', opt_result.terminalAngleErrorDeg);
fprintf('Continuous terminal body-rate error: %.3e deg/s\n\n', opt_result.terminalRateErrorDegPerSec);

tau_total_history = zeros(n_steps, 3);
if isfield(opt_result, 'interpMethod')
    maneuverInterpMethod = opt_result.interpMethod;
else
    maneuverInterpMethod = 'pchip';
end
tau_maneuver_history = interp1(opt_result.timeNodes, opt_result.tauNodes, ...
    tspan(idx_eruption:idx_match) - scenario.t_eruption, maneuverInterpMethod);
tau_total_history(idx_eruption:idx_match, :) = tau_maneuver_history;

if idx_match < idx_end_visible
    tau_total_history(idx_match + 1:idx_end_visible, :) = ...
        trackingHistory.tau_track_total_history(idx_match + 1:idx_end_visible, :);
end

t_control = tspan(idx_eruption:end);
torqueHistoryInterpMethod = 'previous';
maneuverTimeNodes = scenario.t_eruption + opt_result.timeNodes;
tau_total_interp = @(t) interpolate_maneuver_then_tracking_torque(t, maneuverTimeNodes, ...
    opt_result.tauNodes, t_match, tspan, trackingHistory.tau_track_total_history, torqueHistoryInterpMethod);
zero_torque_interp = @(t) zeros(3, 1);
ode_fun_total = @(t, X) vleo.dynamics.sat_dynamics_openloop(t, X, params, tau_total_interp, zero_torque_interp);
[~, X_post_total] = ode45(ode_fun_total, t_control, X_eruption, odeOpts);

X_total_history = X_orbit_history;
X_total_history(idx_eruption:end, :) = X_post_total;
X_total_history(:, 7:10) = vleo.util.normalize_quaternion_rows(X_total_history(:, 7:10));
X_total_history(:, 7:10) = vleo.util.align_quaternion_signs(X_total_history(:, 7:10));

X_desired_history = X_total_history;
X_desired_history(1:idx_eruption - 1, :) = X_orbit_history(1:idx_eruption - 1, :);
X_desired_history(idx_match:idx_end_visible, 7:10) = trackingHistory.q_track_history(idx_match:idx_end_visible, :);
X_desired_history(idx_match:idx_end_visible, 11:13) = trackingHistory.omega_track_body_history(idx_match:idx_end_visible, :);
X_desired_history(:, 7:10) = vleo.util.normalize_quaternion_rows(X_desired_history(:, 7:10));
X_desired_history(:, 7:10) = vleo.util.align_quaternion_signs(X_desired_history(:, 7:10));

R_desired_history = vleo.util.body_to_eci_history_from_quaternion_rows(X_desired_history(:, 7:10));
if params.includeAerodynamics
    fprintf('Evaluating aerodynamic torque history every %.0f seconds...\n', params.aeroSampleStep);
    tau_aero_history = vleo.dynamics.estimate_aero_torque_history(tspan, X_desired_history, R_desired_history, params);
else
    fprintf('Aerodynamic torque estimate disabled; using zero disturbance torque.\n');
    tau_aero_history = zeros(n_steps, 3);
end

tau_control_history = tau_total_history - tau_aero_history;

tau_aero_interp = @(t) interp1(tspan, tau_aero_history, t, 'pchip', 'extrap')';
tau_control_interp = @(t) tau_total_interp(t) - tau_aero_interp(t);
ode_fun_actual = @(t, X) vleo.dynamics.sat_dynamics_openloop(t, X, params, tau_control_interp, tau_aero_interp);
[~, X_post_actual] = ode45(ode_fun_actual, t_control, X_eruption, odeOpts);

X_actual_history = X_orbit_history;
X_actual_history(idx_eruption:end, :) = X_post_actual;
X_actual_history(:, 7:10) = vleo.util.normalize_quaternion_rows(X_actual_history(:, 7:10));
X_actual_history(:, 7:10) = vleo.util.align_quaternion_signs(X_actual_history(:, 7:10));

plotTime = tspan(:);
plotDesiredEulerDeg = vleo.util.quat_history_to_euler_deg(X_desired_history(:, 7:10));
plotActualEulerDeg = vleo.util.quat_history_to_euler_deg(X_actual_history(:, 7:10));
plotDesiredEulerDeg = vleo.util.unwrap_angle_history_deg(plotDesiredEulerDeg);
plotActualEulerDeg = vleo.util.unwrap_angle_history_deg(plotActualEulerDeg);

omega_body_desired_history = X_desired_history(:, 11:13);
omega_body_actual_history = X_actual_history(:, 11:13);

ra_camera_actual_history = zeros(n_steps, 1);
dec_camera_actual_history = zeros(n_steps, 1);
pointing_eci_actual_history = zeros(n_steps, 3);
pointing_error_deg_history = zeros(n_steps, 1);

for i = 1:n_steps
    obs_actual = vleo.analysis.state_to_observation(X_actual_history(i, :)', params);
    ra_camera_actual_history(i) = obs_actual.ra;
    dec_camera_actual_history(i) = obs_actual.dec;
    pointing_eci_actual_history(i, :) = obs_actual.pointing_eci';

    r_sat = X_actual_history(i, 1:3)';
    [rTOO, ~] = vleo.dynamics.too_state_at_time(tspan(i), scenario.r_too_0, params.omega_earth);
    apparentTOOVec = rTOO - r_sat;
    apparentTOOVec = apparentTOOVec / norm(apparentTOOVec);
    pointing_error_deg_history(i) = rad2deg(acos(vleo.util.clamp_scalar(dot(obs_actual.pointing_eci, apparentTOOVec), -1, 1)));
end

tracking_indices = idx_match:idx_end_visible;
ra_error = mod(ra_camera_actual_history(tracking_indices) - ra_TOO_history(tracking_indices) + 180, 360) - 180;
dec_error = dec_camera_actual_history(tracking_indices) - dec_TOO_history(tracking_indices);
body_rate_tracking_error_deg_s = vecnorm(rad2deg(omega_body_actual_history(tracking_indices, :) - ...
    trackingHistory.omega_track_body_history(tracking_indices, :)), 2, 2);

match_attitude_error_deg = vleo.util.rotation_error_angle_deg(X_actual_history(idx_match, 7:10)', q_target);
match_rate_error_deg_s = norm(rad2deg(omega_body_actual_history(idx_match, :)' - omega_target_body));

fprintf('Match-state verification at t = %.2f s:\n', t_match);
fprintf('  Attitude error: %.3e deg\n', match_attitude_error_deg);
fprintf('  Body-rate error: %.3e deg/s\n', match_rate_error_deg_s);
fprintf('Post-match max pointing error: %.3e deg\n', max(pointing_error_deg_history(tracking_indices)));
fprintf('Post-match max RA/Dec error: [%.3e, %.3e] deg\n', max(abs(ra_error)), max(abs(dec_error)));
fprintf('Post-match max body-rate tracking error: %.3e deg/s\n', max(body_rate_tracking_error_deg_s));
fprintf('Max ||tau_total||_2   = %.3e N m\n', max(vecnorm(tau_total_history, 2, 2)));
fprintf('Max ||tau_control||_2 = %.3e N m\n', max(vecnorm(tau_control_history, 2, 2)));
fprintf('Max ||tau_aero||_2    = %.3e N m\n', max(vecnorm(tau_aero_history, 2, 2)));
fprintf('Max required actuator force = %.3e N\n\n', max(max(abs(tau_control_history) ./ params.momentArms')));

simulations_dir = vleo.util.simulations_dir();
anim_data_file = fullfile(simulations_dir, sprintf('control_test4_anim_data_aero_%s_seed_%s.mat', ...
    vleo.util.on_off_text(params.includeAerodynamics), scenarioSeedText));
r_too_0 = scenario.r_too_0;
save(anim_data_file, 'tspan', 'X_actual_history', 'visibility_history', ...
    'pointing_eci_actual_history', 'r_too_0', 'params', 'tau_total_history', ...
    'tau_control_history', 'tau_aero_history', 't_match');
fprintf('Animation data saved to: %s\n', anim_data_file);
fprintf('Simulation outputs directory: %s\n', simulations_dir);
fprintf('Run tools/render_gif.py to generate the GIF from the latest animation data.\n');

torque_csv_path = vleo.util.save_torque_history_csv('control_test4', tspan, visibility_history, ...
    tau_total_history, tau_control_history, tau_aero_history, params.includeAerodynamics, scenarioSeedText);
fprintf('Torque history CSV saved to: %s\n', torque_csv_path);

attitude_csv_path = vleo.util.save_attitude_history_csv('control_test4', plotTime, plotDesiredEulerDeg, ...
    plotActualEulerDeg, omega_body_desired_history, omega_body_actual_history, ...
    params.includeAerodynamics, scenarioSeedText);
fprintf('Attitude history CSV saved to: %s\n', attitude_csv_path);

if generateVisuals
    animationIndices = 1:numel(tspan);
    if ~isnan(t_start_visible) && ~isnan(t_end_visible)
        animationIndices = idx_start_visible:idx_end_visible;
    end

    guiResults = struct();
    guiResults.t = tspan(animationIndices)';
    guiResults.rs = X_actual_history(animationIndices, 1:3);
    guiResults.vs = X_actual_history(animationIndices, 4:6);
    guiResults.betas = X_actual_history(animationIndices, 7:10);
    guiResults.omegas = X_actual_history(animationIndices, 11:13);
    guiResults.torques = tau_control_history(animationIndices, :);
    guiResults.aeroTorques = tau_aero_history(animationIndices, :);
    guiResults.totalTorques = tau_total_history(animationIndices, :);
    guiResults.eulerDeg = plotActualEulerDeg(animationIndices, :);

    if usejava('desktop')
        fprintf('Generating GUI-compatible animated results plot...\n');
        vleo.viz.animate_results(struct('results', guiResults));
    else
        fprintf('Skipping GUI animation because MATLAB desktop is unavailable.\n');
    end

    if ~isempty(idx_start_visible)
        visibleIndices = idx_start_visible:idx_end_visible;
    else
        visibleIndices = 1:idx_end_visible;
    end
    actuatorForceHistory_mN = 1e3 * abs(tau_control_history) ./ params.momentArms';
    peakAxisForceHistory_mN = max(actuatorForceHistory_mN, [], 2);
    peakForceVisible_mN = max(peakAxisForceHistory_mN(visibleIndices));
    forceMarginPercent = 100 * peakForceVisible_mN / 145;

    if usedFallbackMatchScan
        fprintf(['Default handoff needed a local fallback search; the exported scan now ', ...
            'documents the selected robust match time.\n']);
    end

    if ~isscalar(matchScanFigureData.candidateTimes)
        fprintf('Evaluating aerodynamic-aware match-time scan...\n');
        matchScanFigureData = evaluate_match_time_scan_with_aero(matchScanFigureData, ...
            tspan, idx_start_visible, idx_eruption, idx_end_visible, X_orbit_history, X_eruption, ...
            trackingHistory, params, odeOpts);

        fprintf('Generating handoff-time scan figure...\n');
        candidateTimesMin = matchScanFigureData.candidateTimes(:) / 60;
        candidatePeakForceWithAero_mN = 1e3 * matchScanFigureData.peakActuatorForceWithAero_N(:);
        candidatePeakForceManeuverOnly_mN = 1e3 * matchScanFigureData.gammaHistory_N(:);
        feasibleMask = matchScanFigureData.isFeasibleHistory(:);
        selectedMask = abs(matchScanFigureData.candidateTimes(:) - t_match) < 1e-12;
        defaultMask = abs(matchScanFigureData.candidateTimes(:) - matchScanFigureData.defaultSampledTime) < 1e-12;

        figure('Name', 'TOO Handoff Time Scan', 'Color', 'w', 'Position', [120 120 1080 620]);
        axScan = axes();
        plot(axScan, candidateTimesMin, candidatePeakForceWithAero_mN, 'k-', ...
            'LineWidth', 2.0, 'DisplayName', 'Peak force with aero');
        hold(axScan, 'on');
        plot(axScan, candidateTimesMin, candidatePeakForceManeuverOnly_mN, '--', ...
            'Color', [0.55 0.55 0.55], 'LineWidth', 1.3, 'DisplayName', 'Maneuver-only peak force');
        hold(axScan, 'on');
        scatter(axScan, candidateTimesMin(feasibleMask), candidatePeakForceWithAero_mN(feasibleMask), 40, ...
            [0.10 0.55 0.25], 'filled', 'DisplayName', 'Feasible handoff');
        scatter(axScan, candidateTimesMin(~feasibleMask), candidatePeakForceWithAero_mN(~feasibleMask), 55, ...
            [0.85 0.22 0.18], 'x', 'LineWidth', 1.6, 'DisplayName', 'Infeasible handoff');
        if any(defaultMask)
            scatter(axScan, candidateTimesMin(defaultMask), candidatePeakForceWithAero_mN(defaultMask), 85, ...
                [0.15 0.42 0.78], 'd', 'filled', 'DisplayName', 'Default match time');
        end
        scatter(axScan, candidateTimesMin(selectedMask), candidatePeakForceWithAero_mN(selectedMask), 110, ...
            'k', 'p', 'filled', 'DisplayName', 'Selected match time');
        grid(axScan, 'on');
        xlabel(axScan, 'Match time after event (minutes)', 'FontSize', 15);
        ylabel(axScan, 'Peak actuator force equivalent (mN)', 'FontSize', 15);
        title(axScan, 'Robust handoff-time scan including aerodynamic compensation', 'FontSize', 17);
        legend(axScan, 'Location', 'northwest', 'FontSize', 12);
        set(axScan, 'FontSize', 12, 'LineWidth', 1.0, 'YScale', 'log');
        text(axScan, 0.59, 0.88, sprintf([ ...
            'Selected match: %.2f min\n', ...
            'Peak with aero = %.3f mN\n', ...
            'Maneuver only = %.3f mN'], ...
            t_match / 60, matchScanFigureData.peakActuatorForceWithAeroSelected_N * 1e3, ...
            opt_result.gammaForceEquivalent * 1e3), ...
            'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 6);

        scanStem = sprintf('too_minimax_handoff_scan_seed_%s_aero_%s', ...
            scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
        scanExport = vleo.util.export_paper_figure(gcf, scanStem);
        fprintf('Paper figure saved to: %s\n', scanExport.pngPath);
        fprintf('Paper figure saved to: %s\n', scanExport.pdfPath);
    end

    fprintf('Generating attitude summary plots...\n');
    figure('Name', 'Euler Angle Tracking', 'WindowStyle', 'normal', 'Color', 'w', ...
        'Position', [120 120 1500 650]);
    axEuler = axes();
    l1 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 1), 'r-', 'LineWidth', 2); hold(axEuler, 'on');
    l2 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 1), 'r--', 'LineWidth', 2);
    l3 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 2), 'b-', 'LineWidth', 2);
    l4 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 2), 'b--', 'LineWidth', 2);
    l5 = plot(axEuler, plotTime / 60, plotDesiredEulerDeg(:, 3), 'k-', 'LineWidth', 2);
    l6 = plot(axEuler, plotTime / 60, plotActualEulerDeg(:, 3), 'k--', 'LineWidth', 2);
    grid(axEuler, 'on');
    xlabel(axEuler, 'Time (minutes)', 'FontSize', 16);
    ylabel(axEuler, 'Angle (deg)', 'FontSize', 16);
    title(axEuler, 'Desired and Actual Euler Angles', 'FontSize', 18);
    legend(axEuler, [l1, l2, l3, l4, l5, l6], ...
        {'Roll desired', 'Roll actual', 'Pitch desired', 'Pitch actual', 'Yaw desired', 'Yaw actual'}, ...
        'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
    mark_phase_lines(axEuler, t_start_visible, scenario.t_eruption, t_match, t_end_visible, 13);
    vleo.util.set_visible_domain(axEuler, t_start_visible, t_end_visible);
    set(axEuler, 'FontSize', 12, 'LineWidth', 1.0);
    eulerStem = sprintf('too_minimax_euler_tracking_seed_%s_aero_%s', ...
        scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
    eulerExport = vleo.util.export_paper_figure(gcf, eulerStem);
    fprintf('Paper figure saved to: %s\n', eulerExport.pngPath);
    fprintf('Paper figure saved to: %s\n', eulerExport.pdfPath);

    figure('Name', 'Body Rate Tracking', 'WindowStyle', 'normal', 'Color', 'w', ...
        'Position', [140 140 1500 650]);
    axRates = axes();
    r1 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 1)), 'r-', 'LineWidth', 2); hold(axRates, 'on');
    r2 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 1)), 'r--', 'LineWidth', 2);
    r3 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 2)), 'b-', 'LineWidth', 2);
    r4 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 2)), 'b--', 'LineWidth', 2);
    r5 = plot(axRates, plotTime / 60, rad2deg(omega_body_desired_history(:, 3)), 'k-', 'LineWidth', 2);
    r6 = plot(axRates, plotTime / 60, rad2deg(omega_body_actual_history(:, 3)), 'k--', 'LineWidth', 2);
    grid(axRates, 'on');
    xlabel(axRates, 'Time (minutes)', 'FontSize', 16);
    ylabel(axRates, 'Body Rate (deg/s)', 'FontSize', 16);
    title(axRates, 'Desired and Actual Body Rates', 'FontSize', 18);
    legend(axRates, [r1, r2, r3, r4, r5, r6], ...
        {'wx desired', 'wx actual', 'wy desired', 'wy actual', 'wz desired', 'wz actual'}, ...
        'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
    mark_phase_lines(axRates, t_start_visible, scenario.t_eruption, t_match, t_end_visible, 13);
    vleo.util.set_visible_domain(axRates, t_start_visible, t_end_visible);
    set(axRates, 'FontSize', 12, 'LineWidth', 1.0);
    ratesStem = sprintf('too_minimax_body_rate_tracking_seed_%s_aero_%s', ...
        scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
    ratesExport = vleo.util.export_paper_figure(gcf, ratesStem);
    fprintf('Paper figure saved to: %s\n', ratesExport.pngPath);
    fprintf('Paper figure saved to: %s\n', ratesExport.pdfPath);

    fprintf('Generating torque breakdown plots...\n');
    figure('Name', sprintf('Tracking Torque Breakdown (%s)', vleo.util.on_off_text(params.includeAerodynamics)), ...
        'WindowStyle', 'normal', 'Color', 'w', 'Position', [160 160 1450 900]);
    torqueLayout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    axisLabels = {'x', 'y', 'z'};

    for axisIdx = 1:3
        ax = nexttile(torqueLayout, axisIdx);
        plot(ax, tspan / 60, tau_total_history(:, axisIdx), 'k-', 'LineWidth', 2.2); hold(ax, 'on');
        plot(ax, tspan / 60, tau_control_history(:, axisIdx), 'b-', 'LineWidth', 1.8);
        plot(ax, tspan / 60, tau_aero_history(:, axisIdx), 'r--', 'LineWidth', 1.8);
        grid(ax, 'on');
        ylabel(ax, sprintf('\\tau_%s (N\\cdotm)', axisLabels{axisIdx}), 'FontSize', 16);
        title(ax, sprintf('Body-%s Torque Component', upper(axisLabels{axisIdx})), 'FontSize', 17);
        mark_phase_lines(ax, t_start_visible, scenario.t_eruption, t_match, t_end_visible, 13);
        vleo.util.set_visible_domain(ax, t_start_visible, t_end_visible);
        set(ax, 'FontSize', 12, 'LineWidth', 1.0);
        if axisIdx == 1
            legend(ax, {'Total required', 'Actuator required', 'Aerodynamic disturbance'}, ...
                'Location', 'best', 'FontSize', 13, 'AutoUpdate', 'off');
            text(ax, 0.68, 0.86, sprintf('Peak actuator force = %.3f mN\n(%.3f%% of 145 mN budget)', ...
                peakForceVisible_mN, forceMarginPercent), ...
                'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', 'Margin', 6);
        end
    end
    xlabel(torqueLayout, 'Time (minutes)', 'FontSize', 16);
    title(torqueLayout, sprintf('Torque Breakdown with Aerodynamics %s', upper(vleo.util.on_off_text(params.includeAerodynamics))), ...
        'FontSize', 18);
    torqueStem = sprintf('too_minimax_torque_breakdown_seed_%s_aero_%s', ...
        scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
    torqueExport = vleo.util.export_paper_figure(gcf, torqueStem);
    fprintf('Paper figure saved to: %s\n', torqueExport.pngPath);
    fprintf('Paper figure saved to: %s\n', torqueExport.pdfPath);

    figure('Name', 'Post-Match Tracking Error', 'WindowStyle', 'normal', 'Color', 'w', ...
        'Position', [180 180 1250 600]);
    plot(tspan(tracking_indices) / 60, pointing_error_deg_history(tracking_indices), 'k-', 'LineWidth', 2); hold on;
    plot(tspan(tracking_indices) / 60, body_rate_tracking_error_deg_s, 'm--', 'LineWidth', 2);
    xlabel('Time (minutes)', 'FontSize', 18);
    ylabel('Error', 'FontSize', 18);
    title('Post-Match Tracking Errors', 'FontSize', 20);
    legend('Pointing error (deg)', 'Body-rate error (deg/s)', 'Location', 'best', 'FontSize', 14);
    grid on;
    xlim([tspan(tracking_indices(1)) / 60, tspan(tracking_indices(end)) / 60]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.0);
    postMatchStem = sprintf('too_minimax_post_match_error_seed_%s_aero_%s', ...
        scenarioSeedText, vleo.util.on_off_text(params.includeAerodynamics));
    postMatchExport = vleo.util.export_paper_figure(gcf, postMatchStem);
    fprintf('Paper figure saved to: %s\n', postMatchExport.pngPath);
    fprintf('Paper figure saved to: %s\n', postMatchExport.pdfPath);
else
    fprintf('Skipping GUI animation and figure generation (generateVisuals=false).\n');
end

fprintf('\nSimulation complete!\n');

function mark_phase_lines(ax, tStartVisible, tEruption, tMatch, tEndVisible, fontSize)
    hold(ax, 'on');
    if ~isnan(tStartVisible)
        marker = xline(ax, tStartVisible / 60, 'g-', 'LineWidth', 1.5, 'Label', 'Start Visible', 'FontSize', fontSize);
        vleo.util.hide_from_legend(marker);
    end
    marker = xline(ax, tEruption / 60, 'k--', 'LineWidth', 2, 'Label', 'TOO Event', 'FontSize', fontSize);
    vleo.util.hide_from_legend(marker);
    marker = xline(ax, tMatch / 60, 'c-.', 'LineWidth', 2, 'Label', 'State Match', 'FontSize', fontSize);
    vleo.util.hide_from_legend(marker);
    if ~isnan(tEndVisible)
        marker = xline(ax, tEndVisible / 60, 'm-', 'LineWidth', 1.5, 'Label', 'Out of Sight', 'FontSize', fontSize);
        vleo.util.hide_from_legend(marker);
    end
end

function tauBody = interpolate_maneuver_then_tracking_torque(t, maneuverTimeNodes, maneuverTauNodes, ...
        tMatch, trackingTimes, trackingTauHistory, trackingInterpMethod)
    if t <= tMatch
        tauBody = interp1(maneuverTimeNodes, maneuverTauNodes, t, 'previous', 'extrap')';
        return;
    end

    tauBody = interp1(trackingTimes, trackingTauHistory, t, trackingInterpMethod, 'extrap')';
end

function selection = select_robust_match_time(tspan, idxEruption, idxEndVisible, trackingHistory, ...
        qInitial, omegaInitialBody, inertiaBody, momentArms, tMatchDefault)
    dt = max(eps, abs(tspan(2) - tspan(1)));
    candidateBufferSec = 5;
    candidateStepSec = 5;
    candidateBufferSteps = max(1, round(candidateBufferSec / dt));
    candidateStep = max(1, round(candidateStepSec / dt));
    [~, idxDefault] = min(abs(tspan - tMatchDefault));

    idxStart = min(idxEndVisible - 1, idxEruption + candidateBufferSteps);
    idxCandidates = idxStart:candidateStep:(idxEndVisible - candidateBufferSteps);
    idxCandidates = unique([idxCandidates, idxDefault]);
    idxCandidates = idxCandidates(idxCandidates > idxEruption & idxCandidates < idxEndVisible);

    nCandidates = numel(idxCandidates);
    gammaHistory_N = nan(nCandidates, 1);
    angleErrorDegHistory = inf(nCandidates, 1);
    rateErrorDegPerSecHistory = inf(nCandidates, 1);
    isFeasibleHistory = false(nCandidates, 1);
    qTargets = zeros(4, nCandidates);
    omegaTargets = zeros(3, nCandidates);
    results = cell(nCandidates, 1);

    for candidateIdx = 1:nCandidates
        idxMatch = idxCandidates(candidateIdx);
        qTarget = trackingHistory.q_track_history(idxMatch, :)';
        qTarget = qTarget / norm(qTarget);
        if dot(qInitial, qTarget) < 0
            qTarget = -qTarget;
        end
        omegaTargetBody = trackingHistory.omega_track_body_history(idxMatch, :)';

        result = vleo.control.solve_minimax_attitude_slew_fast_legacy(qInitial, omegaInitialBody, qTarget, ...
            omegaTargetBody, inertiaBody, tspan(idxMatch) - tspan(idxEruption), 20, ...
            'refineNonlinear', true, 'nRefinementIntervals', [12, 16, 20], ...
            'momentArms', momentArms);

        qTargets(:, candidateIdx) = qTarget;
        omegaTargets(:, candidateIdx) = omegaTargetBody;
        results{candidateIdx} = result;
        gammaHistory_N(candidateIdx) = result.gammaForceEquivalent;
        angleErrorDegHistory(candidateIdx) = result.terminalAngleErrorDeg;
        rateErrorDegPerSecHistory(candidateIdx) = result.terminalRateErrorDegPerSec;
        isFeasibleHistory(candidateIdx) = is_match_result_feasible(result);
    end

    feasibleCandidates = find(isFeasibleHistory);
    if any(idxCandidates(feasibleCandidates) == idxDefault)
        selectedCandidateIdx = find(idxCandidates == idxDefault, 1, 'first');
        selectionText = 'default 3/4 visibility rule';
    elseif ~isempty(feasibleCandidates)
        feasibleTimes = tspan(idxCandidates(feasibleCandidates));
        sortTable = [abs(feasibleTimes(:) - tMatchDefault), gammaHistory_N(feasibleCandidates)];
        [~, order] = sortrows(sortTable, [1, 2]);
        selectedCandidateIdx = feasibleCandidates(order(1));
        selectionText = 'auto-selected near default 3/4 visibility rule';
    else
        penalty = angleErrorDegHistory + rateErrorDegPerSecHistory;
        [~, selectedCandidateIdx] = min(penalty);
        selectionText = 'fallback lowest-error handoff candidate';
    end

    selection = struct();
    selection.idxCandidates = idxCandidates;
    selection.candidateTimes = tspan(idxCandidates);
    selection.gammaHistory_N = gammaHistory_N;
    selection.angleErrorDegHistory = angleErrorDegHistory;
    selection.rateErrorDegPerSecHistory = rateErrorDegPerSecHistory;
    selection.isFeasibleHistory = isFeasibleHistory;
    selection.idxSelected = idxCandidates(selectedCandidateIdx);
    selection.selectedTime = tspan(idxCandidates(selectedCandidateIdx));
    selection.defaultTime = tMatchDefault;
    selection.defaultSampledTime = tspan(idxDefault);
    selection.selectionText = selectionText;
    selection.qTarget = qTargets(:, selectedCandidateIdx);
    selection.omegaTargetBody = omegaTargets(:, selectedCandidateIdx);
    selection.resultHistory = results;
    selection.resultSelected = results{selectedCandidateIdx};
end

function selection = evaluate_match_time_scan_with_aero(selection, tspan, idxStartVisible, idxEruption, ...
        idxEndVisible, XOrbitHistory, XEruption, trackingHistory, params, odeOpts)
    if ~isempty(idxStartVisible)
        windowIndices = idxStartVisible:idxEndVisible;
    else
        windowIndices = idxEruption:idxEndVisible;
    end

    windowTimes = tspan(windowIndices);
    nCandidates = numel(selection.candidateTimes);
    peakActuatorForceWithAero_N = nan(nCandidates, 1);

    for candidateIdx = 1:nCandidates
        idxMatch = selection.idxCandidates(candidateIdx);
        result = selection.resultHistory{candidateIdx};

        tauTotalWindow = zeros(numel(windowIndices), 3);
        maneuverIndices = idxEruption:idxMatch;
        maneuverInterpMethod = 'pchip';
        if isfield(result, 'interpMethod')
            maneuverInterpMethod = result.interpMethod;
        end
        tauManeuver = interp1(result.timeNodes, result.tauNodes, ...
            tspan(maneuverIndices) - tspan(idxEruption), maneuverInterpMethod);
        localManeuverIndices = maneuverIndices - windowIndices(1) + 1;
        tauTotalWindow(localManeuverIndices, :) = tauManeuver;

        if idxMatch < idxEndVisible
            trackingIndices = idxMatch + 1:idxEndVisible;
            localTrackingIndices = trackingIndices - windowIndices(1) + 1;
            tauTotalWindow(localTrackingIndices, :) = trackingHistory.tau_track_total_history(trackingIndices, :);
        end

        tControlCandidate = tspan(idxEruption:idxEndVisible);
        maneuverTimeNodes = tspan(idxEruption) + result.timeNodes;
        tauTotalInterp = @(t) interpolate_maneuver_then_tracking_torque(t, maneuverTimeNodes, ...
            result.tauNodes, tspan(idxMatch), tspan, trackingHistory.tau_track_total_history, 'previous');
        zeroTorqueInterp = @(t) zeros(3, 1);
        odeFunTotal = @(t, X) vleo.dynamics.sat_dynamics_openloop(t, X, params, tauTotalInterp, zeroTorqueInterp);
        [~, XPostTotal] = ode45(odeFunTotal, tControlCandidate, XEruption, odeOpts);

        XDesiredWindow = XOrbitHistory(windowIndices, :);
        controlIndices = idxEruption:idxEndVisible;
        localControlIndices = controlIndices - windowIndices(1) + 1;
        XDesiredWindow(localControlIndices, :) = XPostTotal;

        trackingIndices = idxMatch:idxEndVisible;
        localTrackingIndices = trackingIndices - windowIndices(1) + 1;
        XDesiredWindow(localTrackingIndices, 7:10) = trackingHistory.q_track_history(trackingIndices, :);
        XDesiredWindow(localTrackingIndices, 11:13) = trackingHistory.omega_track_body_history(trackingIndices, :);
        XDesiredWindow(:, 7:10) = vleo.util.normalize_quaternion_rows(XDesiredWindow(:, 7:10));
        XDesiredWindow(:, 7:10) = vleo.util.align_quaternion_signs(XDesiredWindow(:, 7:10));

        RDesiredWindow = vleo.util.body_to_eci_history_from_quaternion_rows(XDesiredWindow(:, 7:10));
        if params.includeAerodynamics
            tauAeroWindow = vleo.dynamics.estimate_aero_torque_history(windowTimes, XDesiredWindow, RDesiredWindow, params);
        else
            tauAeroWindow = zeros(numel(windowIndices), 3);
        end

        tauControlWindow = tauTotalWindow - tauAeroWindow;
        peakActuatorForceWithAero_N(candidateIdx) = max(max(abs(tauControlWindow) ./ params.momentArms'));
    end

    selection.windowIndices = windowIndices;
    selection.peakActuatorForceWithAero_N = peakActuatorForceWithAero_N;
    selectedMask = abs(selection.candidateTimes(:) - selection.selectedTime) < 1e-12;
    selection.peakActuatorForceWithAeroSelected_N = peakActuatorForceWithAero_N(selectedMask);
end

function isFeasible = is_match_result_feasible(result)
    isFeasible = result.terminalAngleErrorDeg <= 5e-2 && ...
        result.terminalRateErrorDegPerSec <= 5e-3 && ...
        result.maxInequalityViolation <= 1e-10;
end
