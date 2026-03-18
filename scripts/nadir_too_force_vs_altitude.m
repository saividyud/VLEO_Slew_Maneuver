% Nadir-Pass TOO Tracking: Max Actuator Force vs Altitude
% Forces the TOO to pass directly through nadir, computes the tracking
% torque with aerodynamic compensation, converts to actuator force, and
% sweeps altitude to find where the peak force crosses 145 mN.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc;

%% Configuration
momentArms = [0.3; 0.3; 0.2];
forceLimit_N = 0.145;
altitudes_km = 100:10:200;
nAltitudes = numel(altitudes_km);

fprintf('=== NADIR-PASS TOO TRACKING: FORCE vs ALTITUDE ===\n');
fprintf('Moment arms [m]: [%.2f, %.2f, %.2f]\n', momentArms);
fprintf('Force limit: %.1f mN\n', forceLimit_N * 1e3);
fprintf('Altitude range: %d to %d km (%d points)\n', altitudes_km(1), altitudes_km(end), nAltitudes);
fprintf('Orbit: circular, incl = 51.6 deg, fixed RAAN/TA\n');
fprintf('TOO placement: directly at sub-satellite point (nadir) at t=0\n\n');

%% Preallocate
maxForce_mN = zeros(nAltitudes, 1);
maxForceAxis = zeros(nAltitudes, 1);
maxTotalTorque_mNm = zeros(nAltitudes, 1);
maxAeroTorque_mNm = zeros(nAltitudes, 1);
maxControlTorque_mNm = zeros(nAltitudes, 1);
visibleDuration_s = zeros(nAltitudes, 1);

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
sweepTimer = tic;

%% Altitude sweep
for altIdx = 1:nAltitudes
    altitude_m = altitudes_km(altIdx) * 1e3;
    params = vleo.dynamics.default_control_test_params(true);
    scenario = build_nadir_scenario(params, altitude_m);

    %% Propagate orbit (Keplerian, no control)
    odeFun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, scenario.z_body_eci, params, ...
        'useJ2', false, 'useAtmDrag', false, 'useControl', false);
    [tspan, X_history] = vleo.dynamics.propagate_bidirectional(odeFun, scenario.X0, ...
        scenario.t_start, scenario.t_end, scenario.dt, odeOpts);

    %% Compute tracking history (desired torques to follow TOO)
    trackingHistory = vleo.analysis.compute_observation_history(tspan, X_history, scenario.r_too_0, params);
    tau_total_history = trackingHistory.tau_track_total_history;

    %% Estimate aerodynamic torque in tracking attitude
    tau_aero_history = vleo.dynamics.estimate_aero_torque_history( ...
        tspan, X_history, trackingHistory.R_B2E_hist, params);

    %% Actuator (control) torque = total required - aero contribution
    tau_control_history = tau_total_history - tau_aero_history;

    %% Restrict to visible window
    visibilityInfo = vleo.analysis.analyze_visibility(tspan, trackingHistory.visibility_history);
    if ~isempty(visibilityInfo.idx_start_visible) && ~isempty(visibilityInfo.idx_end_visible)
        visIdx = visibilityInfo.idx_start_visible:visibilityInfo.idx_end_visible;
        visibleDuration_s(altIdx) = tspan(visIdx(end)) - tspan(visIdx(1));
    else
        visIdx = find(trackingHistory.visibility_history);
        if ~isempty(visIdx)
            visibleDuration_s(altIdx) = tspan(visIdx(end)) - tspan(visIdx(1));
        end
    end

    if isempty(visIdx)
        maxForce_mN(altIdx) = NaN;
        fprintf('[%02d/%02d] %3d km: NO VISIBILITY\n', altIdx, nAltitudes, altitudes_km(altIdx));
        continue;
    end

    %% Convert control torque to force per axis, find peak
    forceHistory = abs(tau_control_history(visIdx, :)) ./ momentArms';
    [peakForce, peakLinearIdx] = max(forceHistory(:));
    [~, peakAxisIdx] = ind2sub(size(forceHistory), peakLinearIdx);
    maxForce_mN(altIdx) = peakForce * 1e3;
    maxForceAxis(altIdx) = peakAxisIdx;

    maxTotalTorque_mNm(altIdx) = max(vecnorm(tau_total_history(visIdx, :), 2, 2)) * 1e3;
    maxAeroTorque_mNm(altIdx) = max(vecnorm(tau_aero_history(visIdx, :), 2, 2)) * 1e3;
    maxControlTorque_mNm(altIdx) = max(vecnorm(tau_control_history(visIdx, :), 2, 2)) * 1e3;

    axisLabels = {'x', 'y', 'z'};
    fprintf('[%02d/%02d] %3d km: F_max=%.1f mN (axis %s)  tau_ctrl=%.2f mNm  tau_aero=%.2f mNm  vis=%.0f s\n', ...
        altIdx, nAltitudes, altitudes_km(altIdx), maxForce_mN(altIdx), ...
        axisLabels{peakAxisIdx}, maxControlTorque_mNm(altIdx), ...
        maxAeroTorque_mNm(altIdx), visibleDuration_s(altIdx));
end

elapsedMin = toc(sweepTimer) / 60;
fprintf('\nSweep completed in %.1f minutes.\n', elapsedMin);

%% Find crossing altitude — coarse bracket then bisection refinement
validMask = ~isnan(maxForce_mN);
altValid = altitudes_km(validMask);
forceValid = maxForce_mN(validMask);
forceLimit_mN = forceLimit_N * 1e3;

bracketLo = NaN;
bracketHi = NaN;
for k = 1:numel(forceValid) - 1
    if (forceValid(k) - forceLimit_mN) * (forceValid(k + 1) - forceLimit_mN) < 0
        bracketLo = altValid(k);
        bracketHi = altValid(k + 1);
        break;
    end
end

crossingAlt_km = NaN;
crossingForce_mN = NaN;
if isnan(bracketLo)
    if all(forceValid > forceLimit_mN)
        fprintf('Max force exceeds %.0f mN at ALL altitudes in the sweep range.\n', forceLimit_mN);
    elseif all(forceValid < forceLimit_mN)
        fprintf('Max force is below %.0f mN at ALL altitudes in the sweep range.\n', forceLimit_mN);
    else
        fprintf('Could not find a bracket for the crossing.\n');
    end
else
    fprintf('\n--- Bisection refinement between %.0f and %.0f km ---\n', bracketLo, bracketHi);
    bisectionTol_km = 0.1;
    maxBisectionIter = 20;

    fLo = evaluate_max_force_at_altitude(bracketLo, params, momentArms, odeOpts);
    fHi = evaluate_max_force_at_altitude(bracketHi, params, momentArms, odeOpts);
    fprintf('  bracket: [%.1f km -> %.1f mN]  [%.1f km -> %.1f mN]\n', bracketLo, fLo, bracketHi, fHi);

    for bisIter = 1:maxBisectionIter
        altMid = 0.5 * (bracketLo + bracketHi);
        fMid = evaluate_max_force_at_altitude(altMid, params, momentArms, odeOpts);
        fprintf('  iter %2d: %.2f km -> %.2f mN\n', bisIter, altMid, fMid);

        if (fLo - forceLimit_mN) * (fMid - forceLimit_mN) < 0
            bracketHi = altMid;
            fHi = fMid;
        else
            bracketLo = altMid;
            fLo = fMid;
        end

        if (bracketHi - bracketLo) < bisectionTol_km
            break;
        end
    end

    crossingAlt_km = 0.5 * (bracketLo + bracketHi);
    crossingForce_mN = 0.5 * (fLo + fHi);
    fprintf('\n>>> Force crosses %.0f mN at %.2f km altitude (refined to %.1f km bracket) <<<\n', ...
        forceLimit_mN, crossingAlt_km, bracketHi - bracketLo);
end

%% Plot
figure('Name', 'Nadir TOO Tracking Force vs Altitude', ...
    'Color', 'w', 'Position', [100, 100, 1100, 650]);

plot(altitudes_km(validMask), forceValid, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'Max actuator force (with aero)');
hold on;
yline(forceLimit_mN, 'r--', 'LineWidth', 1.8, 'Label', sprintf('%.0f mN limit', forceLimit_mN), ...
    'FontSize', 13, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');

if ~isnan(crossingAlt_km)
    xline(crossingAlt_km, 'k-.', 'LineWidth', 1.5, ...
        'Label', sprintf('%.1f km', crossingAlt_km), 'FontSize', 13, ...
        'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
    plot(crossingAlt_km, forceLimit_mN, 'rp', 'MarkerSize', 14, ...
        'MarkerFaceColor', 'r', 'DisplayName', sprintf('Crossing: %.1f km', crossingAlt_km));
end

hold off;
grid on;
xlabel('Altitude (km)', 'FontSize', 16);
ylabel('Peak Actuator Force (mN)', 'FontSize', 16);
title('Nadir-Pass TOO Tracking: Max Force vs Altitude (with Aerodynamics)', 'FontSize', 18);
legend('Location', 'northeast', 'FontSize', 13);
xlim([altitudes_km(1), altitudes_km(end)]);

fprintf('\nDone.\n');

function maxForce_mN = evaluate_max_force_at_altitude(altitude_km, params, momentArms, odeOpts)
    altitude_m = altitude_km * 1e3;
    scenario = build_nadir_scenario(params, altitude_m);

    odeFun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, scenario.z_body_eci, params, ...
        'useJ2', false, 'useAtmDrag', false, 'useControl', false);
    [tspan, X_history] = vleo.dynamics.propagate_bidirectional(odeFun, scenario.X0, ...
        scenario.t_start, scenario.t_end, scenario.dt, odeOpts);

    trackingHistory = vleo.analysis.compute_observation_history(tspan, X_history, scenario.r_too_0, params);
    tau_aero_history = vleo.dynamics.estimate_aero_torque_history( ...
        tspan, X_history, trackingHistory.R_B2E_hist, params);
    tau_control_history = trackingHistory.tau_track_total_history - tau_aero_history;

    visibilityInfo = vleo.analysis.analyze_visibility(tspan, trackingHistory.visibility_history);
    if ~isempty(visibilityInfo.idx_start_visible) && ~isempty(visibilityInfo.idx_end_visible)
        visIdx = visibilityInfo.idx_start_visible:visibilityInfo.idx_end_visible;
    else
        visIdx = find(trackingHistory.visibility_history);
    end

    if isempty(visIdx)
        maxForce_mN = NaN;
        return;
    end

    forceHistory = abs(tau_control_history(visIdx, :)) ./ momentArms';
    maxForce_mN = max(forceHistory(:)) * 1e3;
end

function scenario = build_nadir_scenario(params, altitude_m)
    rOrbit = params.R_e + altitude_m;
    vCirc = sqrt(params.mu / rOrbit);

    inclination = deg2rad(51.6);
    raan = 0;
    trueAnomaly = 0;

    rRaan = [cos(raan), -sin(raan), 0; sin(raan), cos(raan), 0; 0, 0, 1];
    rInc = [1, 0, 0; 0, cos(inclination), -sin(inclination); 0, sin(inclination), cos(inclination)];
    rTa = [cos(trueAnomaly), -sin(trueAnomaly), 0; sin(trueAnomaly), cos(trueAnomaly), 0; 0, 0, 1];
    rPerifocalToEci = rRaan * rInc * rTa;

    rPerifocal = [rOrbit; 0; 0];
    vPerifocal = [0; vCirc; 0];
    rEci = rPerifocalToEci * rPerifocal;
    vEci = rPerifocalToEci * vPerifocal;

    orbitNormalEci = cross(rEci, vEci) / norm(cross(rEci, vEci));
    zBodyEci = -rEci / norm(rEci);
    xBodyEci = vEci / norm(vEci);
    yBodyEci = cross(zBodyEci, xBodyEci);
    rBodyToEci = [xBodyEci, yBodyEci, zBodyEci];
    q0 = dcm2quat(rBodyToEci')';

    omegaOrbitMag = vCirc / rOrbit;
    omegaEci = omegaOrbitMag * orbitNormalEci;
    omegaBody = rBodyToEci' * omegaEci;
    x0 = [rEci; vEci; q0; omegaBody];

    orbitalPeriod = 2 * pi * sqrt(rOrbit^3 / params.mu);

    rTOO0 = rEci / norm(rEci) * params.R_e;

    scenario = struct();
    scenario.altitude = altitude_m;
    scenario.r_orbit = rOrbit;
    scenario.v_circ = vCirc;
    scenario.r_eci = rEci;
    scenario.v_eci = vEci;
    scenario.z_body_eci = zBodyEci;
    scenario.X0 = x0;
    scenario.orbital_period = orbitalPeriod;
    scenario.r_too_0 = rTOO0;
    scenario.theta_max = acos(params.R_e / rOrbit);
    scenario.theta = 0;
    scenario.phi = 0;
    scenario.visibility_at_eruption = dot(rTOO0, x0(1:3) - rTOO0);
    scenario.t_start = -orbitalPeriod / 2;
    scenario.t_end = orbitalPeriod / 2;
    scenario.t_eruption = 0;
    scenario.dt = 1;
end
