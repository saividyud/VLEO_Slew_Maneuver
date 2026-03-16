function solverCase = build_volcano_minimax_slew_case(scenarioSeed)
    if ~(isscalar(scenarioSeed) && isnumeric(scenarioSeed) && isfinite(scenarioSeed) && ...
            scenarioSeed >= 0 && scenarioSeed == floor(scenarioSeed))
        error('VLEO_Slew_Maneuver:InvalidScenarioSeed', ...
            'scenarioSeed must be a nonnegative integer scalar.');
    end

    rng(double(scenarioSeed), 'twister');
    params = vleo.dynamics.default_control_test_params(false);
    scenario = vleo.analysis.generate_volcano_scenario(params);

    odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    odeFun = @(t, X) vleo.dynamics.sat_dynamics_controlled(t, X, scenario.z_body_eci, params, ...
        'useJ2', false, 'useAtmDrag', false, 'useControl', false);
    [tspan, xOrbitHistory] = vleo.dynamics.propagate_bidirectional(odeFun, scenario.X0, ...
        scenario.t_start, scenario.t_end, scenario.dt, odeOpts);

    xOrbitHistory(:, 7:10) = vleo.util.normalize_quaternion_rows(xOrbitHistory(:, 7:10));
    xOrbitHistory(:, 7:10) = vleo.util.align_quaternion_signs(xOrbitHistory(:, 7:10));

    trackingHistory = vleo.analysis.compute_observation_history(tspan, xOrbitHistory, scenario.r_volcano_0, params);
    visibilityInfo = vleo.analysis.analyze_visibility(tspan, trackingHistory.visibility_history);
    idxEruption = find(abs(tspan - scenario.t_eruption) < 1e-12, 1, 'first');

    if isempty(visibilityInfo.idx_end_visible) || visibilityInfo.t_end_visible <= scenario.t_eruption
        error('control_test4:InvalidVisibilityWindow', ...
            'The volcano is not visible long enough after eruption to define a solver case.');
    end

    tMatchRequest = scenario.t_eruption + 0.75 * (visibilityInfo.t_end_visible - scenario.t_eruption);
    [~, idxMatch] = min(abs(tspan - tMatchRequest));
    if idxMatch <= idxEruption || idxMatch >= visibilityInfo.idx_end_visible
        error('control_test4:SnappedMatchTimeOutOfRange', ...
            'The sampled match time is outside the valid post-eruption visibility interval.');
    end

    xEruption = xOrbitHistory(idxEruption, :)';
    qInitial = xEruption(7:10) / norm(xEruption(7:10));
    omegaInitialBody = xEruption(11:13);
    qTarget = trackingHistory.q_track_history(idxMatch, :)';
    qTarget = qTarget / norm(qTarget);
    if dot(qInitial, qTarget) < 0
        qTarget = -qTarget;
    end
    omegaTargetBody = trackingHistory.omega_track_body_history(idxMatch, :)';

    solverCase = struct();
    solverCase.caseId = sprintf('volcano_minimax_slew_seed_%d', scenarioSeed);
    solverCase.scenarioSeed = double(scenarioSeed);
    solverCase.tStartVisibleSec = visibilityInfo.t_start_visible;
    solverCase.tEndVisibleSec = visibilityInfo.t_end_visible;
    solverCase.tEruptionSec = scenario.t_eruption;
    solverCase.tMatchSec = tspan(idxMatch);
    solverCase.maneuverDurationSec = tspan(idxMatch) - scenario.t_eruption;
    solverCase.qInitial = qInitial;
    solverCase.qTarget = qTarget;
    solverCase.eulerInitialDeg = vleo.util.quat_to_euler_deg(qInitial');
    solverCase.eulerTargetDeg = vleo.util.quat_to_euler_deg(qTarget');
    solverCase.omegaInitialBody = omegaInitialBody;
    solverCase.omegaTargetBody = omegaTargetBody;
    solverCase.omegaInitialDegPerSec = rad2deg(omegaInitialBody);
    solverCase.omegaTargetDegPerSec = rad2deg(omegaTargetBody);
    solverCase.relativeAngleDeg = vleo.util.rotation_error_angle_deg(qInitial, qTarget);
    solverCase.targetRaDeg = trackingHistory.ra_volcano_history(idxMatch);
    solverCase.targetDecDeg = trackingHistory.dec_volcano_history(idxMatch);
end
