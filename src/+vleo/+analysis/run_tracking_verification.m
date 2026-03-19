function verification = run_tracking_verification(tspan, xHistory, trackingHistory, visibilityInfo, ...
        rTOO0, tauControlHistory, tauAeroHistory, params, odeOpts)
    verification = struct( ...
        'plotTime', tspan(:), ...
        'plotDesiredEulerDeg', trackingHistory.euler_track_history_deg, ...
        'plotActualEulerDeg', vleo.util.quat_history_to_euler_deg(xHistory(:, 7:10)), ...
        'ra_error', [], ...
        'dec_error', [], ...
        'pointing_error_deg', [], ...
        'tspan_verif', [], ...
        'X_verif', []);

    if isnan(visibilityInfo.t_start_visible) || isnan(visibilityInfo.t_end_visible)
        return;
    end

    idxStart = visibilityInfo.idx_start_visible;
    idxEnd = visibilityInfo.idx_end_visible;
    t0Verif = tspan(idxStart);
    rSat0 = xHistory(idxStart, 1:3)';
    vSat0 = xHistory(idxStart, 4:6)';

    [rTOOAtTime0, vTOOAtTime0] = vleo.dynamics.too_state_at_time(t0Verif, rTOO0, params.omega_earth);
    rRel0 = rTOOAtTime0 - rSat0;
    vRel0 = vTOOAtTime0 - vSat0;

    zBody0 = rRel0 / norm(rRel0);
    dist0 = norm(rRel0);
    zBodyDot0 = (vRel0 / dist0) - zBody0 * (dot(rRel0, vRel0) / dist0^2);

    yBody0 = cross(zBody0, zBodyDot0);
    if norm(yBody0) < 1e-10
        yBody0 = cross(zBody0, vSat0);
        if norm(yBody0) < 1e-10
            yBody0 = cross(zBody0, rSat0);
        end
    end
    yBody0 = yBody0 / norm(yBody0);

    xBody0 = cross(yBody0, zBody0);
    rEciToBody0 = [xBody0, yBody0, zBody0]';
    q0Verif = dcm2quat(rEciToBody0)';
    omega0VerifBody = trackingHistory.omega_track_body_history(idxStart, :)';
    x0Verif = [rSat0; vSat0; q0Verif; omega0VerifBody];

    tauControlInterp = @(t) interp1(tspan, tauControlHistory, t, 'pchip', 'extrap')';
    tauAeroInterp = @(t) interp1(tspan, tauAeroHistory, t, 'pchip', 'extrap')';
    odeFunVerif = @(t, X) vleo.dynamics.sat_dynamics_openloop(t, X, params, tauControlInterp, tauAeroInterp);
    verification.tspan_verif = tspan(idxStart:idxEnd);
    [~, verification.X_verif] = ode45(odeFunVerif, verification.tspan_verif, x0Verif, odeOpts);

    verification.plotTime = verification.tspan_verif(:);
    verification.plotDesiredEulerDeg = trackingHistory.euler_track_history_deg(idxStart:idxEnd, :);
    verification.plotActualEulerDeg = vleo.util.quat_history_to_euler_deg(verification.X_verif(:, 7:10));

    raVerifHistory = zeros(numel(verification.tspan_verif), 1);
    decVerifHistory = zeros(numel(verification.tspan_verif), 1);
    verification.pointing_error_deg = zeros(numel(verification.tspan_verif), 1);
    for k = 1:numel(verification.tspan_verif)
        obsVerif = vleo.analysis.state_to_observation(verification.X_verif(k, :)', params);
        raVerifHistory(k) = obsVerif.ra;
        decVerifHistory(k) = obsVerif.dec;

        rSatVerif = verification.X_verif(k, 1:3)';
        [rTOO, ~] = vleo.dynamics.too_state_at_time(verification.tspan_verif(k), rTOO0, params.omega_earth);
        apparentTOOVec = rTOO - rSatVerif;
        apparentTOOVec = apparentTOOVec / norm(apparentTOOVec);
        verification.pointing_error_deg(k) = rad2deg(acos(vleo.util.clamp_scalar( ...
            dot(obsVerif.pointing_eci, apparentTOOVec), -1, 1)));
    end

    raTargetHistory = trackingHistory.ra_TOO_history(idxStart:idxEnd);
    decTargetHistory = trackingHistory.dec_TOO_history(idxStart:idxEnd);
    verification.ra_error = mod(raVerifHistory - raTargetHistory + 180, 360) - 180;
    verification.dec_error = decVerifHistory - decTargetHistory;
end
