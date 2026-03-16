function history = compute_observation_history(tspan, xHistory, rVolcano0, params)
    nSteps = size(xHistory, 1);
    history = struct();
    history.ra_camera_history = zeros(nSteps, 1);
    history.dec_camera_history = zeros(nSteps, 1);
    history.ra_volcano_history = zeros(nSteps, 1);
    history.dec_volcano_history = zeros(nSteps, 1);
    history.visibility_history = false(nSteps, 1);
    history.zenith_distance_history = zeros(nSteps, 1);
    history.pointing_eci_history = zeros(nSteps, 3);
    history.R_B2E_hist = zeros(3, 3, nSteps);
    history.q_track_history = zeros(nSteps, 4);
    history.euler_track_history_deg = zeros(nSteps, 3);

    for i = 1:nSteps
        obs = vleo.analysis.state_to_observation(xHistory(i, :)', params);
        history.ra_camera_history(i) = obs.ra;
        history.dec_camera_history(i) = obs.dec;
        history.pointing_eci_history(i, :) = obs.pointing_eci';

        rSat = xHistory(i, 1:3)';
        vSat = xHistory(i, 4:6)';
        [rVolcano, vVolcano] = vleo.dynamics.volcano_state_at_time(tspan(i), rVolcano0, params.omega_earth);
        apparentVolcanoVec = rVolcano - rSat;
        apparentVolcanoVec = apparentVolcanoVec / norm(apparentVolcanoVec);

        decRad = asin(apparentVolcanoVec(3));
        raRad = atan2(apparentVolcanoVec(2), apparentVolcanoVec(1));
        if raRad < 0
            raRad = raRad + 2 * pi;
        end
        history.ra_volcano_history(i) = rad2deg(raRad);
        history.dec_volcano_history(i) = rad2deg(decRad);

        history.visibility_history(i) = dot(rVolcano, rSat - rVolcano) > 0;

        zenithAtVolcano = rVolcano / norm(rVolcano);
        volcanoToSat = rSat - rVolcano;
        volcanoToSatUnit = volcanoToSat / norm(volcanoToSat);
        cosZenith = dot(zenithAtVolcano, volcanoToSatUnit);
        history.zenith_distance_history(i) = rad2deg(acos(vleo.util.clamp_scalar(cosZenith, -1, 1)));

        history.R_B2E_hist(:, :, i) = vleo.dynamics.build_tracking_frame(rSat, vSat, rVolcano, vVolcano);
        history.q_track_history(i, :) = dcm2quat(history.R_B2E_hist(:, :, i)')';
        history.euler_track_history_deg(i, :) = vleo.util.rotation_matrix_to_euler_deg(history.R_B2E_hist(:, :, i));
    end

    history.q_track_history = vleo.util.normalize_quaternion_rows(history.q_track_history);
    history.q_track_history = vleo.util.align_quaternion_signs(history.q_track_history);

    dt = 1;
    if numel(tspan) > 1
        dt = abs(tspan(2) - tspan(1));
    end
    history.omega_track_body_history = vleo.dynamics.compute_body_rate_history_from_dcm_history(history.R_B2E_hist, dt);
    history.alpha_track_body_history = vleo.dynamics.compute_angular_acceleration_history(history.omega_track_body_history, dt);
    history.tau_track_total_history = zeros(nSteps, 3);
    history.omega_track_eci_history = zeros(nSteps, 3);

    for i = 1:nSteps
        omegaBody = history.omega_track_body_history(i, :)';
        alphaBody = history.alpha_track_body_history(i, :)';
        history.tau_track_total_history(i, :) = vleo.dynamics.rigid_body_total_torque(omegaBody, alphaBody, params.I_CB)';
        history.omega_track_eci_history(i, :) = (history.R_B2E_hist(:, :, i) * omegaBody)';
    end
end
