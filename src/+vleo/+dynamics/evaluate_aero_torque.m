function tauAero = evaluate_aero_torque(t, X, rEciToBody, params)
    tauAero = [0; 0; 0];

    rEci = X(1:3);
    vEci = X(4:6);
    if ~(all(isfinite(rEci)) && all(isfinite(vEci)) && norm(rEci) > 0)
        return;
    end

    vBody = rEciToBody * vEci;
    if norm(vBody) > 1
        % Compute unit vector in wind frame (vBody direction)
        % +X_wind is aligned with the spacecraft velocity vector relative to the air
        u_vel_body = vBody / norm(vBody);
        wind_X = u_vel_body;
        
        % Construct a wind frame
        % Default +Z wind is loosely +Z body, but orthogonalized
        if abs(wind_X(3)) < 0.99
            wind_Y = cross(wind_X, [0; 0; 1]);
        else
            wind_Y = cross(wind_X, [1; 0; 0]);
        end
        wind_Y = wind_Y / norm(wind_Y);
        wind_Z = cross(wind_X, wind_Y);
        
        R_body2wind = [wind_X, wind_Y, wind_Z]'; % DCM from body to wind
        attitude.quaternion = dcm2quat(R_body2wind)'; % [q0; q1; q2; q3]
    else
        attitude.quaternion = [1; 0; 0; 0];
    end

    location.positionECI = rEci;
    totalSeconds = params.aero.secondsOfDay + t;
    time.year = params.aero.year;
    time.dayOfYear = params.aero.dayOfYear + floor(totalSeconds / 86400);
    time.UTseconds = mod(totalSeconds, 86400);

    try
        aeroResults = vleo.aero.compute_aero_forces(params.aero.objFilePath, location, attitude, time, params.aero.options);
        tauCandidate = reshape(aeroResults.M_aero, [], 1);
        if isfield(aeroResults, 'M_solar')
            tauCandidate = tauCandidate + reshape(aeroResults.M_solar, [], 1);
        end
        if numel(tauCandidate) == 3 && all(isfinite(tauCandidate))
            tauAero = tauCandidate;
        end
    catch
        tauAero = [0; 0; 0];
    end
end
