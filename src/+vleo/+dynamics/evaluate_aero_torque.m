function tauAero = evaluate_aero_torque(t, X, rEciToBody, params)
    tauAero = [0; 0; 0];

    rEci = X(1:3);
    vEci = X(4:6);
    if ~(all(isfinite(rEci)) && all(isfinite(vEci)) && norm(rEci) > 0)
        return;
    end

    vBody = rEciToBody * vEci;
    if norm(vBody) > 1
        attitude.aoa = atan2d(vBody(3), vBody(1));
        attitude.aos = atan2d(vBody(2), vBody(1));
    else
        attitude.aoa = 0;
        attitude.aos = 0;
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
