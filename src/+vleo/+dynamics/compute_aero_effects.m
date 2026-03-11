function [aAeroEci, mAeroBody, aeroResults] = compute_aero_effects(t, X, config)
    aAeroEci = [0; 0; 0];
    mAeroBody = [0; 0; 0];
    aeroResults = struct();

    location.positionECI = X(1:3);

    qScalarFirst = reshape(X(7:10), 1, []);
    vEci = X(4:6);
    rEciToBody = quat2dcm(qScalarFirst);
    vBody = rEciToBody * vEci;

    if norm(vBody) > 1
        attitude.aoa = atan2d(vBody(3), vBody(1));
        attitude.aos = atan2d(vBody(2), vBody(1));
    else
        attitude.aoa = 0;
        attitude.aos = 0;
    end

    totalSeconds = config.secondsOfDay + t;
    time.year = config.year;
    time.dayOfYear = config.dayOfYear + floor(totalSeconds / 86400);
    time.UTseconds = mod(totalSeconds, 86400);

    try
        aeroResults = vleo.aero.compute_aero_forces( ...
            config.objFilePath, location, attitude, time, config.options);

        forceAeroWind = aeroResults.F_aero;
        mAeroBody = reshape(aeroResults.M_aero, [], 1);

        aoaRad = deg2rad(attitude.aoa);
        aosRad = deg2rad(attitude.aos);
        windToBody = [cos(aosRad) * cos(aoaRad), -sin(aosRad) * cos(aoaRad), -sin(aoaRad); ...
            sin(aosRad), cos(aosRad), 0; ...
            cos(aosRad) * sin(aoaRad), -sin(aosRad) * sin(aoaRad), cos(aoaRad)];

        forceAeroBody = windToBody * forceAeroWind;
        forceAeroEci = rEciToBody' * forceAeroBody;
        aAeroEci = forceAeroEci / config.mass;

        if isfield(aeroResults, 'F_solar')
            forceSolarBody = windToBody * aeroResults.F_solar;
            aAeroEci = aAeroEci + (rEciToBody' * forceSolarBody) / config.mass;
        end

        if isfield(aeroResults, 'M_solar')
            mAeroBody = mAeroBody + reshape(aeroResults.M_solar, [], 1);
        end
    catch
        aAeroEci = [0; 0; 0];
        mAeroBody = [0; 0; 0];
        aeroResults = struct();
    end
end
