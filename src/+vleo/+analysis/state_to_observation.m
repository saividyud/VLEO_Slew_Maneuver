function obs = state_to_observation(X, ~)
    if size(X, 2) > 1
        X = X';
    end

    qScalarFirst = reshape(X(7:10), 1, []);
    qScalarFirst = qScalarFirst / norm(qScalarFirst);
    cameraBody = [0; 0; 1];
    rBodyToEci = quat2dcm(qScalarFirst)';
    pointingEci = rBodyToEci * cameraBody;
    pointingEci = pointingEci / norm(pointingEci);

    decRad = asin(pointingEci(3));
    raRad = atan2(pointingEci(2), pointingEci(1));
    if raRad < 0
        raRad = raRad + 2 * pi;
    end

    obs.dec = rad2deg(decRad);
    obs.ra = rad2deg(raRad);
    obs.dec_rad = decRad;
    obs.ra_rad = raRad;
    obs.pointing_eci = pointingEci;
end
