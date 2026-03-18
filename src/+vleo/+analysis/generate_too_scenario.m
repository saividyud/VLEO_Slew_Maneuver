function scenario = generate_too_scenario(params)
    altitude = 200e3;
    rOrbit = params.R_e + altitude;
    vCirc = sqrt(params.mu / rOrbit);

    inclination = rand() * pi;
    raan = rand() * 2 * pi;
    trueAnomaly = rand() * 2 * pi;

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
    rProjected = x0(1:3) / norm(x0(1:3)) * params.R_e;
    thetaMax = acos(params.R_e / norm(x0(1:3)));
    theta = rand() * thetaMax;
    phi = rand() * 2 * pi;

    axis1 = vEci / norm(vEci);
    rTemp = rProjected * cos(theta) + ...
        cross(axis1, rProjected) * sin(theta) + ...
        axis1 * dot(axis1, rProjected) * (1 - cos(theta));

    axis2 = rProjected / norm(rProjected);
    rTOO0 = rTemp * cos(phi) + ...
        cross(axis2, rTemp) * sin(phi) + ...
        axis2 * dot(axis2, rTemp) * (1 - cos(phi));
    rTOO0 = rTOO0 / norm(rTOO0) * params.R_e;

    tooDirection = rTOO0 / norm(rTOO0);
    tooLat = rad2deg(asin(tooDirection(3)));
    tooLon = rad2deg(atan2(tooDirection(2), tooDirection(1)));
    visibilityAtEruption = dot(rTOO0, x0(1:3) - rTOO0);

    zenithAtTOO = rTOO0 / norm(rTOO0);
    tooToSat0 = x0(1:3) - rTOO0;
    tooToSat0Unit = tooToSat0 / norm(tooToSat0);
    cosZenith0 = dot(zenithAtTOO, tooToSat0Unit);
    zenithDistance0 = rad2deg(acos(vleo.util.clamp_scalar(cosZenith0, -1, 1)));

    scenario = struct();
    scenario.altitude = altitude;
    scenario.r_orbit = rOrbit;
    scenario.v_circ = vCirc;
    scenario.r_eci = rEci;
    scenario.v_eci = vEci;
    scenario.z_body_eci = zBodyEci;
    scenario.X0 = x0;
    scenario.orbital_period = orbitalPeriod;
    scenario.r_too_0 = rTOO0;
    scenario.too_lat = tooLat;
    scenario.too_lon = tooLon;
    scenario.theta_max = thetaMax;
    scenario.theta = theta;
    scenario.phi = phi;
    scenario.visibility_at_eruption = visibilityAtEruption;
    scenario.zenith_distance_0 = zenithDistance0;
    scenario.t_start = -orbitalPeriod / 2;
    scenario.t_end = orbitalPeriod / 2;
    scenario.t_eruption = 0;
    scenario.dt = 1;
    scenario.tspan = scenario.t_start:scenario.dt:scenario.t_end;
end
