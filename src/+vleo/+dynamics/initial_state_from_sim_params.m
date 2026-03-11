function X0 = initial_state_from_sim_params(simParams)
    a = simParams.initParams.Orbit.semiMajorAxis;
    e = simParams.initParams.Orbit.eccentricity;
    inc = simParams.initParams.Orbit.inclination;
    raan = simParams.initParams.Orbit.RAAN;
    aop = simParams.initParams.Orbit.argPeriapse;
    ta = simParams.initParams.Orbit.trueAnomaly;
    muOrbit = vleo.util.get_env_value(simParams, 'mu', 3.986004418e14);

    [rEci, vEci] = vleo.analysis.keplerian_to_eci_safe( ...
        a, e, inc, raan, aop, ta, ...
        'GravitationalParameter', muOrbit, 'Action', 'None');

    r0 = reshape(rEci, 1, []);
    v0 = reshape(vEci, 1, []);

    b1 = v0' / norm(v0);
    b3 = -r0' / norm(r0);
    b2 = cross(b3, b1);
    b2 = b2 / norm(b2);
    rBodyFromInertialRef = [b1'; b2'; b3'];

    roll = deg2rad(simParams.initParams.Attitude.roll);
    pitch = deg2rad(simParams.initParams.Attitude.pitch);
    yaw = deg2rad(simParams.initParams.Attitude.yaw);

    cr = cos(roll);
    sr = sin(roll);
    cp = cos(pitch);
    sp = sin(pitch);
    cy = cos(yaw);
    sy = sin(yaw);

    rOffset = [cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr; ...
        sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr; ...
        -sp, cp * sr, cp * cr];

    rBodyFromInertial = rOffset * rBodyFromInertialRef;
    q0 = dcm2quat(rBodyFromInertial);
    if q0(1) < 0
        q0 = -q0;
    end

    omega0 = deg2rad([simParams.initParams.Attitude.rollRate, ...
        simParams.initParams.Attitude.pitchRate, ...
        simParams.initParams.Attitude.yawRate]);

    X0 = [r0, v0, q0, omega0]';
end
