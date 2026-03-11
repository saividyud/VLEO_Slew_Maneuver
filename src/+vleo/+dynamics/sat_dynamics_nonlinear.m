function Xd = sat_dynamics_nonlinear(t, X)
    vehicle = vleo.dynamics.default_vehicle_config();
    mu = 3.986e14;
    rNorm = norm(X(1:3));

    aeroConfig = struct();
    aeroConfig.objFilePath = vehicle.objFilePath;
    aeroConfig.mass = vehicle.mass;
    aeroConfig.year = 2002;
    aeroConfig.dayOfYear = 106;
    aeroConfig.secondsOfDay = 0;
    aeroConfig.options = default_aero_options();
    [aAeroEci, mAeroBody] = vleo.dynamics.compute_aero_effects(t, X, aeroConfig);

    Xd = zeros(13, 1);
    gravityModelMu = 3.986e14;
    [gx, gy, gz] = gravityzonal(X(1:3)', 'Custom', 6378.14e3, gravityModelMu, [1082e-6, -2.53e-6, 0], 'None');
    aJ2 = [gx; gy; gz] + gravityModelMu * X(1:3) / rNorm^3;

    Xd(1:3) = X(4:6);
    Xd(4:6) = -mu * X(1:3) / rNorm^3 + aJ2 + aAeroEci;

    bMatrix = [X(7), -X(8), -X(9), -X(10); ...
        X(8), X(7), -X(10), X(9); ...
        X(9), X(10), X(7), -X(8); ...
        X(10), -X(9), X(8), X(7)];
    Xd(7:10) = 0.5 * bMatrix * [0; X(11); X(12); X(13)];

    tauControl = vleo.control.attitude_pd_controller(t, X);
    tauTotal = tauControl + mAeroBody;
    wx = [0, -X(13), X(12); X(13), 0, -X(11); -X(12), X(11), 0];
    Xd(11:13) = vehicle.inertiaBody \ (tauTotal - wx * vehicle.inertiaBody * X(11:13));
end

function options = default_aero_options()
    options.f107Average = 150;
    options.f107Daily = 150;
    options.magneticIndex = ones(1, 7) * 3;
    options.anomalousOxygen = 0;
    options.gsi_model = 'cook';
    options.alpha = 1;
    options.Tw = 300;
    options.enableShadow = 1;
    options.enableSolar = 1;
    options.sol_cR = 0.15;
    options.sol_cD = 0.25;
end
