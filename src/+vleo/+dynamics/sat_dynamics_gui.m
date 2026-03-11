function Xd = sat_dynamics_gui(t, X, subFigHandle)
    simParams = guidata(subFigHandle);
    vehicle = vleo.dynamics.default_vehicle_config();
    mu = vleo.util.get_env_value(simParams, 'mu', 3.986004418e14);

    aeroOptions.f107Average = vleo.util.get_env_value(simParams, 'F107Average', 150);
    aeroOptions.f107Daily = vleo.util.get_env_value(simParams, 'F107Daily', 150);
    aeroOptions.magneticIndex = vleo.util.get_env_value(simParams, 'magneticIndices', ones(1, 7) * 3);
    aeroOptions.anomalousOxygen = vleo.util.get_env_value(simParams, 'enableAnomalousOxygen', false);
    aeroOptions.gsi_model = vleo.util.normalize_gsi_model(vleo.util.get_env_value(simParams, 'gasSurfaceInteractionModel', 'cook'));
    aeroOptions.alpha = vleo.util.get_env_value(simParams, 'accommodationCoefficient', 1);
    aeroOptions.Tw = vleo.util.get_env_value(simParams, 'wallTemperature', 300);
    aeroOptions.enableShadow = vleo.util.get_env_value(simParams, 'enableShadowAnalysis', true);
    aeroOptions.enableSolar = vleo.util.get_env_value(simParams, 'enableSolarRadiationPressure', true);
    aeroOptions.sol_cR = vleo.util.get_env_value(simParams, 'specularReflectivity', 0.15);
    aeroOptions.sol_cD = vleo.util.get_env_value(simParams, 'diffuseReflectivity', 0.25);

    aAeroEci = [0; 0; 0];
    mAeroBody = [0; 0; 0];
    if vleo.util.get_mode_flag(simParams, 'enableAero', true)
        aeroConfig = struct();
        aeroConfig.objFilePath = vehicle.objFilePath;
        aeroConfig.mass = vehicle.mass;
        aeroConfig.year = vleo.util.get_env_value(simParams, 'year', 2002);
        aeroConfig.dayOfYear = vleo.util.get_env_value(simParams, 'dayOfYear', 1);
        aeroConfig.secondsOfDay = vleo.util.get_env_value(simParams, 'secondsOfDay', 0);
        aeroConfig.options = aeroOptions;
        [aAeroEci, mAeroBody] = vleo.dynamics.compute_aero_effects(t, X, aeroConfig);
    end

    Xd = zeros(13, 1);
    rNorm = norm(X(1:3));
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

    tauControl = vleo.dynamics.evaluate_controller_torque(simParams, t, X);
    tauTotal = tauControl + mAeroBody;
    wx = [0, -X(13), X(12); X(13), 0, -X(11); -X(12), X(11), 0];
    Xd(11:13) = vehicle.inertiaBody \ (tauTotal - wx * vehicle.inertiaBody * X(11:13));
end
