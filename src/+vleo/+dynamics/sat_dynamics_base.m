function Xd = sat_dynamics_base(t, X, mu, aeroConfig, tauControl, inertiaBody, useJ2)
    if nargin < 4 || isempty(aeroConfig)
        aeroConfig = [];
    end
    if nargin < 5 || isempty(tauControl)
        tauControl = [0; 0; 0];
    end
    if nargin < 6 || isempty(inertiaBody)
        vehicle = vleo.dynamics.default_vehicle_config();
        inertiaBody = vehicle.inertiaBody;
    end
    if nargin < 7 || isempty(useJ2)
        useJ2 = true;
    end

    constants = vleo.util.constants();
    aAeroEci = [0; 0; 0];
    mAeroBody = [0; 0; 0];
    if ~isempty(aeroConfig)
        [aAeroEci, mAeroBody] = vleo.dynamics.compute_aero_effects(t, X, aeroConfig);
    end

    rEci = X(1:3);
    vEci = X(4:6);
    omegaBody = X(11:13);
    rNorm = norm(X(1:3));
    aJ2 = [0; 0; 0];
    if useJ2
        [gx, gy, gz] = gravityzonal(rEci', 'Custom', constants.R_earth, mu, constants.gravity_zonal_coeffs, 'None');
        aJ2 = [gx; gy; gz] + mu * rEci / rNorm^3;
    end

    Xd = zeros(13, 1);
    Xd(1:3) = vEci;
    Xd(4:6) = -mu * rEci / rNorm^3 + aJ2 + aAeroEci;

    bMatrix = [X(7), -X(8), -X(9), -X(10); ...
        X(8), X(7), -X(10), X(9); ...
        X(9), X(10), X(7), -X(8); ...
        X(10), -X(9), X(8), X(7)];
    Xd(7:10) = 0.5 * bMatrix * [0; omegaBody];

    tauTotal = reshape(tauControl, [], 1) + mAeroBody;
    Xd(11:13) = inertiaBody \ (tauTotal - cross(omegaBody, inertiaBody * omegaBody));
end
