function [a_J23] = J2(X, params)
% J2 computes perturbation acceleration due to Earth's oblateness
% (Version without warnings for repeated calls)

    % Input validation
    if length(X) < 3
        a_J23 = zeros(3, 1);
        return;
    end

    % Define constants
    mu = params.mu;

    if isfield(params, 'R_e')
        R_e = params.R_e;
    else
        R_e = 6378137;
    end

    if isfield(params, 'J2')
        J2_coeff = params.J2;
    else
        J2_coeff = 1.08263e-3;
    end

    if isfield(params, 'J3')
        J3_coeff = params.J3;
    else
        J3_coeff = -2.53e-6;
    end

    % Extract position
    r_vec = X(1:3);
    x = X(1);
    y = X(2);
    z = X(3);

    r = norm(r_vec);

    % If below surface, return zero (integration will stop)
    if r < R_e
        a_J23 = zeros(3, 1);
        return;
    end

    % Calculate spherical coordinates
    phi = atan2(z, sqrt(x^2 + y^2));
    lambda = atan2(y, x);

    % Rotation matrix
    R_IS = [
        cos(phi) * cos(lambda), -sin(phi) * cos(lambda), -sin(lambda);
        cos(phi) * sin(lambda), -sin(phi) * sin(lambda),  cos(lambda);
        sin(phi),                cos(phi),                0
    ];

    % J2 acceleration
    factor_J2 = (3 * mu * J2_coeff * R_e^2) / (2 * r^4);
    a_J2_sph = factor_J2 * [
        (1 - 3 * sin(phi)^2);
        sin(2 * phi);
        0
    ];

    % J3 acceleration
    factor_J3 = (5 * mu * J3_coeff * R_e^3) / (2 * r^5);
    a_J3_sph = factor_J3 * [
        3 * sin(phi) * (1 - (7/5) * sin(phi)^2);
        (6/5) * cos(phi) * (1 - (7/3) * sin(phi)^2);
        0
    ];

    % Total and transform
    a_total_sph = a_J2_sph + a_J3_sph;
    a_J23 = R_IS * a_total_sph;

end