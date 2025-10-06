function [F_J23] = J2(X)
% J2 takes in the state vector of a satellite and outputs a force vector
% due to the blateness and nonuniformity of the Earth.
% 
% Arguments
% ---------
% X : 13x1 vector
%   Current state of the satellite
%
% Returns
% -------
% F_J23 : 3x1 vector
%   Perturbation force due to J2 effect and a higher order effect: J3

    % Defining constants
    mu = 3.986e14; % Gravitational parameter of Earth [m^3 / s^2]
    R_e = earthRadius; % Average radius of the Earth [m]
    J2 = 1082e-6; % J2 zonal harmonic coefficient
    J3 = -2.53e-6; % J3 zonal harmonic coefficient
    
    % Extracting current inertial position of satellite
    r_vec = X(1:3);
    
    x = X(1);
    y = X(2);
    z = X(3);
    
    % Calculating spherical coordinates (r, phi, lambda), or distance,
    % latitude, longitude
    r = norm(r_vec);
    
    phi = atan2(z, sqrt(x^2 + y^2));
    
    lambda = atan2(y, x);
    
    % Defining rotation matrix from spherical coordinate frame to inertial
    % frame    
    R_IS = [
        cos(phi) * cos(lambda), -sin(phi) * cos(lambda), -sin(lambda);
        cos(phi) * sin(lambda), -sin(phi) * sin(lambda), cos(lambda);
        sin(phi)              , cos(phi)               , 0
    ];
    
    % Defining J2 force in spherical frame
    F_J2_coeff = -mu * R_e^2/r^4 * J2;
    F_J2 = F_J2_coeff * [
        3/4 * (3*cos(phi) - 1);
        3/2 * sin(2*phi);
        0
    ];
    
    % Calculating the J3 force in spherical frame
    F_J3_coeff = -mu * R_e^3/r^5 * J3;
    F_J3 = F_J3_coeff * [
        1/2 * (5*sin(3*phi));
        3/8 * (cos(phi) - 5*cos(3*phi));
        0
    ];
    
    F_J23 = R_IS * (F_J2 + F_J3);

end