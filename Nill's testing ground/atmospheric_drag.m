function a_drag = atmospheric_drag(X, params)
% atmospheric_drag calculates atmospheric drag acceleration
%
% This function computes drag acceleration using TLE B* parameter.
% Uses the standard conversion: ballistic_coeff = 12.741621 * B*
%
% Parameters
% ----------
% X : 13x1 vector
%     State vector with position (1:3) [m] and velocity (4:6) [m/s]
% params : struct
%     Required fields:
%     - R_e: Earth radius [m]
%     - bstar: B* drag term from TLE [dimensionless, 1/earth_radii]
%     Optional fields:
%     - H: Scale height [m] (default: 8500 m for LEO)
%     - rho_0: Reference density [kg/m^3] (default: 1.225)
%     - h_0: Reference altitude [m] (default: 0)
%
% Returns
% -------
% a_drag : 3x1 vector
%     Drag acceleration in ECI frame [m/s^2]
%
% Reference:
%   B* = (C_D * rho_0_sgp4 * A) / (2 * m)
%   Ballistic coefficient = (C_D * A / m) = 12.741621 * B* [m^2/kg]
%   This factor comes from 2 * R_e / rho_0_sgp4 conversion

% Extract position and velocity
r_vec = X(1:3);
v_vec = X(4:6);
r = norm(r_vec);

% Check for required parameters
if ~isfield(params, 'R_e')
    error('params must contain R_e (Earth radius)');
end

if ~isfield(params, 'bstar')
    error('params must contain bstar (B* from TLE)');
end

% Get altitude
altitude = r - params.R_e; % [m]

% Return zero drag if altitude is too high (above 1000 km) or B* is zero
if altitude > 1e6 || abs(params.bstar) < 1e-20
    a_drag = zeros(3, 1);
    return;
end

% Atmospheric parameters
if isfield(params, 'H')
    H = params.H;
else
    H = 8500; % Scale height [m]
end

if isfield(params, 'rho_0')
    rho_0 = params.rho_0;
else
    rho_0 = 1.225; % Sea level density [kg/m^3]
end

if isfield(params, 'h_0')
    h_0 = params.h_0;
else
    h_0 = 0; % Reference altitude [m]
end

% Compute atmospheric density using exponential model
rho = rho_0 * exp(-(altitude - h_0) / H); % [kg/m^3]

% Account for Earth's rotation (atmospheric co-rotation)
omega_earth = 7.2921159e-5; % [rad/s]
omega_vec = [0; 0; omega_earth];

% Velocity relative to atmosphere
v_rel = v_vec - cross(omega_vec, r_vec);
v_rel_mag = norm(v_rel);

if v_rel_mag < 1e-6
    a_drag = zeros(3, 1);
    return;
end

% Convert B* to ballistic coefficient (C_D * A / m) [m^2/kg]
% Standard conversion from CNES/ISAE aerospace libraries
ballistic_coeff = 12.741621 * params.bstar; % [m^2/kg]

% Drag acceleration: a = -0.5 * (C_D*A/m) * rho * v_rel^2 * v_hat
a_drag_mag = 0.5 * ballistic_coeff * rho * v_rel_mag^2;
a_drag = -a_drag_mag * (v_rel / v_rel_mag);

end