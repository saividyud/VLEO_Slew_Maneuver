function [a_srp] = solar_radiation_pressure(X, t, params)
% solar_radiation_pressure calculates SRP acceleration
%
% Simplified model assuming constant solar flux and cylindrical shadow
%
% Parameters
% ----------
% X : 13x1 vector
%     Current state
% t : float
%     Current time [s]
% params : struct
%     Parameters with fields:
%     - mass: Satellite mass [kg]
%     - Cr: Reflectivity coefficient (default: 1.3)
%     - A_srp: Cross-sectional area [m^2] (default: 1.0)
%
% Returns
% -------
% a_srp : 3x1 vector
%     SRP acceleration [m/s^2]

    % Solar flux constant at 1 AU
    P_sun = 4.56e-6;  % [N/m^2]

    % Get SRP parameters or use defaults
    if isfield(params, 'Cr')
        Cr = params.Cr;
    else
        Cr = 1.3;
    end

    if isfield(params, 'A_srp')
        A = params.A_srp;
    else
        A = 1.0;  % [m^2]
    end

    % Simplified: assume sun direction is +X
    % In reality, should calculate from time
    sun_dir = [1; 0; 0];

    % Check if satellite is in Earth's shadow (simplified cylindrical model)
    r_vec = X(1:3);
    if isfield(params, 'R_e')
        R_e = params.R_e;
    else
        R_e = 6378137;
    end

    % Simple shadow check
    if dot(r_vec, sun_dir) < 0 && norm(r_vec - dot(r_vec, sun_dir)*sun_dir) < R_e
        % In shadow
        a_srp = zeros(3, 1);
    else
        % In sunlight
        a_srp = Cr * A / params.mass * P_sun * sun_dir;
    end

end
