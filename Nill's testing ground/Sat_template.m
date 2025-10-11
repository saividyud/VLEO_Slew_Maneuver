function Xd = Sat_template(t, X, u, params, varargin)
% Sat_template calculates the time rate of change of the state X at a time t.
%
% This function implements 6-DOF satellite dynamics including orbital
% mechanics and attitude dynamics in the ECI (Earth-Centered Inertial) frame.
%
% Parameters
% ----------
% t : float
%     Time [s]
% X : 13x1 vector
%     State of the system with the following attributes:
%     - Position in ECI (1:3) [m]
%     - Velocity in ECI (4:6) [m/s]
%     - Quaternion (7:10) [scalar-last convention: qx, qy, qz, qw]
%     - Angular velocity in body frame (11:13) [rad/s]
% u : 6x1 vector
%     Control inputs:
%     - Force in ECI frame (1:3) [N]
%     - Torque in body frame (4:6) [Nm]
% params : struct
%     Physical parameters with fields:
%     - mu: Gravitational parameter [m^3/s^2]
%     - R_e: Earth radius [m]
%     - mass: Satellite mass [kg]
%     - I_CB: 3x3 moment of inertia tensor [kg*m^2]
% varargin : optional name-value pairs
%     - 'useJ2': boolean (default: true) - Include J2 perturbations
%     - 'useAtmDrag': boolean (default: true) - Include atmospheric drag
%     - 'useControl': boolean (default: true) - Include control forces
%     - 'useSRP': boolean (default: false) - Include solar radiation pressure
%
% Returns
% -------
% Xd : 13x1 vector
%     Time derivative of state vector
%
% Example
% -------
%   % Use all forces (default)
%   Xd = Sat_template(t, X, u, params);
%
%   % Disable J2 and atmospheric drag
%   Xd = Sat_template(t, X, u, params, 'useJ2', false, 'useAtmDrag', false);

    % Parse optional inputs with defaults
    p = inputParser;
    addParameter(p, 'useJ2', true, @islogical);
    addParameter(p, 'useAtmDrag', true, @islogical);
    addParameter(p, 'useControl', true, @islogical);
    addParameter(p, 'useSRP', false, @islogical);
    parse(p, varargin{:});

    flags = p.Results;

    % Input validation
    assert(length(X) == 13, 'State vector must be 13x1');
    assert(length(u) == 6, 'Control input must be 6x1');
    assert(isfield(params, 'mass'), 'params must contain mass');

    % Ensure X is column vector
    if size(X, 2) > 1
        X = X';
    end

    % Initialize output
    Xd = zeros(13, 1);

    % Extract physical parameters
    mu = 3.986004418e14;
    params.J2 = 1.08263e-3;
    params.R_e = 6378137;
    mass = params.mass;

    % Get moment of inertia or use default
    if isfield(params, 'I_CB')
        I_CB = params.I_CB;
    else
        I_CB = mass * eye(3);  % Default to unit sphere
    end

    % Extract position and velocity
    r_vec = X(1:3);
    v_vec = X(4:6);
    r = norm(r_vec);

    % Validate physical state (optional check)
    if isfield(params, 'R_e')
        if r < params.R_e
            warning('Satellite position below Earth surface: r = %.1f km', r/1e3);
        end
    end

    % Extract quaternion and normalize
    q = X(7:10);
    q_norm = norm(q);
    if abs(q_norm - 1.0) > 1e-6
        q = q / q_norm;
    end

    % Extract angular velocity
    omega = X(11:13);

    % === Translational Dynamics ===

    % Position derivative
    Xd(1:3) = v_vec;

    % Two-body acceleration (always included)
    a_2body = -mu * r_vec / r^3;

    % Initialize total acceleration
    a_total = a_2body;

    % J2 perturbation
    if flags.useJ2
        try
            a_J2 = J2(X, params);
            a_total = a_total + a_J2;
        catch
            % If J2 function not available, skip
        end
    end

    % Atmospheric drag
    if flags.useAtmDrag
        try
            a_drag = atmospheric_drag(X, params);
            a_total = a_total + a_drag;
        catch
            % If atmospheric_drag function not available, skip
        end
    end

    % Solar radiation pressure (if enabled)
    if flags.useSRP
        try
            a_srp = solar_radiation_pressure(X, t, params);
            a_total = a_total + a_srp;
        catch
            % If SRP function not available, skip
        end
    end

    % Control force acceleration
    if flags.useControl
        a_control = u(1:3) / mass;
        a_total = a_total + a_control;
    end

    % Total acceleration
    Xd(4:6) = a_total;

    % === Attitude Dynamics ===

    % Quaternion kinematics (scalar-last convention)
    % dq/dt = 0.5 * Omega(omega) * q
    Omega = [    0,      omega(3), -omega(2),  omega(1);
            -omega(3),     0,      omega(1),  omega(2);
             omega(2), -omega(1),     0,      omega(3);
            -omega(1), -omega(2), -omega(3),     0    ];

    Xd(7:10) = 0.5 * Omega * q;

    % Rotational dynamics: I*omega_dot + omega x (I*omega) = tau
    if flags.useControl
        tau_control = u(4:6);
    else
        tau_control = zeros(3, 1);
    end

    I_omega = I_CB * omega;
    omega_cross_I_omega = cross(omega, I_omega);

    Xd(11:13) = I_CB \ (tau_control - omega_cross_I_omega);

end
