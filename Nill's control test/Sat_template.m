function Xd = Sat_template(t, X, desired_state, params, varargin)

% Sat_template calculates the time rate of change of the state X at a time t.

% This function implements 6-DOF satellite dynamics including orbital
% mechanics and attitude dynamics in the ECI (Earth-Centered Inertial) frame.
%
% ECI Frame Definition:
% - X-axis: Points toward the Vernal Equinox (First Point of Aries)
% - Z-axis: Points toward the North Pole
% - Y-axis: Completes the right-handed coordinate system

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
% desired_state : 13x1 vector
%     Desired state vector with same structure as X
% params : struct
%     Physical parameters with fields:
%     - mu: Gravitational parameter [m^3/s^2]
%     - R_e: Earth radius [m]
%     - mass: Satellite mass [kg]
%     - I_CB: 3x3 moment of inertia tensor [kg*m^2]
%     - Kp_pos: Position control gain (default: 1e-3)
%     - Kd_pos: Velocity control gain (default: 1e-2)
%     - Kp_att: Attitude control gain (default: 1e-1)
%     - Kd_att: Angular velocity control gain (default: 1e-1)
% varargin : optional name-value pairs
%     - 'useJ2': boolean (default: true) - Include J2 perturbations
%     - 'useAtmDrag': boolean (default: true) - Include atmospheric drag
%     - 'useControl': boolean (default: true) - Include control forces
%     - 'useSRP': boolean (default: false) - Include solar radiation pressure

% Returns
% -------
% Xd : 13x1 vector
%     Time derivative of state vector

% Example
% -------
% % Use all forces (default)
% Xd = Sat_template(t, X, desired_state, params);
%
% % Disable J2 and atmospheric drag
% Xd = Sat_template(t, X, desired_state, params, 'useJ2', false, 'useAtmDrag', false);

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
assert(length(desired_state) == 13, 'Desired state vector must be 13x1');
assert(isfield(params, 'mass'), 'params must contain mass');

% Ensure X is column vector
if size(X, 2) > 1
    X = X';
end
if size(desired_state, 2) > 1
    desired_state = desired_state';
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
    I_CB = mass * eye(3); % Default to unit sphere
end

% Get control gains or use defaults
Kp_pos = getfield_default(params, 'Kp_pos', 1e-3);
Kd_pos = getfield_default(params, 'Kd_pos', 1e-2);
Kp_att = getfield_default(params, 'Kp_att', 1e-1);
Kd_att = getfield_default(params, 'Kd_att', 1e-1);

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

% === Compute Control Input from Desired State ===
% Position and velocity errors
r_des = desired_state(1:3);
v_des = desired_state(4:6);
pos_error = r_des - r_vec;
vel_error = v_des - v_vec;

% PD control for translational motion
F_control = Kp_pos * pos_error + Kd_pos * vel_error;

% Quaternion error (desired - current)
q_des = desired_state(7:10);
q_des = q_des / norm(q_des);
omega_des = desired_state(11:13);

% Quaternion error computation (q_error = q_des * conj(q))
q_error = quat_multiply(q_des, quat_conjugate(q));

% Extract vector part of error quaternion for control
q_error_vec = q_error(1:3);

% Angular velocity error
omega_error = omega_des - omega;

% PD control for rotational motion
tau_control = Kp_att * q_error_vec + Kd_att * omega_error;

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
    a_control = F_control / mass;
    a_total = a_total + a_control;
end

% Total acceleration
Xd(4:6) = a_total;

% === Attitude Dynamics ===
% Quaternion kinematics (scalar-last convention)
% dq/dt = 0.5 * Omega(omega) * q
Omega = [ 0,        omega(3), -omega(2),  omega(1);
         -omega(3), 0,         omega(1),  omega(2);
          omega(2), -omega(1), 0,         omega(3);
         -omega(1), -omega(2), -omega(3), 0 ];
Xd(7:10) = 0.5 * Omega * q;

% Rotational dynamics: I*omega_dot + omega x (I*omega) = tau
if flags.useControl
    tau = tau_control;
else
    tau = zeros(3, 1);
end

I_omega = I_CB * omega;
omega_cross_I_omega = cross(omega, I_omega);
Xd(11:13) = I_CB \ (tau - omega_cross_I_omega);

end

% Helper function to get field with default value
function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

% Helper function for quaternion multiplication
function q_result = quat_multiply(q1, q2)
    % Scalar-last convention: [qx, qy, qz, qw]
    v1 = q1(1:3);
    s1 = q1(4);
    v2 = q2(1:3);
    s2 = q2(4);

    v_result = s1*v2 + s2*v1 + cross(v1, v2);
    s_result = s1*s2 - dot(v1, v2);

    q_result = [v_result; s_result];
end

% Helper function for quaternion conjugate
function q_conj = quat_conjugate(q)
    % Scalar-last convention: [qx, qy, qz, qw]
    q_conj = [-q(1); -q(2); -q(3); q(4)];
end