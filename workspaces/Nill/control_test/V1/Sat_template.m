% Sat_template calculates the time derivative of the satellite's state vector.

% This function implements 6-DOF satellite dynamics with proportional attitude control.
% All vector quantities in the input and output state vectors (position, velocity, 
% angular velocity) are expressed in the Earth-Centered Inertial (ECI) reference frame.

% State Vector X (13x1):
% - X(1:3): Position vector in ECI frame [m].
% - X(4:6): Velocity vector in ECI frame [m/s].
% - X(7:10): Quaternion [qx, qy, qz, qw] for the rotation from ECI to the body frame.
% - X(11:13): Angular velocity of the body relative to ECI, expressed in ECI frame [rad/s].

function Xd = Sat_template(t, X, target_pointing_eci, params, varargin)
% Input parser
p = inputParser;
addParameter(p, 'useJ2', true, @islogical);
addParameter(p, 'useAtmDrag', true, @islogical);
addParameter(p, 'useControl', true, @islogical);
addParameter(p, 'useSRP', false, @islogical);
parse(p, varargin{:});
flags = p.Results;

% Initialize output
Xd = zeros(13, 1);

% Extract parameters
mu = params.mu;
mass = params.mass;
I_CB = params.I_CB; % Moment of inertia in the body frame
Kp_att = params.Kp_att; % Proportional gain for attitude control
radius = params.radius; % Satellite radius for thrust calculation

% Extract state components
r_eci = X(1:3);
v_eci = X(4:6);
q_eci_to_body = X(7:10);
omega_eci = X(11:13);

% Normalize quaternion to prevent numerical drift
q_eci_to_body = q_eci_to_body / norm(q_eci_to_body);

% === Frame Conversion ===
% Reconstruct the rotation matrix from the quaternion
R_Body_to_ECI = quaternion_to_dcm(q_eci_to_body);
R_ECI_to_Body = R_Body_to_ECI';

% Convert angular velocity from ECI to the body frame for dynamics calculation
omega_body = R_ECI_to_Body * omega_eci;

% === Translational Dynamics (in ECI frame) ===
Xd(1:3) = v_eci;
r_norm = norm(r_eci);
a_2body = -mu * r_eci / r_norm^3;
a_total = a_2body; % Initialize with gravity

Xd(4:6) = a_total;

% === Attitude Dynamics ===
% 1. Quaternion Kinematics (requires body-frame angular velocity)
Omega_matrix = [ 0, omega_body(3), -omega_body(2), omega_body(1);
                -omega_body(3), 0, omega_body(1), omega_body(2);
                 omega_body(2), -omega_body(1), 0, omega_body(3);
                -omega_body(1), -omega_body(2), -omega_body(3), 0 ];
Xd(7:10) = 0.5 * Omega_matrix * q_eci_to_body;

% 2. Rotational Dynamics (solved in body frame)
tau_body = zeros(3, 1);

if flags.useControl
    % Current camera pointing direction in ECI (body +Z axis)
    camera_body = [0; 0; 1];
    current_pointing_eci = R_Body_to_ECI * camera_body;
    current_pointing_eci = current_pointing_eci / norm(current_pointing_eci);

    % Target pointing direction (already in ECI, normalized)
    target_pointing_eci = target_pointing_eci / norm(target_pointing_eci);

    % Compute pointing error vector (in ECI frame)
    pointing_error_eci = cross(current_pointing_eci, target_pointing_eci);

    % Convert pointing error to body frame
    pointing_error_body = R_ECI_to_Body * pointing_error_eci;

    % Proportional control law for torque (in body frame)
    % Torque is proportional to the pointing error
    tau_body = Kp_att * pointing_error_body;

    % Store torque in params for logging (convert to ECI for consistency)
    params.tau_eci = R_Body_to_ECI * tau_body;
end

% Euler's equation of motion (in body frame)
alpha_body = I_CB \ (tau_body - cross(omega_body, I_CB * omega_body));

% Convert angular acceleration from body frame back to ECI frame for the state output
alpha_eci = R_Body_to_ECI * alpha_body;
Xd(11:13) = alpha_eci;

end

% --- Helper Functions ---
function R = quaternion_to_dcm(q)
% Converts a scalar-last quaternion [qx, qy, qz, qw] to a DCM (Body to ECI)
qx = q(1); qy = q(2); qz = q(3); qw = q(4);
R = [1 - 2*(qy^2 + qz^2),   2*(qx*qy - qw*qz),   2*(qx*qz + qw*qy);
     2*(qx*qy + qw*qz),   1 - 2*(qx^2 + qz^2),   2*(qy*qz - qw*qx);
     2*(qx*qz - qw*qy),   2*(qy*qz + qw*qx),   1 - 2*(qx^2 + qy^2)];
end
