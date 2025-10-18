% state_to_observation calculates the inertial pointing direction of the satellite's camera.
%
% This function takes the satellite's full state vector and computes the
% celestial coordinates (Right Ascension and Declination) corresponding to the
% direction the camera is pointing.
%
% It assumes the camera is rigidly mounted to the satellite's body and is aligned
% with the body's +Z axis, which is defined to be nadir-pointing.
%
% Parameters
% ----------
% X : 13x1 vector
%     State of the system with the following attributes (all vectors in ECI):
%     - Position in ECI (1:3) [m]
%     - Velocity in ECI (4:6) [m/s]
%     - Quaternion (7:10) [scalar-last convention: qx, qy, qz, qw]
%     - Angular velocity in ECI frame (11:13) [rad/s]
% params : struct (optional)
%     Parameters (not used)
%
% Returns
% -------
% obs : struct
%     Observation output with fields:
%     - dec: Declination [degrees]
%     - ra: Right Ascension [degrees]
%     - pointing_eci: 3x1 unit vector of camera pointing direction in ECI frame
%
function obs = state_to_observation(X, params)

% Ensure X is column vector
if size(X, 2) > 1
    X = X';
end

% Extract quaternion
q = X(7:10);
q = q / norm(q);

% --- CORRECTED CAMERA DEFINITION ---
% Define the camera's pointing vector in the body frame.
% The body frame is now defined with the +Z axis pointing to nadir.
% Therefore, a nadir-pointing camera is aligned with the +Z axis.
camera_body = [0; 0; 1];

% Reconstruct the rotation matrix from Body to ECI coordinates.
qx = q(1); qy = q(2); qz = q(3); qw = q(4);
R_Body_to_ECI = [
    1 - 2*(qy^2 + qz^2),   2*(qx*qy - qw*qz),   2*(qx*qz + qw*qy);
    2*(qx*qy + qw*qz),     1 - 2*(qx^2 + qz^2), 2*(qy*qz - qw*qx);
    2*(qx*qz - qw*qy),     2*(qy*qz + qw*qx),   1 - 2*(qx^2 + qy^2)
];

% Transform the camera's body-fixed pointing vector into the ECI frame.
pointing_eci = R_Body_to_ECI * camera_body;
pointing_eci = pointing_eci / norm(pointing_eci);

% Calculate RA and Dec from the ECI pointing vector.
dec_rad = asin(pointing_eci(3));
ra_rad = atan2(pointing_eci(2), pointing_eci(1));

if ra_rad < 0
    ra_rad = ra_rad + 2*pi;
end

% Package output
obs.dec = rad2deg(dec_rad);
obs.ra = rad2deg(ra_rad);
obs.dec_rad = dec_rad;
obs.ra_rad = ra_rad;
obs.pointing_eci = pointing_eci;

end