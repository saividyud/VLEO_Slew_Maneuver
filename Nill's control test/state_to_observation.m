function obs = state_to_observation(X, params)

% state_to_observation transforms the satellite state to camera observation output
%
% This function computes the Right Ascension (RA) and Declination (Dec) that
% the camera is pointing, assuming the camera is nadir-pointing in the body frame.
%
% ECI Frame Definition:
% - X-axis: Points toward the Vernal Equinox (First Point of Aries)
% - Z-axis: Points toward the North Pole
% - Y-axis: Completes the right-handed coordinate system
%
% Camera Configuration:
% - In the initial/nominal state, the camera is nadir-pointing (pointing toward Earth center)

% Parameters
% ----------
% X : 13x1 vector
%     State of the system with the following attributes:
%     - Position in ECI (1:3) [m]
%     - Velocity in ECI (4:6) [m/s]
%     - Quaternion (7:10) [scalar-last convention: qx, qy, qz, qw]
%     - Angular velocity in body frame (11:13) [rad/s]
% params : struct (optional)
%     Parameters (currently not used, but included for future extensions)

% Returns
% -------
% obs : struct
%     Observation output with fields:
%     - dec: Declination [degrees] - angle above/below celestial equator (-90 to +90)
%     - ra: Right Ascension [degrees] - angle measured eastward from Vernal Equinox (0 to 360)
%     - dec_rad: Declination [radians]
%     - ra_rad: Right Ascension [radians]
%     - pointing_eci: 3x1 unit vector of camera pointing direction in ECI frame

% Example
% -------
% obs = state_to_observation(X);
% fprintf('Camera pointing: RA = %.2f deg, Dec = %.2f deg\n', obs.ra, obs.dec);

% Ensure X is column vector
if size(X, 2) > 1
    X = X';
end

% Extract position and quaternion
r_vec = X(1:3);
q = X(7:10);

% Normalize quaternion
q = q / norm(q);

% Define nadir direction in body frame (pointing toward Earth)
% In body frame, nadir is typically along -Z axis (assuming standard spacecraft convention)
% This means the camera points in the -Z_body direction
nadir_body = [0; 0; -1];

% Convert quaternion to rotation matrix (ECI to Body)
% For scalar-last convention: q = [qx, qy, qz, qw]
qx = q(1); qy = q(2); qz = q(3); qw = q(4);

% Rotation matrix from ECI to Body
R_ECI_to_Body = [
    1 - 2*(qy^2 + qz^2),     2*(qx*qy - qw*qz),     2*(qx*qz + qw*qy);
    2*(qx*qy + qw*qz),       1 - 2*(qx^2 + qz^2),   2*(qy*qz - qw*qx);
    2*(qx*qz - qw*qy),       2*(qy*qz + qw*qx),     1 - 2*(qx^2 + qy^2)
];

% Rotation matrix from Body to ECI (transpose)
R_Body_to_ECI = R_ECI_to_Body';

% Transform nadir direction from body frame to ECI frame
pointing_eci = R_Body_to_ECI * nadir_body;

% Normalize (should already be unit vector, but ensure numerical stability)
pointing_eci = pointing_eci / norm(pointing_eci);

% Convert pointing vector to Right Ascension and Declination
% Declination: angle above the celestial equator
% Dec = arcsin(z_component)
dec_rad = asin(pointing_eci(3));
dec_deg = rad2deg(dec_rad);

% Right Ascension: angle in the equatorial plane from Vernal Equinox
% RA = atan2(y_component, x_component)
ra_rad = atan2(pointing_eci(2), pointing_eci(1));

% Convert RA to [0, 2*pi) range
if ra_rad < 0
    ra_rad = ra_rad + 2*pi;
end

ra_deg = rad2deg(ra_rad);

% Package output
obs.dec = dec_deg;
obs.ra = ra_deg;
obs.dec_rad = dec_rad;
obs.ra_rad = ra_rad;
obs.pointing_eci = pointing_eci;

end
