function setup_volcano_scenario()
% SETUP_VOLCANO_SCENARIO
% Initialize scenario parameters from control_test2
% Defines: X0, params, r_volcano_0, orbital_period

global X0 params r_volcano_0 orbital_period

% Parameters
params.mu = 3.986004418e14;
params.R_e = 6378137;
params.mass = 83.6;
params.radius = 0.58/2;
params.J2 = 1.08263e-3;
params.I_CB = 2/5 * params.mass * params.radius^2 * eye(3);
params.Kp_att = 10;
params.omega_earth = 7.2921159e-5;

% Orbital parameters
altitude = 200e3;
r_orbit = params.R_e + altitude;
v_circ = sqrt(params.mu / r_orbit);

% Random orbital state
rng(42);  % Fixed seed for reproducibility
inclination = rand() * pi;
RAAN = rand() * 2 * pi;
true_anomaly = rand() * 2 * pi;

R_RAAN = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
R_inc = [1, 0, 0; 0, cos(inclination), -sin(inclination); 0, sin(inclination), cos(inclination)];
R_ta = [cos(true_anomaly), -sin(true_anomaly), 0; sin(true_anomaly), cos(true_anomaly), 0; 0, 0, 1];

R_perifocal_to_ECI = R_RAAN * R_inc * R_ta;

r_perifocal = [r_orbit; 0; 0];
v_perifocal = [0; v_circ; 0];

r_eci = R_perifocal_to_ECI * r_perifocal;
v_eci = R_perifocal_to_ECI * v_perifocal;

% Initial attitude - nadir pointing
orbit_normal_eci = cross(r_eci, v_eci) / norm(cross(r_eci, v_eci));
z_body_eci = -r_eci / norm(r_eci);
x_body_eci = v_eci / norm(v_eci);
y_body_eci = cross(z_body_eci, x_body_eci);

R_Body_to_ECI = [x_body_eci, y_body_eci, z_body_eci];
q0 = dcm_to_quaternion(R_Body_to_ECI);

omega_orbit_mag = v_circ / r_orbit;
omega_eci = omega_orbit_mag * orbit_normal_eci;

X0 = [r_eci; v_eci; q0; omega_eci];

% Orbital period
orbital_period = 2 * pi * sqrt(r_orbit^3 / params.mu);

% Volcano position at t=0
r_sat_0 = X0(1:3);
r_projected = r_sat_0 / norm(r_sat_0) * params.R_e;
theta_max = acos(params.R_e / norm(r_sat_0));

theta = 0.3 * theta_max;  % Fixed for reproducibility
phi = rand() * 2 * pi;

axis1 = v_eci / norm(v_eci);
r_temp = r_projected * cos(theta) + ...
         cross(axis1, r_projected) * sin(theta) + ...
         axis1 * dot(axis1, r_projected) * (1 - cos(theta));

axis2 = r_projected / norm(r_projected);
r_volcano_0 = r_temp * cos(phi) + ...
              cross(axis2, r_temp) * sin(phi) + ...
              axis2 * dot(axis2, r_temp) * (1 - cos(phi));

r_volcano_0 = r_volcano_0 / norm(r_volcano_0) * params.R_e;

fprintf('Initial scenario setup complete\\n');
fprintf('Orbital period: %.2f minutes\\n', orbital_period/60);

end

function q = dcm_to_quaternion(R)
% Convert DCM to quaternion [qx, qy, qz, qw]
trace_R = trace(R);

if trace_R > 0
    s = 0.5 / sqrt(trace_R + 1.0);
    qw = 0.25 / s;
    qx = (R(3,2) - R(2,3)) * s;
    qy = (R(1,3) - R(3,1)) * s;
    qz = (R(2,1) - R(1,2)) * s;
elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
    s = 2.0 * sqrt(1.0 + R(1,1) - R(2,2) - R(3,3));
    qw = (R(3,2) - R(2,3)) / s;
    qx = 0.25 * s;
    qy = (R(1,2) + R(2,1)) / s;
    qz = (R(1,3) + R(3,1)) / s;
elseif R(2,2) > R(3,3)
    s = 2.0 * sqrt(1.0 + R(2,2) - R(1,1) - R(3,3));
    qw = (R(1,3) - R(3,1)) / s;
    qx = (R(1,2) + R(2,1)) / s;
    qy = 0.25 * s;
    qz = (R(2,3) + R(3,2)) / s;
else
    s = 2.0 * sqrt(1.0 + R(3,3) - R(1,1) - R(2,2));
    qw = (R(2,1) - R(1,2)) / s;
    qx = (R(1,3) + R(3,1)) / s;
    qy = (R(2,3) + R(3,2)) / s;
    qz = 0.25 * s;
end

q = [qx; qy; qz; qw];
q = q / norm(q);
end