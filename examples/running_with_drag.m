clear
clc
close all

initialize;

%% Initial state (in Keplerian orbital elements)
a = 250e3 + earthRadius; % 250 km above Earth semimajor axis
e = 0; % Eccentricity
i = 0; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]

% Computing intial body frame such that b_1 is in the direction of motion,
% b_3 is pointing towards surface (camera pointing), and b_2 is normal to
% both b_3 and b_1
b_1_i = v_i' / norm(v_i);
b_3_i = -r_i' / norm(r_i);
b_2_i = cross(b_3_i, b_1_i); % Calculate b_2 as the cross product of b_3 and b_1

% Computing initial quaternion from these initial body axes
R_BI_i = [
    b_1_i';
    b_2_i';
    b_3_i'
];

beta_i = QfromDCM(R_BI_i)'; % Initial quaternion

% Zero initial angular spin
omega_i = [0, 0, 0]; % Initial angular rate

X_i = [r_i, v_i, beta_i, omega_i]';

%% Calculating stuff to put into the aerodynamics
shadow = 0;
inparam.alpha = 1; % Accommodation (altitude dependent)
inparam.Tw = 300; % Wall Temperature [K]

solar = 1;
inparam.sol_cR = 0.15; % Specular Reflectivity
inparam.sol_cD = 0.25; % Diffuse Reflectivity

verb = 1;
del = 0;
