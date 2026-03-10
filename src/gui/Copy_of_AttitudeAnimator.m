clear
clc
close all



%% Initial state (in Keplerian orbital elements)
earthRadius = 6378.14e3; % Earth equatorial radius [m]
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

%% Simulating
% Simualation bounds
t0 = 0;
t_span = 5*60; % 0.5 hour
dt = 1;

ts = t0 : dt : t_span;

opts = odeset('RelTol', 1e-12,'AbsTol', 1e-12);
[t, X] = ode45(@Sat_template, ts, X_i, opts);

% Extract position and velocity from the state vector
rs = X(:, 1:3);
vs = X(:, 4:6);
betas = X(:, 7:10);
omegas = X(:, 11:13);

% Extract torques
torques = zeros(length(rs), 3);
for i = 1 : 1 : length(rs)
    torques(i, :) = control_torques(ts(i), X(i, :)')';
end

[yaws, pitches, rolls] = quat2angle(betas, "ZYX");
rolls = rad2deg(rolls);
pitches = rad2deg(pitches);
yaws = rad2deg(yaws);

