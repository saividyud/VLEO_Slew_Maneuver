clear
clc
close all

%% Initial state (in Keplerian orbital elements)
a = 250e3 + earthRadius; % 250 km above Earth semimajor axis
e = 0; % Eccentricity
i = 20; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]

beta_i = [1, 0, 0, 0]; % Initial quaternion
omega_i = [0, 0, 0]; % Initial angular rate

X_i = [r_i, v_i, beta_i, omega_i]';

%% Simulating
% Simualation bounds
t0 = 0;
t_span = 10*60; % 0.5 hour
dt = 1;

ts = t0 : dt : t_span;

opts = odeset('RelTol', 1e-12,'AbsTol', 1e-12);
[t, X] = ode45(@Sat_template, ts, X_i, opts);

% Extract position and velocity from the state vector
rs = X(:, 1:3);
vs = X(:, 4:6);
betas = X(:, 7:10);
omegas = X(:, 11:13);

%% Plotting
% Plotting position
fig = figure(1);

E = wgs84Ellipsoid;
[x,y,z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x, y, z, FaceAlpha="texturemap", FaceColor="texturemap", EdgeAlpha="texturemap"); % Plot the Earth

hold on;
plot3(rs(:, 1), rs(:, 2), rs(:, 3), LineWidth=1, Color='r')

hold off;

axis equal;

view(30, 30)

xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Satellite Trajectory');

bounds = 2 * earthRadius;
xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])

grid on;