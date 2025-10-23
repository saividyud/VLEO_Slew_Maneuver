clear
clc
close all

%% Initial state (in Keplerian orbital elements)
a = 350e3 + earthRadius; % 350 km above Earth semimajor axis
e = 0; % Eccentricity
i = 0; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]

beta_i = [0, 0, 0, 0]; % Initial quaternion
omega_i = [0, 0, 0]; % Initial angular rate

X_i = [r_i, v_i, beta_i, omega_i]';

%% Simulating
% Simualation bounds
t0 = 0;
t_span = 1 * 24 * 3600; % 1 day
dt = 1;

ts = t0 : dt : t_span;

opts = odeset('RelTol', 1e-12,'AbsTol', 1e-12);
[t, X] = ode45(@Sat_template, ts, X_i, opts);

% Extract position and velocity from the state vector
rs = X(:, 1:3);
vs = X(:, 4:6);
betas = X(:, 6:10);
omegas = X(:, 10:13);

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

xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Satellite Trajectory');

view = 2 * earthRadius;
xlim([-view, view])
ylim([-view, view])
zlim([-view, view])

grid on;