clear
clc
close all

% Constants
mu = 3.986e14;

% Defining sample state
X = zeros(13, 1);
X(1) = earthRadius + 350e3; % [m]

r_vec = X(1:3);
r = norm(r_vec);

% Acceleration of gravity
a_grav = -mu*r_vec/r^3

% J2 acceleration
a_J23 = J2(X)