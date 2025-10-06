clear
clc
close all

% Constants
mu = 3.986e14;

% Defining sample state
X = zeros(13, 1);
X(1) = 350e3; % [m]

r_vec = X(1:3);
r = norm(r_vec);

% Force of gravity
F_grav = -mu*r_vec/r^3

% J2 force
F_J23 = J2(X)