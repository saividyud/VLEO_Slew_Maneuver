projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(projectRoot);
setup_project();

X0 = [0, 0, 6678e3, -7789, 0, 0, 0.7860666, 0.1675188, 0.5709415, 0.1675188, 0, 0, 0]';
Xr = [1; 0; 0; 0; 0; 0; 0];
[tOut, XOut] = ode45(@(t, X) vleo.dynamics.sat_dynamics_linearized(t, X, Xr), linspace(0, 500, 500), X0);

q0 = XOut(:, 7);
q1 = XOut(:, 8);
q2 = XOut(:, 9);
q3 = XOut(:, 10);
wx = XOut(:, 11);
wy = XOut(:, 12);
wz = XOut(:, 13);

figure(1);
phi = atan2(2 * (q0 .* q1 + q2 .* q3), 1 - 2 * (q1.^2 + q2.^2));
theta = asin(2 * (q0 .* q2 - q1 .* q3));
psi = atan2(2 * (q0 .* q3 + q1 .* q2), 1 - 2 * (q2.^2 + q3.^2));
plot(tOut, phi, tOut, theta, tOut, psi);
legend('Phi', 'Theta', 'Psi');

figure(2);
plot(tOut, wx, tOut, wy, tOut, wz);
legend('omega x', 'omega y', 'omega z');

figure(3);
plot3(XOut(:, 1), XOut(:, 2), XOut(:, 3));
