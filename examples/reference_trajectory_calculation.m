projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

clc;
clear;
close all;

t = 0:1:3600 * 0.25;
longitudeDeg = 30;
latitudeDeg = 45;
refPoint0 = [6378e3 * cosd(longitudeDeg) * cosd(latitudeDeg), ...
    6378e3 * sind(longitudeDeg) * cosd(latitudeDeg), ...
    6378e3 * sind(latitudeDeg)];

earthRate = 7.29e-5;
refPoint = refPoint0;
for idx = 2:length(t)
    rEcefToEci = [cos(earthRate * t(idx)), -sin(earthRate * t(idx)), 0; ...
        sin(earthRate * t(idx)), cos(earthRate * t(idx)), 0; ...
        0, 0, 1];
    refPoint(end + 1, :) = (rEcefToEci * refPoint0')'; %#ok<AGROW>
end

earthRadius = 6378.14e3;
[rEci, vEci] = vleo.analysis.keplerian_to_eci_safe(earthRadius + 250e3, 0, 45, 0, 0, 0, ...
    'GravitationalParameter', 3.986004e14, 'Action', 'None');
X0 = [rEci(:); vEci(:); 0.5; -0.5; -0.5; 0.5; 0; 0; 0];

[tOut, XOut] = ode45(@vleo.dynamics.sat_dynamics_nonlinear, t, X0);

figure(1);
plot3(XOut(:, 1), XOut(:, 2), XOut(:, 3), 'o');
hold on;
plot3(refPoint(:, 1), refPoint(:, 2), refPoint(:, 3), '*');
E = wgs84Ellipsoid;
[x, y, z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x, y, z, 'FaceAlpha', 'texturemap', 'FaceColor', 'texturemap', 'EdgeAlpha', 'texturemap');
legend('Satellite', 'Reference Location', 'Earth');
title('Location of Satellite and Reference Location With Rotating Earth');
axis equal;
hold off;

vPoint = XOut(:, 1:3) - refPoint(:, 1:3);
for idx = 1:length(vPoint)
    vPoint(idx, :) = vPoint(idx, :) / norm(vPoint(idx, :));
end

figure(2);
plot3(vPoint(:, 1), vPoint(:, 2), vPoint(:, 3), '*');
hold on;
plot3(0, 0, 0, 'o');
title('Pointing Vector From Satellite to Reference Point in ECI');

q = zeros(4, length(vPoint));
ang = zeros(3, length(vPoint));
for idx = 1:length(vPoint)
    v1 = vPoint(idx, 1:3);
    up = [0, 0, 1];
    v2 = cross(up, v1);
    v2 = v2 / norm(v2);
    v3 = cross(v1, v2);
    v3 = v3 / norm(v3);
    dcm = [v1', v2', v3'];
    q(:, idx) = dcm2quat(dcm)';
    [ang(1, idx), ang(2, idx), ang(3, idx)] = quat2angle(q(:, idx)');
end

figure(3);
plot(tOut, rad2deg(ang(1, :)), 'o', tOut, rad2deg(ang(2, :)), 'o', tOut, rad2deg(ang(3, :)), 'o');
legend('Roll', 'Pitch', 'Yaw');
title('Euler Angles of Satellite When Pointing to Reference Location');

figure(4);
plot(tOut, q(1, :), 'o', tOut, q(2, :), 'o', tOut, q(3, :), 'o', tOut, q(4, :), 'o');
legend('q0', 'q1', 'q2', 'q3');
title('Quaternion History When Pointing to Reference Location');
