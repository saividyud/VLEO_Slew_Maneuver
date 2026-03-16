projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

clc;
clear;
close all;

t = 0:1:3600 * 0.25;
c = vleo.util.constants();
longitudeDeg = 15;
latitudeDeg = 0;
refPoint0 = [c.R_earth * cosd(longitudeDeg) * cosd(latitudeDeg), ...
    c.R_earth * sind(longitudeDeg) * cosd(latitudeDeg), ...
    c.R_earth * sind(latitudeDeg)];

earthRate = c.omega_earth;
refPoint = zeros(numel(t), 3);
refPoint(1, :) = refPoint0;
for idx = 2:numel(t)
    rEcefToEci = [cos(earthRate * t(idx)), -sin(earthRate * t(idx)), 0; ...
        sin(earthRate * t(idx)), cos(earthRate * t(idx)), 0; ...
        0, 0, 1];
    refPoint(idx, :) = (rEcefToEci * refPoint0')';
end

earthRadius = c.R_earth;
[rEci, vEci] = vleo.analysis.keplerian_to_eci_safe(earthRadius + 250e3, 0, 0, 0, 0, 0);
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

vPoint = refPoint(:,1:3)-XOut(:,1:3);

figure(2);
plot3(vPoint(:, 1), vPoint(:, 2), vPoint(:, 3), '*');
hold on;
plot3(0, 0, 0, 'o');
title('Pointing Vector From Satellite to Reference Point in ECI');

q = zeros(4, length(vPoint));
ang = zeros(3, length(vPoint));
for idx = 1:length(vPoint)
    v1 = vPoint(idx, 1:3);
    v1 = v1/norm(v1);
    up = [0, 0, 1];
    v2 = cross(up, v1);
    v2 = v2 / norm(v2);
    v3 = cross(v1, v2);
    v3 = v3 / norm(v3);
    dcm = [v1', v2', v3'];
    q(:, idx) = dcm2quat(dcm)';
    if idx > 131
        q(:, idx) = -q(:,idx);
    end
    [ang(1,idx) , ang(2,idx), ang(3,idx) ] = quat2angle(q(:, idx)');
end

figure(3);
plot(tOut, rad2deg(ang(1, :)), 'o', tOut, rad2deg(ang(2, :)), 'o', tOut, rad2deg(ang(3, :)), 'o');
legend('Roll', 'Pitch', 'Yaw');
title('Euler Angles of Satellite When Pointing to Reference Location');

figure(4);
plot(tOut, q(1, :), 'o', tOut, q(2, :), 'o', tOut, q(3, :), 'o', tOut, q(4, :), 'o');
legend('q0', 'q1', 'q2', 'q3');
title('Quaternion History When Pointing to Reference Location');
