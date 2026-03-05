clc;clear; close all;

t = 0:1:3600*.25;
%initial location of reference point on the ground
long = 30;
lat = 45;
Ref_point = [6378e3*cos(long/57.3)*cos(lat/57.3) 6378e3*sin(long/57.3)*cos(lat/57.3) 6378e3*sin(lat/57.3)];

%angular velocity of Earth
om = 7.29e-5;

%solve for reference point as Earth rotates
    for i = 2:length(t)
        r_ECEF_ECI = [cos(om .*t(i)) -sin(om .*t(i)) 0 ; sin(om .*t(i)) cos(om.*t(i)) 0; 0 0 1];
        Ref_point(end+1,:) = r_ECEF_ECI * [6378e3*cos(long/57.3)*cos(lat/57.3) 6378e3*sin(long/57.3)*cos(lat/57.3) 6378e3*sin(lat/57.3)]';
    end

%initial conditions of satellite
earthRadius = 6378.14e3; % Earth equatorial radius [m]
a = 250e3 + earthRadius; % 250 km above Earth semimajor axis
e = 0; % Eccentricity
i = 45; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]
Xi = [r_i(1);r_i(2);r_i(3);v_i(1);v_i(2);v_i(3); 0.5; -0.5; -0.5; 0.5; 0; 0; 0]';

%integrate to find satellite position
[t, X] = ode45(@Sat_template, t, Xi);

%plotting
figure(1)
plot3(X(:,1),X(:,2),X(:,3),'o');
hold on;
plot3(Ref_point(:,1),Ref_point(:,2),Ref_point(:,3),'*');
E = wgs84Ellipsoid;
[x,y,z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x, y, z, FaceAlpha="texturemap", FaceColor="texturemap", EdgeAlpha="texturemap"); % Plot the Earth
legend('Satellite','Reference Location','Earth')
title('Location of Satellite and Reference Location With Rotating Earth')
axis equal;
hold off;

%calculate pointing vector
v_point = X(:,1:3) - Ref_point(:,1:3);
for i = 1:length(v_point)
    v_point(i,1:3) = v_point(i,1:3)/norm(v_point(i,1:3));
end

%plot pointing vector
figure(2)
plot3(v_point(:,1),v_point(:,2),v_point(:,3),'*');
hold on;
title('Pointing Vector From Satellite to Reference Point in ECI');
plot3(0,0,0,'o')

q = zeros(4,length(v_point));
ang = zeros(3,length(v_point));
for i = 1:length(v_point)
    v_1 = v_point(i,1:3);
    up  = [0 0 1];
    v_2 = cross(up, v_1);
    v_2 = v_2 / norm(v_2);
    v_3 = cross(v_1, v_2);
    v_3 = v_3 / norm(v_3);
    DCM = [v_1', v_2', v_3'];
    q(:,i) = QfromDCM(DCM);
    [ang(1,i), ang(2,i), ang(3,i)] = quat2angle(q(:,i)');
end

figure(3)
ang = ang.*57.3;
plot(t,ang(1,:),'o',t,ang(2,:),'o',t,ang(3,:),'o');
legend('Roll','Pitch','Yaw');
title('Euler Angles of Satellite When Pointing to Reference Location')

figure(4)
q= q.*57.3;
plot(t,q(1,:),'o',t,q(2,:),'o',t,q(3,:),'o',t,q(4,:),'o');
legend('q0','q1','q2','q3');
title('Euler Angles of Satellite When Pointing to Reference Location')



