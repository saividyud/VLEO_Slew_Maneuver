clc;clear;

t = 0:1:3600*.25;
%initial location of reference point on the ground
Ref_point = [6378e3*cos(30/57.3) 6378e3*sin(30/57.3) 4000];

%angular velocity of Earth
om = 7.29e-5;

%solve for reference point as Earth rotates
    for i = 2:length(t)
        r_ECEF_ECI = [cos(om .*t(i)) -sin(om .*t(i)) 0 ; sin(om .*t(i)) cos(om.*t(i)) 0; 0 0 1];
        Ref_point(end+1,:) = r_ECEF_ECI * [6378e3*cos(30/57.3) 6378e3*sin(30/57.3) 4000]';
    end

%initial conditions of satellite
Xi = [6378e3+400e3; 0; 0; 0; 7.672e3; 0; 0.5; -0.5; -0.5; 0.5; 0; 0; 0]';

%integrate to find satellite position
[t, X] = ode45(@Sat_template, t, Xi);

%plotting
figure(1)
hold on;
plot3(X(:,1),X(:,2),X(:,3),'o');
plot3(Ref_point(:,1),Ref_point(:,2),Ref_point(:,3),'*');
axis equal;
hold off;

%calculate pointing vector
v_point = X(:,1:3) - Ref_point(:,1:3);
for i = 1:length(v_point)
    v_point(i,1:3) = v_point(i,1:3)/norm(v_point(i,1:3));
end

%plot pointing vector
figure(2)
hold on
plot3(v_point(:,1),v_point(:,2),v_point(:,3),'*');
plot3(0,0,0,'o')

q = zeros(4,length(v_point));
ang = zeros(3,length(v_point));
for i = 1:length(v_point)
    v_1 = v_point(i,1:3);
    v_3 = [0 0 1];
    v_2 = cross(v_3,v_1);
    v_3 = cross(v_1,v_2);
    DCM = [v_1' , v_2', v_3' ];
    q(:,i) = QfromDCM(DCM);
    [ang(1,i) ,ang(2,i),ang(3,i)] = quat2angle(q(:,i)');
end

figure(3)
plot(t,ang(1,:),'o',t,ang(2,:),'o',t,ang(3,:),'o');

