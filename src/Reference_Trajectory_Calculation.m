t = 0:1:3600*.25;
%initial location of reference point on the ground
Ref_point = [6378e3*cos(15/57.3) 6378e3*sin(15/57.3) 0];

%angular velocity of Earth
om = -7.29e-5;

%solve for reference point as Earth rotates
    for i = 2:length(t)
        r_ECEF_ECI = [cos(om .*t(i)) -sin(om .*t(i)) 0 ; sin(om .*t(i)) cos(om.*t(i)) 0; 0 0 1];
        Ref_point(end+1,:) = r_ECEF_ECI * [6378e3*cos(15/57.3) 6378e3*sin(15/57.3) 0]';
    end

%initial conditions of satellite
Xi = [6628140; 0; 0; 0; 7754.84333577692; 0; 0.5; -0.5; -0.5; 0.5; 0; 0; 0]';

%integrate to find satellite position
opts = odeset('RelTol', 1e-12,'AbsTol', 1e-12);
[t, X] = ode45(@Sat_template, t, Xi, opts);

%plotting
figure(1)
hold on;
plot3(X(:,1),X(:,2),X(:,3),'o');
plot3(Ref_point(:,1),Ref_point(:,2),Ref_point(:,3),'*');
axis equal;
hold off;

%calculate pointing vector
v_point = X(:,1:3) - Ref_point(:,1:3);

%plot pointing vector
figure(2)
plot3(v_point(:,1),v_point(:,2),v_point(:,3),'*');
axis equal;