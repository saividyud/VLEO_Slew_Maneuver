%Testing the sat template
Xi = [0, 0, 6678e3, -7789, 0, 0,0.7860666, 0.1675188, 0.5709415, 0.1675188,0,0,0]';
Xr = [1 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
y = zeros(13,1);
options = odeset('MaxStep',1);
[t,y] = ode45(@(t,X) Sat_template2_linear(t,X,Xr,y),linspace(0,500,500),Xi,options);

%quaternions to euler angles
q0 = y(:,7);
q1 = y(:,8);
q2 = y(:,9);
q3 = y(:,10);
wx = y(:,11);
wy = y(:,12);
wz = y(:,13);
%coversion to roll pitch yaw

%plotting
%orientation euler angles
figure(1)
phi = atan2(2*(q0.*q1 + q2.*q3),1-2*(q1.^2+q2.^2));
theta = asin(2*(q0.*q2-q1.*q3));
psi = atan2(2*(q0.*q3 + q1.*q2),1-2*(q2.^2 + q3.^2));
plot(t,phi,t,theta,t,psi);
legend('Phi','Theta','Psi')

%angular velocities
figure(2)
plot(t,wx,t,wy,t,wz);
legend('omega x','omega y','omega z')

%Position
figure(3)
plot3(y(:,1),y(:,2),y(:,3))
