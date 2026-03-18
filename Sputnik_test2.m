clc;clear;close all;
%Testing the sat template
X = [0, 0, 6678e3, -7789, 0, 0, 1, 0.0523491, 0.0473595, 0.0523491,0,0,0,0,0,0,0]';

%initial quaternions are the starting location of the tracking
q_i = vleo.dynamics.Ref(0, X);

%Final quaternions are the reference for linearizing around
q_f = vleo.dynamics.Ref(1000,[-6178127.10093010	0	2702313.19215046]');

%for the initial quaternions of the simulation, need to subtract the
%initial quaternions and the final quaternions since the state that is
%linearized around becomes the origin.
Xi = [0, 0, 6678e3, -7789, 0, 0, q_i' - q_f' ,0,0,0,0,0,0,0]';
Xr = [q_f',0,0,0]';
[t,y] = ode45(@(t,X) vleo.dynamics.sat_dynamics_linearized_augmented(t,X,Xr) , [0 1000],Xi);

%need to add the reference for linearization back to the quaternions
q1 = y(:,8)+q_f(2);
q2 = y(:,9)+q_f(3);
q3 = y(:,10) +q_f(4);

%solve for q0 from q1 q2 and q3
q0 = sqrt(1 -q1.^2 - q2.^2 - q3.^2);
wx = y(:,11);
wy = y(:,12);
wz = y(:,13);

%Plotting the reference
figure(1)
q_ref = zeros(4,length(t));
for i = 1:length(t)
    q_ref(:,i) = Ref(t(i),y(i,:)');
end
plot(t,q_ref(1,:),t,q_ref(2,:),t,q_ref(3,:),t,q_ref(4,:));

%Plotting orientation in quatenions
figure(2)
plot(t,q0,t,q1,t,q2,t,q3);
legend('q0','q1','q2','q3')

%error of quaternions
e0 = q_ref(1,:) - q0';
e1 = q_ref(2,:) - q1';
e2 = q_ref(3,:) - q2';
e3 = q_ref(4,:) - q3';
figure(3)
plot(t,e0,t,e1,t,e2,t,e3);
legend('q0','q1','q2','q3')

%error of  euler angles
phi = atan2(2*(q0.*q1 + q2.*q3),1-2*(q1.^2+q2.^2));
theta = asin(2*(q0.*q2-q1.*q3));
psi = atan2(2*(q0.*q3 + q1.*q2),1-2*(q2.^2 + q3.^2));

q1r = q_ref(2,:)';
q2r = q_ref(3,:)';
q3r = q_ref(4,:)';
q0r = q_ref(1,:)';

phir = atan2(2*(q0r.*q1r + q2r.*q3r),1-2*(q1r.^2+q2r.^2));
thetar = asin(2*(q0r.*q2r-q1r.*q3r));
psir = atan2(2*(q0r.*q3r + q1r.*q2r),1-2*(q2r.^2 + q3r.^2));

phie = phi - phir;
thetae = theta - thetar;
psie = psi - psir;
figure(4)
plot(t,phie,t,thetae,t,psie);

% %angular velocities
% figure(2)
% plot(t,wx,t,wy,t,wz);
% legend('omega x','omega y','omega z')
% 
% %Position
% figure(3)
% plot3(y(:,1),y(:,2),y(:,3))


