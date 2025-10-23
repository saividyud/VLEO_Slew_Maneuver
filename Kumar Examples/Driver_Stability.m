clc;
clear all;
close all;

global Ip sig IW Ww;
Ip = diag([350; 500; 400]);
IW = 25;
%%%
% draw ellipsoid
%%%
a = 1.5;
b = 2.5;
c = 5;
maxr = max([a,b,c]);
dtr = pi/180;
[xe ye ze] = ellipsoid(0,0,0,a,b,c, 50);
[xs ys zs] = sphere(50);
surf(xe, ye, ze);
axis([-maxr maxr -maxr maxr -maxr maxr]);

% shading interp;
sig = Ip(3,3)/Ip(1,1);
w0 = 4;
Wt1 = (Ip(1,1) - Ip(3,3))/IW;
Wt2 = (Ip(3,3) - Ip(2,2))/IW;
Ww = 500*w0;%w0*(-Wt2 + 1000);
icw = [0, 0.02, w0]';
ica = [0, 0, 0]'*dtr;
RB = rotation(1, ica(3))*rotation(2, ica(2))*rotation(3, ica(1));
PHI = acos(0.5*(trace(RB) - 1));
%
if abs(PHI) > eps
    EHAT = (1/2/sin(PHI))*[RB(2,3) - RB(3,2); RB(3,1) - RB(1,3); RB(1,2) - RB(2,1)];
else
    EHAT = [0; 0; 0];
end
CP = cos(PHI/2);
SP = sin(PHI/2);
icq = [CP; EHAT*SP];
this = linspace(0,60,400);
[t x] = ode45('eulerdynusym', this, [icq; icw]);

I = [7; 0; 0];
J = [0; 7; 0];
K = [0; 0; 6]; 

mctr = 1;
for i = 1:length(t)
%     rotmat = rotation(3, x(i,3))*rotation(1, x(i,2))*rotation(3, x(i,1));
    q4t = x(i,1); q1t = x(i,2); q2t = x(i,3); q3t = x(i,4);
    q13t = [q1t; q2t; q3t];
    qcr = [0 -q3t q2t; q3t 0 -q1t; -q2t q1t 0];
    rotmat = (q4t^2 - q13t'*q13t)*eye(3) + 2*q13t*q13t' - 2*q4t*qcr;
    for j = 1:size(xe,1)
        for k = 1:size(xe,2)
            vold = [xe(j,k); ye(j,k); ze(j,k)];
            vnew = rotmat'*vold;
            xen(j,k) = vnew(1); yen(j,k) = vnew(2); zen(j,k) = vnew(3);
        end
    end
    In = rotmat'*I; Jn = rotmat'*J; Kn = rotmat'*K;
    Khis(i,:) = Kn';
    Jhis(i,:) = Jn';
    Ihis(i,:) = In';    
    figure(1);
    surf(xen, yen, zen); hold on;
    plot3([0 In(1)], [0 In(2)], [0 In(3)], 'm', 'linewidth', 3);
    plot3([0 Jn(1)], [0 Jn(2)], [0 Jn(3)], 'm', 'linewidth', 3);
    plot3([0 Kn(1)], [0 Kn(2)], [0 Kn(3)], 'm', 'linewidth', 3);    
    plot3([0 K(1)], [0 K(2)], [0 K(3)], 'k', 'linewidth', 3);
    plot3(Khis(:,1), Khis(:,2), Khis(:,3), 'g.');
    text(In(1) + 0.2, In(2) + 0.2, In(3) + 0.2, 'b_1', 'fontweight', 'bold');
    text(Jn(1) + 0.2, Jn(2) + 0.2, Jn(3) + 0.2, 'b_2', 'fontweight', 'bold');
    text(Kn(1) + 0.2, Kn(2) + 0.2, Kn(3) + 0.2, 'b_3', 'fontweight', 'bold');
    title(['Time: ' num2str(t(i)), 's']);

%     plot3(Jhis(:,1), Jhis(:,2), Jhis(:,3), 'g.');
%     plot3(Ihis(:,1), Ihis(:,2), Ihis(:,3), 'g.');    
%     plot3(xen(37,37), yen(37,37), zen(37,37), 'kp', 'MarkerFaceColor', 'k')
    axis([-maxr maxr -maxr maxr -maxr maxr]);
    hold off
    getframe();
    pause(0.01);
%     movi(mctr) = getframe;
%     mctr = mctr + 1;
end

figure(2);
plot(this, x(:,5));
hold on
plot(this, x(:,6), 'r');
plot(this, x(:,7), 'k');
axis([0, this(end), (1.5*min(min(x(:,5:7))) -0.1), (1.1*max(max(x(:,5:7))) +0.1)]);
legend('\omega_1', '\omega_2', '\omega_3')
