%Principal body frame
clear; clc; close all;

global g dtr ftm ftm2 stk lbtk lbftn psftpa psitpa mps2mph;
dtr = pi/180;
ftm = 0.3048;
ftm2 = ftm^2;
mps2mph = 2.23694;
stk = 14.5939;
lbtk = 0.453592;
lbftn = 4.44822;
psftpa = 47.88026;
psitpa = 6.89476e3;

m = 630000*lbftn;
Ixx = 18.2e6*stk*ftm2;
Iyy = 33.1e6*stk*ftm2;
Izz = 49.7e6*stk*ftm2;
Ixz = -1.56e6*stk*ftm2;
Ixy = -5e6; %-0.04e6*stk*ftm2;
Iyz = 1e6*stk*ftm2;
Imat = [Ixx -Ixy -Ixz;
    -Ixy Iyy -Iyz;
    -Ixz -Iyz Izz];



%%%PLOT
O = [0; 0; 0];  %common origin
F = eye(3); %generic x, y, z axes as columns
ps = 20*dtr; th = 65*dtr; ph = 30*dtr;
RBI = FRE(1, ph)*FRE(2, th)*FRE(3, ps);
b1 = RBI(1,:);
b2 = RBI(2,:);
b3 = RBI(3,:);

%%%obtain principal body axis
[eve, eva] = eig(Imat);
for i = 1:3
    [mx, mp] = max(abs(eve(i,:)));
    if eve(i,mp) < 0;
        eve(i,:) = -eve(i,:);
    end
end
% if det(eve) < 0
%     eve(2,:) = -eve(2,:);
%     eve(3,:) = -eve(3,:);
% end
RBBstar = eve;
RBstarI = RBBstar'*RBI;
b1star = RBstarI(1,:);
b2star = RBstarI(2,:);
b3star = RBstarI(3,:);

%%% plot original
drawline(O, F(1,:), 'k', 2); hold on;
drawline(O, F(2,:), 'k', 2);
drawline(O, F(3,:), 'k', 2);
%%%plot final
drawline(O, b1, 'r', 2);
drawline(O, b2, 'r', 2);
drawline(O, b3, 'r', 2);
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
%view(145, 35);
set(gca, 'fontsize', 16, 'fontweight', 'bold');
text(b1(1) + 0.05,b1(2) + 0.05,b1(3), 'b_1', 'color', 'r');
text(b2(1) + 0.05,b2(2) + 0.05,b2(3), 'b_2', 'color', 'r');
text(b3(1) + 0.05,b3(2) + 0.05,b3(3), 'b_3', 'color', 'r');
text(0.2, 0.2, 0.0, 'O', 'fontsize', 12, 'fontweight', 'bold', 'color', 'k');
text(0., 0., 1.05, 'i_3', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
text(1.05, 0., 0., 'i_1', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
text(0., 1.05, 0., 'i_2', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
title('Demonstration: Euler Principal Axis Theorem');
%%%
% get airplane
%%%
aclen = .75;
acwid = .12;
achgt = .12;
wspan = .75;
wchord = .12;
tchord = .06;
tspan = .40;
I = [0.65*aclen; 0; 0];
J = [0; -0.65*wspan; 0];
K = [0; 0; -0.3*wspan]; 
BCM = [0, 0, 0];
[xa, ya, za] = ellipsoid(BCM(1), BCM(2), BCM(3), aclen/2, acwid/2, achgt/2, 50);    %fuselage
[xw, yw, zw] = ellipsoid(BCM(1) + aclen/20, BCM(2), BCM(3), wchord/2, wspan/2, wchord/10, 50);   %wing
[xt, yt, zt] = ellipsoid(BCM(1) - aclen/2.5, BCM(2), BCM(3), tchord/2, tspan/2, tchord/5, 50);  %tail
[xv, yv, zv] = ellipsoid(BCM(1) - aclen/2.5, BCM(2), BCM(3), tchord/2, tchord/5, tspan/2);  %vertical fin
fz = find(zv < BCM(3));
xv(fz) = BCM(1); yv(fz) =  BCM(2); zv(fz) = BCM(3);
[xan, yan, zan] = rotateobj(xa, ya, za, RBI');
[xwn, ywn, zwn] = rotateobj(xw, yw, zw, RBI');
[xtn, ytn, ztn] = rotateobj(xt, yt, zt, RBI');
[xvn, yvn, zvn] = rotateobj(xv, yv, zv, RBI');
surf(xan, yan, zan);
hold on
surf(xwn, ywn, zwn);
surf(xtn, ytn, ztn);
%surf(xvn, yvn, -zvn);
colormap('jet')
shading interp
alpha('direct');
alphamap([.1; 1])
camlight;
lighting gouraud;

drawline(O, b1star, 'b', 2);
drawline(O, b2star, 'b', 2);
drawline(O, b3star, 'b', 2);
text(b1star(1) + 0.05,b1star(2) + 0.05,b1star(3), 'b_1^*', 'color', 'b');
text(b2star(1) + 0.05,b2star(2) + 0.05,b2star(3), 'b_2^*', 'color', 'b');
text(b3star(1) + 0.05,b3star(2) + 0.05,b3star(3), 'b_3^*', 'color', 'b');

axis equal;
