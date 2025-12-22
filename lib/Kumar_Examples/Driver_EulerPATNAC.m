%%% demonstration of euler's principal axis theorem
clear; clc; close all;
dtr = pi/180;
O = [0; 0; 0];  %common origin
F = eye(3); %generic x, y, z axes as columns
RS = [3; 2; 1];   %get rotation sequence

% %%%
% % second randomly generated set of unit vectors
% %%%
% b1 = randn(3,1);
% b1 = b1/norm(b1);
% b2 = randn(2,1);
% mis2 = dot(b1(1:2), b2)/- b1(end);
% b2 = [b2; mis2];
% b2 = b2/norm(b2);
% %b3 = randn(1);
% %rhs = -[b3*b1(end); b3*b2(end)];
% %A = [b1(1:2)'; b2(1:2)'];
% %b3_12 = inv(A)*rhs;
% %b3 = [b3_12; b3];
% b3 = cross(b1, b2);
% b3 = b3/norm(b3);
% b = [b1, b2, b3];
% RBI = b';
ps = 30*dtr; th = 65*dtr; ph = 130*dtr;
RBI = FRE(1, ph)*FRE(2, th)*FRE(3, ps);
b1 = RBI(1,:);
b2 = RBI(2,:);
b3 = RBI(3,:);

%%% determine initial conditions:
%%% short rotation (PHI < pi)
PHI = acos(0.5*(trace(RBI) - 1));
Be = 0.5*[(RBI(2,3) - RBI(3,2));
    (RBI(3,1) - RBI(1,3));
    (RBI(1,2) - RBI(2,1))]/sin(PHI);

%%%draw transition
%%% plot original
drawline(O, F(1,:), 'k', 2); hold on;
drawline(O, F(2,:), 'k', 2);
drawline(O, F(3,:), 'k', 2);
%%%plot final
drawline(O, b1, 'r', 2);
drawline(O, b2, 'r', 2);
drawline(O, b3, 'r', 2);
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
view(145, 35);
set(gca, 'fontsize', 16, 'fontweight', 'bold');
text(b1(1) + 0.05,b1(2) + 0.05,b1(3), 'b_1', 'color', 'r');
text(b2(1) + 0.05,b2(2) + 0.05,b2(3), 'b_2', 'color', 'r');
text(b3(1) + 0.05,b3(2) + 0.05,b3(3), 'b_3', 'color', 'r');
text(0.2, 0.2, 0.0, 'O', 'fontsize', 12, 'fontweight', 'bold', 'color', 'k');
text(0., 0., 1.05, 'i_3', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
text(1.05, 0., 0., 'i_1', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
text(0., 1.05, 0., 'i_2', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
title('Demonstration: Euler Principal Axis Theorem');
pause;

%%%
% get airplane
%%%
aclen = 75;
acwid = 12;
achgt = 12;
wspan = 75;
wchord = 12;
tchord = 6;
tspan = 40;
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

%%%
% rotate the aircraft to align with the body frame
%%%
[xan, yan, zan] = rotateobj(xa, ya, za, RBI);
[xwn, ywn, zwn] = rotateobj(xw, yw, zw, RBI);
[xtn, ytn, ztn] = rotateobj(xt, yt, zt, RBI);
[xvn, yvn, zvn] = rotateobj(xv, yv, zv, RBI);
In = RBI*I; Jn = RBI*J; Kn = RBI*K;


figure(1);
tlen = 50;
PHIH = linspace(0, PHI, tlen)';
Fhis = struct();
Fhis.x = zeros(tlen, 3);
Fhis.y = zeros(tlen, 3);
Fhis.z = zeros(tlen, 3);
zendhis = zeros(tlen, 3);
yendhis = zeros(tlen, 3);
xendhis = zeros(tlen, 3);
for i = 1:tlen
    Rcurr = GetREP(Be, PHIH(i));
    FC = Rcurr'*F;
    Fhis.x(i,:) = FC(:,1)';
    Fhis.y(i,:) = FC(:,2)';
    Fhis.z(i,:) = FC(:,3)';
    %%% plot original
    drawline(O, F(1,:), 'k', 2); hold on;
    drawline(O, F(2,:), 'k', 2);
    drawline(O, F(3,:), 'k', 2);
    %%%plot final
    drawline(O, b1, 'r', 2);
    drawline(O, b2, 'r', 2);
    drawline(O, b3, 'r', 2);
    %%% plot intermediate
    drawline(O, FC(:,1), 'g', 2);
    drawline(O, FC(:,2), 'g', 2);
    drawline(O, FC(:,3), 'g', 2);
    zendhis(i,:) = FC(:,3);
    yendhis(i,:) = FC(:,2);
    xendhis(i,:) = FC(:,1);
    plot3(zendhis(1:i,1), zendhis(1:i,2), zendhis(1:i,3), 'c', 'linewidth', 2);
    %plot3(yendhis(1:i,1), yendhis(1:i,2), yendhis(1:i,3), 'c', 'linewidth', 2);
    %plot3(xendhis(1:i,1), xendhis(1:i,2), xendhis(1:i,3), 'c', 'linewidth', 2);
    drawline(O, Be, 'm-.', 1);
    drawline(O, -Be, 'm-.', 1);
    axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
    set(gca, 'fontsize', 16, 'fontweight', 'bold');
    view(145, 35);
    %%%write text
    text(b1(1) + 0.05,b1(2) + 0.05,b1(3), 'b_1', 'color', 'r');
    text(b2(1) + 0.05,b2(2) + 0.05,b2(3), 'b_2', 'color', 'r');
    text(b3(1) + 0.05,b3(2) + 0.05,b3(3), 'b_3', 'color', 'r');
    text(0.2, 0.2, 0.0, 'O', 'fontsize', 12, 'fontweight', 'bold', 'color', 'k');
    text(0., 0., 1.05, 'i_3', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
    text(1.05, 0., 0., 'i_1', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');
    text(0., 1.05, 0., 'i_2', 'fontsize', 8, 'fontweight', 'bold', 'color', 'k');    
    title('Demonstration: Euler Principal Axis Theorem');
    hold off;
    pause(0.1);
end
