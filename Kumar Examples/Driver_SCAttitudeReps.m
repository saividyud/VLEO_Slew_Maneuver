clear; close all; clc;
global w;
O = [0; 0; 0];
i1 = [1; 0; 0];
i2 = [0; 1; 0];
i3 = [0; 0; 1];

plot3([0 1], [0 0], [0 0], 'linewidth', 2); hold on
plot3([0 0], [0 1], [0 0], 'linewidth', 2);
plot3([0 0], [0 0], [0 1], 'linewidth', 2);
text(1.05,0,0, 'i_1', 'color', 'b');
text(0,1.05,0, 'i_2', 'color', 'b');
text(0,0,1.05, 'i_3', 'color', 'b');
% b1 = randn(3,1);
% b1 = b1/norm(b1);
% b2 = randn(2,1);
% mis2 = dot(b1(1:2), b2)/-b1(end);
% b2 = [b2; mis2];
% b2 = b2/norm(b2);
% b3 = randn(1);
% rhs = -[b3*b1(end); b3*b2(end)];
% A = [b1(1:2)'; b2(1:2)'];
% b3_12 = inv(A)*rhs;
% b3 = [b3_12; b3];
% b3 = b3/norm(b3);
dtr = pi/180;
th = 19.7*dtr;
ps = 35*dtr;
ph = 0.*dtr;
Rbi = FRE(1, ph)*FRE(2, th)*FRE(3, ps);
b1 = Rbi(1,:)';
b2 = Rbi(2,:)';
b3 = Rbi(3,:)';

drawline(O, b1, 'k', 2);
drawline(O, b2, 'k', 2);
drawline(O, b3, 'k', 2);
text(b1(1) + 0.05,b1(2) + 0.05,b1(3), 'b_1', 'color', 'k');
text(b2(1) + 0.05,b2(2) + 0.05,b2(3), 'b_2', 'color', 'k');
text(b3(1)+ 0.05,b3(2)+ 0.05,b3(3), 'b_3', 'color', 'k');
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);

w = [-0.; 1; 0];
% Rbi = [dot(b1,i1), dot(b1, i2), dot(b1, i3);
%     dot(b2,i1), dot(b2, i2), dot(b2, i3);
%     dot(b3,i1), dot(b3, i2), dot(b3, i3)];
% th2 = asin(-Rbi(1,3));
% cth2 = cos(th2);    %has GOT to be positive
% th1 = atan2(Rbi(1,2), Rbi(1,1));
% th3 = atan2(Rbi(2,3), Rbi(3,3));
% Rm = R1(th3)*R2(th2)*R3(th1);
% if norm(Rm - Rbi) < 1e-9
% else
%     th2 = pi - th2;
%     cth2 = cos(th2);
%     %     if cth2 < 0
%     %         th1 = atan2(-Rbi(1,2), -Rbi(1,1));
%     %         th3 = atan2(-Rbi(2,3), -Rbi(3,3));
%     %     else
%     th1 = atan2(Rbi(1,2)/cth2, Rbi(1,1)/cth2);
%     th3 = atan2(Rbi(2,3), Rbi(3,3));
%     %     end
% end
% Rm = R1(th3)*R2(th2)*R3(th1);
% if norm(Rm - Rbi) < 1e-9
%     fprintf('good\n');
% else
%     fprintf('returned\n');
%     return;
% end

%%%
% 3-2-1 Euler angle representation of the rotation matrix
%%%
Bth0 = [-sin(th), 0, 1;
    cos(th)*sin(ph), cos(ph), 0;
    cos(th)*cos(ph), -sin(ph), 0];
icthdot = inv(Bth0)*w;
icth = [ps; th; ph];
t = linspace(0, 10, 100);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[tt xth] = ode45('eulerdynamics', t, icth, options);
figure(2);
subplot(311)
plot(tt, xth(:,1));
title('Time history of 3-2-1 Euler Angle Sequence')
subplot(312)
plot(tt, xth(:,2));
hold on
plot(tt, pi/2*ones(length(tt),1), 'r');
subplot(313)
plot(tt, xth(:,3));
for i = 1:length(tt)
    Rth = FRE(1, xth(i,3))*FRE(2, xth(i,2))*FRE(3, xth(i,1));
    alpthis(i,1) = acos(0.5*(trace(Rth) - 1));
    detthis(i,1) = det(Rth);
end

%%%
% Quaternion representation of the rotation matrix
%%%
alph = acos(0.5*(trace(Rbi) - 1));
salph = sin(alph);
calph = cos(alph);
ehat = [Rbi(2,3) - Rbi(3,2); Rbi(3,1) - Rbi(1,3); Rbi(1,2) - Rbi(2,1)]/(2*salph);
ecross = [0 -ehat(3) ehat(2); ehat(3) 0 -ehat(1); -ehat(2) ehat(1) 0];

Reuler = eye(3)*calph + (1 - calph)*ehat*ehat' - salph*ecross;

q1eu = ehat(1)*sin(alph/2);
q2eu = ehat(2)*sin(alph/2);
q3eu = ehat(3)*sin(alph/2);
q4eu = cos(alph/2);
%%% Stanley's method:
trR = trace(Rbi);
q4s = 0.25*(1 + trR);
q1s = 0.25*(1 + 2*Rbi(1,1) - trR);
q2s = 0.25*(1 + 2*Rbi(2,2) - trR);
q3s = 0.25*(1 + 2*Rbi(3,3) - trR);
qms = max([q1s q2s q3s q4s]);
if abs(qms - q1s) < 1e-6
    q1 = sqrt(q1s);
    q2 = 0.25*(Rbi(1,2) + Rbi(2,1))/q1;
    q3 = 0.25*(Rbi(1,3) + Rbi(3,1))/q1;
    q4 = 0.25*(Rbi(2,3) - Rbi(3,2))/q1;
elseif abs(qms - q2s) < eps
    q2 = sqrt(q2s);
    q1 = 0.25*(Rbi(1,2) + Rbi(2,1))/q2;
    q3 = 0.25*(Rbi(2,3) + Rbi(3,2))/q2;
    q4 = 0.25*(Rbi(3,1) - Rbi(1,3))/q2;
elseif abs(qms - q3s) < eps
    q3 = sqrt(q3s);
    q1 = 0.25*(Rbi(1,3) + Rbi(3,1))/q3;
    q2 = 0.25*(Rbi(2,3) + Rbi(3,2))/q3;
    q4 = 0.25*(Rbi(1,2) - Rbi(2,1))/q3;
else
    q4 = sqrt(q4s);
    q1 = 0.25*(Rbi(2,3) - Rbi(3,2))/q4;
    q2 = 0.25*(Rbi(3,1) - Rbi(1,3))/q4;
    q3 = 0.25*(Rbi(1,2) - Rbi(2,1))/q4;
end

icq = [q4; q1; q2; q3];
[tq, xq] = ode45('quaterniondynamics', t, icq, options);
figure(3)
subplot(411)
plot(tq, xq(:,2)); ylabel('q_1');
title('Time history of Quaternions (Euler Parameters)')
subplot(412)
plot(tq, xq(:,3)); ylabel('q_2');
% hold on
% plot(tq, pi/2*ones(length(tq),1), 'r');
subplot(413)
plot(tq, xq(:,4)); ylabel('q_3');
subplot(414)
plot(tq, xq(:,1)); ylabel('q_4'); %xlabel('Time, s');

for i = 1:length(tq)
    q4h = xq(i,1); q1h = xq(i,2); q2h = xq(i,3); q3h = xq(i,4);
    Rt = [(q1h^2 - q2h^2 - q3h^2 + q4h^2), 2*(q1h*q2h + q3h*q4h), 2*(q1h*q3h - q2h*q4h);
        2*(q1h*q2h - q3h*q4h), (-q1h^2 + q2h^2 - q3h^2 + q4h^2), 2*(q2h*q3h + q1h*q4h);
        2*(q1h*q3h + q2h*q4h), 2*(q2h*q3h - q1h*q4h), (-q1h^2 - q2h^2 + q3h^2 + q4h^2)];
    detqhis(i,1) = det(Rt);
    alpqhis(i,1) = acos(0.5*(trace(Rt) - 1));
    qconshis(i,1) = q1h^2 + q2h^2 + q3h^2 + q4h^2;
    th2_1(i,1) = asin(-Rt(1,3));
    th2_2(i,1) = pi - th2_1(i,1);
    th1_1(i,1) = atan2(Rt(1,2), Rt(1,1));
    th3_1(i,1) = atan2(Rt(2,3), Rt(3,3));
    th1_2(i,1) = atan2(-Rt(1,2), -Rt(1,1));
    th3_2(i,1) = atan2(-Rt(2,3), -Rt(3,3));    
end
figure(4)
plot(tq, alpqhis/dtr);
set(gca, 'fontsize', 14);
title('Time history of Principal Angle, \alpha');
figure(5)
plot(t, abs(detqhis - 1), 'k');
hold on
plot(t, abs(qconshis - 1), 'b');
set(gca, 'fontsize', 14);
legend('Determinant(R)', 'Quaternion Constraint', 'Location', 'NorthWest');
title('Sanity Checks for Quaternions');

%%%%
% CRP 
%%%
s1 = q1/q4;
s2 = q2/q4;
s3 = q3/q4;

ics = [s1; s2; s3];
[ts xs] = ode45('crpdynamics', t, ics);
figure(6);
subplot(311)
plot(ts, xs(:,1)); ylabel('s_1');
title('Time history of Classical Rodrigues Parameters')
subplot(312)
plot(ts, xs(:,2)); ylabel('s_2');
hold on
% plot(tt, pi/2*ones(length(tt),1), 'r');
subplot(313)
plot(ts, xs(:,3)); ylabel('s_3');

for i = 1:length(ts)
    st = [xs(i,1); xs(i,2); xs(i,3)];
    stcross = [0 -st(3) st(2); st(3) 0 -st(1); -st(2) st(1) 0]; 
    Rt = ((1 - st'*st)*eye(3) + 2*st*st' - 2*stcross)/(1 + st'*st);
    alpshis(i,1) = acos(0.5*(trace(Rt) - 1));
    detshis(i,1) = det(Rt);
end


%%%
% MRP representation
%%%
sg1 = q1/(1 + q4);
sg2 = q2/(1 + q4);
sg3 = q3/(1 + q4);
icsg = [sg1; sg2; sg3];
[tsg, xsg] = ode45('mrpdynamics', t, icsg);
figure(7);
subplot(311)
plot(tsg, xsg(:,1)); ylabel('\sigma_1');
set(gca, 'fontsize', 14);
title('Time history of Modified Rodrigues Parameters')
subplot(312)
plot(tsg, xsg(:,2));  ylabel('\sigma_2');
hold on
% plot(tt, pi/2*ones(length(tt),1), 'r');
subplot(313)
plot(tsg, xsg(:,3));  ylabel('\sigma_3');

for i = 1:length(tsg)
    sgt = [xsg(i,1); xsg(i,2); xsg(i,3)];
    sgtcross = [0 -sgt(3) sgt(2); sgt(3) 0 -sgt(1); -sgt(2) sgt(1) 0]; 
    Rtsg = eye(3) + (8*sgtcross^2 - 4*(1 - sgt'*sgt)*sgtcross)/(1 + sgt'*sgt)^2;
    alpsghis(i,1) = acos(0.5*(trace(Rtsg) - 1));
    detsghis(i,1) = det(Rtsg);
end

figure
plot(tq, abs(detqhis - 1), 'b');
hold on;
plot(ts, abs(detshis - 1), 'k-.', 'linewidth', 2);
hold on
plot(tsg, abs(detsghis - 1), 'r');
set(gca, 'fontsize', 14);
legend('Quaternions', 'CRPs', 'MRPs', 'Location', 'NorthWest');
title('Time history of matrix determinant');

figure
plot(tq, th1_1, 'k');
hold on;
% plot(tq, alpqhis);
% hold on;
% plot(tsg, alpsghis, 'r');
plot(tq, th2_1, 'b', 'linewidth', 2);
plot(tq, th3_1, 'r');
plot(tq, pi/2*ones(length(tq),1), 'm');
legend('\psi', '\theta', '\phi');
