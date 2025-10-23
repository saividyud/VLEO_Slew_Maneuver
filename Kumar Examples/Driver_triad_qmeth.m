clear; clc; close all;

sb = [0.8273 0.5541  -0.0920]';
mb = [-0.8285 0.5522  -0.0955]';

si = [-0.1517  -0.9669 0.2050]';
mi = [-0.8393  0.4494 -0.3044]';

t1b = sb;
t2b = cross(sb, mb);
t2b = t2b/norm(t2b);    %renormalize to ensure unit vector, i.e. norm = 1;
t3b = cross(t1b, t2b);
t3b = t3b/norm(t3b);    %renormalize to ensure unit vector, i.e. norm = 1;


t1i = si;
t2i = cross(si, mi);
t2i = t2i/norm(t2i);    %renormalize to ensure unit vector (norm = 1)
t3i = cross(t1i, t2i);
t3i = t3i/norm(t3i);    %renormalize to unit norm = 1;

RBT = [t1b t2b t3b];
RIT = [t1i t2i t3i];

%%%
% Final step
RBI = RBT*RIT';

%%%
% plot
%%%
O = [0; 0; 0];
i1 = [1; 0; 0];
i2 = [0; 1; 0];
i3 = [0; 0; 1];
figure(1)
plot3([0 1], [0 0], [0 0], 'linewidth', 2); hold on
plot3([0 0], [0 1], [0 0], 'linewidth', 2);
plot3([0 0], [0 0], [0 1], 'linewidth', 2);
text(1.05,0,0, 'i_1', 'color', 'b');
text(0,1.05,0, 'i_2', 'color', 'b');
text(0,0,1.05, 'i_3', 'color', 'b');
drawline(O, si, 'y', 3);
drawline(O, mi, 'r', 3);
text(si(1), si(2),  si(3) - 0.1, 's_i', 'color', 'm');
text(mi(1), mi(2),  mi(3) - 0.1, 'm_i', 'color', 'r');
drawline(O, RBI(:,1), 'g', 3);
drawline(O, RBI(:,2), 'g', 3);
drawline(O, RBI(:,3), 'g', 3);
text(RBI(1,1), RBI(2,1), RBI(3,1) - 0.1, 'b_1', 'color', 'g'); 
text(RBI(1,2), RBI(2,2), RBI(3,2) - 0.1, 'b_2', 'color', 'g'); 
text(RBI(1,3), RBI(2,3), RBI(3,3) - 0.1, 'b_3', 'color', 'g'); 
title('Triad Algorithm');


figure(2)
plot3([0 1], [0 0], [0 0], 'g', 'linewidth', 2); hold on
plot3([0 0], [0 1], [0 0], 'g', 'linewidth', 2);
plot3([0 0], [0 0], [0 1], 'g', 'linewidth', 2);
text(1.05,0,0, 'b_1', 'color', 'g');
text(0,1.05,0, 'b_2', 'color', 'g');
text(0,0,1.05, 'b_3', 'color', 'g');
drawline(O, sb, 'y', 3);
drawline(O, mb, 'r', 3);
text(sb(1), sb(2),  sb(3) - 0.1, 's_b', 'color', 'm');
text(mb(1), mb(2),  mb(3) - 0.1, 'm_b', 'color', 'r');
drawspacecraft(eye(3));
figure(1)
drawline(O, t1i, 'k', 2);
drawline(O, t2i, 'k', 2);
drawline(O, t3i, 'k', 2);
text(t1i(1), t1i(2),  t1i(3) - 0.1, 't1_i', 'color', 'k');
text(t2i(1), t2i(2),  t2i(3) - 0.1, 't2_i', 'color', 'k');
text(t3i(1), t3i(2),  t3i(3) - 0.1, 't3_i', 'color', 'k');
drawspacecraft(RBI);

% figure(2)
% drawline(O, t1b, 'k', 2);
% drawline(O, t2b, 'k', 2);
% drawline(O, t3b, 'k', 2);
% text(t1b(1), t1b(2),  t1b(3) - 0.1, 't1_i', 'color', 'k');
% text(t2b(1), t2b(2),  t2b(3) - 0.1, 't2_i', 'color', 'k');
% text(t3b(1), t3b(2),  t3b(3) - 0.1, 't3_i', 'color', 'k');

%%%
% sensor error checks
%%%
err_s = norm(sb - RBI*si);
err_m = norm(mb - RBI*mi);

%%%
% q-method
%%%
N = 2;
w(1) = 10;   %corresp. sun vector
w(2) = 1;   %corresp. magnetic field vector
u(1).bo = sb; 
u(1).in = si;
u(2).bo = mb; 
u(2).in = mi;

B = zeros(3);
for i = 1:N
    B = B + w(i)*(u(i).bo)*(u(i).in)';
end
S = B + B';
Z = [B(2,3) - B(3,2), B(3,1) - B(1,3), B(1,2) - B(2,1)]';
sig = trace(B);
K = [(S - sig*eye(3)), Z];
K = [K; [Z' sig]];
[eigv eigd] = eig(K);
deigd = diag(eigd);
[lam_opt pos] = max(deigd);
qopt = eigv(:, pos);

skewq = [0 -qopt(3) qopt(2);
    qopt(3) 0 -qopt(1);
    -qopt(2) qopt(1) 0];
RBI_qmethod = (qopt(4)^2 - norm(qopt(1:3))^2)*eye(3) + 2*qopt(1:3)*qopt(1:3)' - 2*qopt(4)*skewq;


figure(3)
plot3([0 1], [0 0], [0 0], 'linewidth', 2); hold on
plot3([0 0], [0 1], [0 0], 'linewidth', 2);
plot3([0 0], [0 0], [0 1], 'linewidth', 2);
text(1.05,0,0, 'i_1', 'color', 'b');
text(0,1.05,0, 'i_2', 'color', 'b');
text(0,0,1.05, 'i_3', 'color', 'b');
drawline(O, si, 'y', 3);
drawline(O, mi, 'r', 3);
text(si(1), si(2),  si(3) - 0.1, 's_i', 'color', 'm');
text(mi(1), mi(2),  mi(3) - 0.1, 'm_i', 'color', 'r');
drawline(O, RBI_qmethod(:,1), 'g', 3);
drawline(O, RBI_qmethod(:,2), 'g', 3);
drawline(O, RBI_qmethod(:,3), 'g', 3);
text(RBI_qmethod(1,1), RBI_qmethod(2,1), RBI_qmethod(3,1) - 0.1, 'b_1', 'color', 'g'); 
text(RBI_qmethod(1,2), RBI_qmethod(2,2), RBI_qmethod(3,2) - 0.1, 'b_2', 'color', 'g'); 
text(RBI_qmethod(1,3), RBI_qmethod(2,3), RBI_qmethod(3,3) - 0.1, 'b_3', 'color', 'g'); 
title('q-Method');
drawspacecraft(RBI_qmethod);




%%%
% QUEST Approximation
%%%
lam_qst = sum(w);
rod = inv((lam_qst + sig)*eye(3) - S)*Z;
skewr = [0 -rod(3) rod(2);
    rod(3) 0 -rod(1);
    -rod(2) rod(1) 0];
RBI_qst = ((1 - rod'*rod)*eye(3) + 2*rod*rod' - 2*skewr)/(1 + rod'*rod);
qqst = [rod; 1]/(sqrt(1 + norm(rod)^2));
skewqq = [0 -qqst(3) qqst(2);
    qqst(3) 0 -qqst(1);
    -qqst(2) qqst(1) 0];
RBI_qqst = (qqst(4)^2 - norm(qqst(1:3))^2)*eye(3) + 2*qqst(1:3)*qqst(1:3)' - 2*qqst(4)*skewq; 
figure(4)
plot3([0 1], [0 0], [0 0], 'linewidth', 2); hold on
plot3([0 0], [0 1], [0 0], 'linewidth', 2);
plot3([0 0], [0 0], [0 1], 'linewidth', 2);
text(1.05,0,0, 'i_1', 'color', 'b');
text(0,1.05,0, 'i_2', 'color', 'b');
text(0,0,1.05, 'i_3', 'color', 'b');
drawline(O, si, 'y', 3);
drawline(O, mi, 'r', 3);
text(si(1), si(2),  si(3) - 0.1, 's_i', 'color', 'm');
text(mi(1), mi(2),  mi(3) - 0.1, 'm_i', 'color', 'r');
drawline(O, RBI_qqst(:,1), 'g', 3);
drawline(O, RBI_qqst(:,2), 'g', 3);
drawline(O, RBI_qqst(:,3), 'g', 3);
text(RBI_qqst(1,1), RBI_qqst(2,1), RBI_qqst(3,1) - 0.1, 'b_1', 'color', 'g'); 
text(RBI_qqst(1,2), RBI_qqst(2,2), RBI_qqst(3,2) - 0.1, 'b_2', 'color', 'g'); 
text(RBI_qqst(1,3), RBI_qqst(2,3), RBI_qqst(3,3) - 0.1, 'b_3', 'color', 'g'); 
title('Quest Approximation')


