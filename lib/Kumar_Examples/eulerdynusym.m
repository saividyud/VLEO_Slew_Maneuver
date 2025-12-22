function dots = eulerdynusym(t, x)
global Ip sig IW Ww;
q4 = x(1);
q1 = x(2); q2 = x(3); q3 = x(4);
w = x(5:end);

% ps = x(1); th = x(2); ph = x(3);
% cph = cos(ph); sph = sin(ph);
% w = x(4:end);
% sth = sin(th);
% if abs(sth) <= eps
%     psdot = 0;
% else
%     psdot = (w(1)*sph + w(2)*cph)/sth;
% end
% thdot = (w(1)*cph - w(2)*sph);
% phdot = w(3) - psdot*cos(th);
Bq = [q4 -q1 -q2 -q3;
    q1 q4 -q3 q2;
    q2 q3  q4 -q1;
    q3 -q2  q1 q4];

qdots = 0.5*Bq*[0; w];
% wcr = [0 -w(3) w(2);
%     w(3) 0 -w(1);
%     -w(2) w(1) 0];
% wp = w(3)*(sig - 1);
% wdots = inv(Ip)*(-wcr*Ip*w);
% wdots(1,1) = -w(2)*wp;
% wdots(2,1) = w(1)*wp;
% wdots(3,1) = 0;
if t < 10
    Wwt = 0;
    k = [0;0;0];
else
    Wwt = Ww;
%     k = [.5; .5; .5];
    k = [1; 1; 1];
end
wdots(1,1) = ( -w(2)*w(3)*(Ip(3,3) - Ip(2,2)) - IW*Wwt*w(2) )/Ip(1,1) - k(1)*w(1);
wdots(2,1) = ( -w(1)*w(3)*(Ip(1,1) - Ip(3,3)) + IW*Wwt*w(1) )/Ip(2,2) - k(2)*w(2);
wdots(3,1) = -w(1)*w(2)*(Ip(2,2) - Ip(1,1))/Ip(3,3) - k(3)*w(3);
dots = [qdots; wdots];
% dots = [psdot; thdot; phdot; wdots];