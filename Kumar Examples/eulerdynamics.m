function dots = eulerdynamics(t,x)
global w;
ps = x(1);
th = x(2);
ph = x(3);

Bth = [-sin(th), 0, 1;
    cos(th)*sin(ph), cos(ph), 0;
    cos(th)*cos(ph), -sin(ph), 0];
dots = inv(Bth)*w;
