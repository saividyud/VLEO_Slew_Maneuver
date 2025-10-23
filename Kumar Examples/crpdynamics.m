function dots = crpdynamics(t, x)
global w;

s1 = x(1); s2 = x(2); s3 = x(3);
s = [s1; s2; s3];
scross = [0 -s3 s2; s3 0 -s1; -s2 s1 0];
dots = 0.5*(eye(3) + scross + s*s')*w;