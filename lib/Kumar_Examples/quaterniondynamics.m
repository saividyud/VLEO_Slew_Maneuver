function dots = quaterniondynamics(t,x)
global w;

q4 = x(1); q1 = x(2); q2 = x(3); q3 = x(4); 
Bq = [q4 -q1 -q2 -q3;
    q1 q4 -q3 q2;
    q2 q3 q4 -q1;
    q3 -q2 q1 q4];

dots = 0.5*Bq*[0; w];