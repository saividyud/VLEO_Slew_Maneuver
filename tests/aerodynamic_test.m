clc

pos = X_i(1:3);
vel = X_i(4:6);

% o_1 is the radially outward direction
o_1 = pos' / norm(pos);

% o_2 is the direction of velocity
o_2 = vel' / norm(vel);

% o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
o_3 = cross(o_1, o_2) / norm(cross(o_1, o_2));

% o_2 actual
o_2 = cross(o_3, o_1);

R_OI = [o_1; o_2; o_3];

R_BI = quat2dcm(X_i(7:10)');

R_OB = R_OI * R_BI';

[aoa, aos] = dcm2alphabeta(R_OB)
