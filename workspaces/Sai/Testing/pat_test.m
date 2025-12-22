clear
clc
close all

%% Defining orbital reference frame axes
step = 1;

pos = rs(step, :);
vel = vs(step, :);

% o_1 is the radially outward direction
o_1 = pos' / norm(pos);

% o_2 is the direction of velocity
o_2 = vel' / norm(vel);

% o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
o_3 = cross(o_1, o_2);

% Normalize o_3 to ensure it is a unit vector
o_3 = o_3 / norm(o_3);

