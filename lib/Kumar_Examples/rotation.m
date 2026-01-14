function R = rotation(direc, ang)

cang = cos(ang);
sang = sin(ang);

if direc == 1
    R = [1 0 0;
        0 cang sang;
        0 -sang cang];
elseif direc == 2
    R = [cang 0 -sang;
        0 1 0;
        sang 0 cang];
elseif direc == 3
    R = [cang sang 0;
        -sang cang 0;
        0 0 1];
else R = eye(3);
end