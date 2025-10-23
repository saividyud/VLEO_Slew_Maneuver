function q = stanley(R)

tr = trace(R);


%%%first compute squares of each quaternion
bb(1) = 0.25*(1 + tr);
bb(2) = 0.25*(1 + 2*R(1,1) - tr);
bb(3) = 0.25*(1 + 2*R(2,2) - tr);
bb(4) = 0.25*(1 + 2*R(3,3) - tr);

%%%find the max!
[~, mxpos] = max(bb);

%%%identify the quaternion with max sq-value
q(mxpos) = sqrt(bb(mxpos));

stan(1) = (R(2,3) - R(3,2))/4;
stan(2) = (R(3,1) - R(1,3))/4;
stan(3) = (R(1,2) - R(2,1))/4;
stan(4) = (R(1,2) + R(2,1))/4;
stan(5) = (R(1,3) + R(3,1))/4;
stan(6) = (R(2,3) + R(3,2))/4;

if mxpos == 1
    sindex = [1;2;3];
    qindex = [2;3;4];
elseif mxpos == 2
    sindex = [1;4;5];
    qindex = [1;3;4];
elseif mxpos == 3
    sindex = [2;4;6];
    qindex = [1;2;4];
else
    sindex = [3;5;6];
    qindex = [1;2;3];
end

q(qindex) = stan(sindex)/q(mxpos);
if q(1) < 0
    q = -q;
end
q = q';
