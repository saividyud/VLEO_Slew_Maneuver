function beta = QfromDCM(R_BI)
% QfromDCM calculates the quaternion given a DCM from inertial to body 
% frame.
%
% Parameters
% ----------
% R_BI : 3x3 matrix
%   Rotation DCM for converting from inertial frame to body frame
%
% Returns
% -------
% beta : 4x1 quaternion
%   Output quaternion

tr = trace(R_BI);


%%%first compute squares of each quaternion
bb(1) = 0.25*(1 + tr);
bb(2) = 0.25*(1 + 2*R_BI(1,1) - tr);
bb(3) = 0.25*(1 + 2*R_BI(2,2) - tr);
bb(4) = 0.25*(1 + 2*R_BI(3,3) - tr);

%%%find the max!
[~, mxpos] = max(bb);

%%%identify the quaternion with max sq-value
beta(mxpos) = sqrt(bb(mxpos));

stan(1) = (R_BI(2,3) - R_BI(3,2))/4;
stan(2) = (R_BI(3,1) - R_BI(1,3))/4;
stan(3) = (R_BI(1,2) - R_BI(2,1))/4;
stan(4) = (R_BI(1,2) + R_BI(2,1))/4;
stan(5) = (R_BI(1,3) + R_BI(3,1))/4;
stan(6) = (R_BI(2,3) + R_BI(3,2))/4;

if mxpos == 1
    sindex = [1;2;3];
    betaindex = [2;3;4];
elseif mxpos == 2
    sindex = [1;4;5];
    betaindex = [1;3;4];
elseif mxpos == 3
    sindex = [2;4;6];
    betaindex = [1;2;4];
else
    sindex = [3;5;6];
    betaindex = [1;2;3];
end

beta(betaindex) = stan(sindex)/beta(mxpos);
if beta(1) < 0
    beta = -beta;
end
beta = beta';
