function R_BI = DCMfromQ(beta)
% DCMfromQ calculates the DCM from inertial to body frame given a 
% quaternion.
% Parameters
% ----------
% beta : 4x1 quaternion
%   Input quaternion

% Returns
% -------
% R_BI : 3x3 matrix
%   Rotation DCM for converting from inertial frame to body frame
    
    % Extracting quaternions
    b0 = beta(1);
    b1 = beta(2);
    b2 = beta(3);
    b3 = beta(4);

    % Calculating DCM
    R_BI = [
        b0^2 + b1^2 - b2^2 - b3^2, 2*(b1*b2 + b0*b3), 2*(b1*b3 - b0*b2);
        2*(b1*b2 - b0*b3), b0^2 - b1^2 + b2^2 - b3^2, 2*(b2*b3 + b0*b1);
        2*(b1*b3 + b0*b2), 2*(b2*b3 - b0*b1), b0^2 - b1^2 - b2^2 + b3^2
    ];
    
end