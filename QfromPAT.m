function [beta] = QfromPAT(e_hat, phi)
% PATfromQ calculates the principal axis and rotation about the principal
% axis according to Euler's Principal Axis Theorem given a quaternion.
% Parameters
% ----------
% e_hat : 3x1 unit vector
%   Rotation axis vector
% phi : float
%   Rotation angle about axis vector

% Returns
% -------
% beta : 4x1 quaternion
%   Input quaternion

    beta = [
        cosd(phi/2);
        e_hat(1) * sind(phi/2);
        e_hat(2) * sind(phi/2);
        e_hat(3) * sind(phi/2)
    ];
    
    % Negate the quaternion if the first element is negative
    if beta(1) < 0
        beta = -beta;
    end

    % Ensure the output quaternion is normalized
    beta = beta / norm(beta);

end