function [e_hat, phi] = PATfromQ(beta)
% PATfromQ calculates the principal axis and rotation about the principal
% axis according to Euler's Principal Axis Theorem given a quaternion.
% Parameters
% ----------
% beta : 4x1 quaternion
%   Input quaternion

% Returns
% -------
% e_hat : 3x1 unit vector
%   Rotation axis vector
% phi : float
%   Rotation angle about axis vector
    
    % Negate the quaternion if the first element is negative
    if beta(1) < 0
        beta = -beta;
    end

    % Checking magnitude of quaternion
    if abs(norm(beta) - 1) > 1e-6
        beta = beta / norm(beta); % Normalize the quaternion
    end

    % First element is simply cos(phi/2)
    phi = 2 * acosd(beta(1));

    % Ensure phi is not zero to avoid division by zero
    if phi == 0
        e_hat = [1; 0; 0]; % Default to x-axis if no rotation
    else
        % Next three elements
        e_hat = [
            beta(2) / sind(phi/2);
            beta(3) / sind(phi/2);
            beta(4) / sind(phi/2);
        ];
    end

end