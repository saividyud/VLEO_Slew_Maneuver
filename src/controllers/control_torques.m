function LC = control_torques(t, X)
% CONTROL_TORQUES calculates the torque on the satellite in the inertial
% reference frame based off of different actuators.
% 
% Parameters
% ----------
% t : float
%   Time of simulation
% X : 13x1 vector
%   State of satellite
% 
% Returns
% -------
% LC : 3x1 vector
%   Controlling torques

    ICB = 2/5*83*(.58/2)^2*[1 0 0 ; 0 1 0 ;0 0 1]; % [kg m^2]

    %calculating u
    q_r = [1; 0; 0; 0];
    w_r = [0;0;0];
    wdot_r = [0;0;0];
    P = [.15 0 0; 0 .15 0; 0 0 .15];
    Kp = [.0025 0 0 ; 0 .0025 0 ; 0 0 .0025];
    delw = X(11:13) - w_r;
    u = -Kp * (q_r(2:4) - X(8:10));
    u = u- P * delw;
    u = u + ICB * wdot_r;
    u = u - cross(X(11:13),w_r);
    % u = u + X(11:13)' * ICB * X(11:13)
    % if abs(norm(u)) >= 1
    %     disp("nonsense")
    %     disp(norm(u))
    % end

    LC = u;
    % LC = [0, 0, 0]';

    % t_start = 2.5 * 60; % Start 5 minutes into simulation
    % duration = 5 * 60; % Lasts for 5 minutes
    % 
    % if (t_start < t) && (t <= t_start + duration/2)
    %     LC = [0.0001; 0; 0];
    % elseif (t_start + duration/2 < t) && (t <= t_start + duration)
    %     LC = [-0.0001; 0; 0];
    % else
    %     LC = [0; 0; 0]; % No control torques outside the specified time range
    % end

end
