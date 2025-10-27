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

    t_start = 2.5 * 60; % Start 5 minutes into simulation
    duration = 5 * 60; % Lasts for 5 minutes

    if (t_start < t) && (t <= t_start + duration/2)
        LC = [0.0001; 0; 0];
    elseif (t_start + duration/2 < t) && (t <= t_start + duration)
        LC = [-0.0001; 0; 0];
    else
        LC = [0; 0; 0]; % No control torques outside the specified time range
    end

end