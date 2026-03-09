function tau = myController(t, X)
% testController calculates the torques on the satellite
% given the time and state.
%
% Parameters
% ----------
% t : float
%   Time [s]
% X : 13x1 vector
%   State of the system with the following attributes:
%       - Position (1:3)
%       - Velocity (4:6)
%       - Quaternion (7:10)
%       - Angular velocity (11:13)
%
% Returns
% -------
% tau : 3x1 vector
%   Torques applied by controller

    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

tau = zeros(3, 1);

end
