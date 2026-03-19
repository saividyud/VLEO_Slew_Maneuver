function F_history = allocate_thruster_forces(wrench_history, params, layoutOrMode)
    % allocate_thruster_forces Maps a 6-DOF wrench to 8 unilateral thruster forces.
    %
    % Inputs:
    %   wrench_history : N x 6 matrix [Fx, Fy, Fz, Tx, Ty, Tz] (Forces in N, Torques in Nm)
    %   params         : Simulation parameters struct containing bodyDims
    %   layoutOrMode   : (Optional) thruster layout struct, 'standard',
    %                    'optimized', or legacy logical.
    %
    % Outputs:
    %   F_history      : N x 8 matrix of thruster forces [N], where F >= 0.

    if nargin < 3 || isempty(layoutOrMode)
        layoutOrMode = "standard";
    end

    if size(wrench_history, 2) ~= 6
        error('VLEO_Slew_Maneuver:InvalidWrenchHistory', ...
            'wrench_history must be an N x 6 matrix [Fx Fy Fz Tx Ty Tz].')
    end

    if isstruct(layoutOrMode)
        W = layoutOrMode.wrench_matrix;
    else
        W = vleo.control.corner_thruster_layout(params, layoutOrMode).wrench_matrix;
    end

    nSamples = size(wrench_history, 1);
    F_history = zeros(nSamples, size(W, 2));
    options = optimset('Display', 'off');

    for idx = 1:nSamples
        w = wrench_history(idx, :)';
        % Using lsqnonneg to find F >= 0 such that W*F ~ w
        [F, ~] = lsqnonneg(W, w, options);
        F_history(idx, :) = F';
    end
end
