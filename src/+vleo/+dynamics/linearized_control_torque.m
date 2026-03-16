function tauBody = linearized_control_torque(X, Xr, simParams)
% LINEARIZED_CONTROL_TORQUE Small-angle PD control torque about a fixed trim.

    tauBody = zeros(3, 1);
    if nargin < 3 || isempty(simParams)
        simParams = struct();
    end
    if isstruct(simParams) && ~vleo.util.get_mode_flag(simParams, 'enableControl', true)
        return;
    end

    qBodyFromEci = reshape(X(7:10), [], 1);
    qRef = reshape(Xr(1:4), [], 1);
    qBodyFromEci = qBodyFromEci / norm(qBodyFromEci);
    qRef = qRef / norm(qRef);

    omegaBody = reshape(X(11:13), [], 1);
    omegaRefBody = reshape(Xr(5:7), [], 1);
    attitudeErrorBody = attitude_error_body(qBodyFromEci, qRef);
    omegaErrorBody = omegaBody - omegaRefBody;

    [kpAttitude, kdRate] = default_linearized_gains();
    tauBody = -kpAttitude * attitudeErrorBody - kdRate * omegaErrorBody;
end

function attitudeErrorBody = attitude_error_body(qBodyFromEci, qRef)
    rBodyFromEci = quat2dcm(qBodyFromEci');
    rRefFromEci = quat2dcm(qRef');
    rError = rBodyFromEci * rRefFromEci';
    attitudeErrorBody = 0.5 * [ ...
        rError(2, 3) - rError(3, 2); ...
        rError(3, 1) - rError(1, 3); ...
        rError(1, 2) - rError(2, 1) ...
    ];
end

function [kpAttitude, kdRate] = default_linearized_gains()
    kpAttitude = 1e-3 * eye(3);
    kdRate = 0.1 * eye(3);
end
