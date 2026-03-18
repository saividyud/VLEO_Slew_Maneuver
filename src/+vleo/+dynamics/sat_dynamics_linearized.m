function Xd = sat_dynamics_linearized(~, X, Xr, simParams)
% SAT_DYNAMICS_LINEARIZED Small-angle attitude model about a fixed reference trim.

    if nargin < 3 || numel(Xr) ~= 7
        error('VLEO_Slew_Maneuver:InvalidLinearReference', ...
            'Xr must be a 7x1 reference state [q_ref; omega_ref_body].');
    end
    if nargin < 4 || isempty(simParams)
        simParams = struct();
    end

    c = vleo.util.constants();
    mu = vleo.util.get_env_value(simParams, 'mu', c.mu_earth);
    inertiaBody = vleo.dynamics.default_vehicle_config().inertiaBody;

    qBodyFromEci = reshape(X(7:10), [], 1);
    qNorm = norm(qBodyFromEci);
    if qNorm < 1e-10 || ~isfinite(qNorm)
        Xd = zeros(13, 1);
        return;
    end
    X(7:10) = qBodyFromEci / qNorm;

    Xd = zeros(13, 1);

    Xd(1:3) = X(4:6);

    rNorm = norm(X(1:3));
    Xd(4:6) = -mu * X(1:3) / rNorm^3;

    omegaBody = X(11:13);
    omegaRefBody = Xr(5:7);
    Xd(7:10) = quaternion_kinematics(X(7:10), omegaBody);

    tauBody = vleo.dynamics.linearized_control_torque(X, Xr, simParams);
    omegaErrorBody = omegaBody - omegaRefBody;
    Xd(11:13) = inertiaBody \ (tauBody ...
        - cross(omegaErrorBody, inertiaBody * omegaRefBody) ...
        - cross(omegaRefBody, inertiaBody * omegaErrorBody));
end

function qDot = quaternion_kinematics(qScalarFirst, omegaBody)
    bMatrix = [qScalarFirst(1), -qScalarFirst(2), -qScalarFirst(3), -qScalarFirst(4); ...
        qScalarFirst(2), qScalarFirst(1), -qScalarFirst(4), qScalarFirst(3); ...
        qScalarFirst(3), qScalarFirst(4), qScalarFirst(1), -qScalarFirst(2); ...
        qScalarFirst(4), -qScalarFirst(3), qScalarFirst(2), qScalarFirst(1)];
    qDot = 0.5 * bMatrix * [0; omegaBody];
end
