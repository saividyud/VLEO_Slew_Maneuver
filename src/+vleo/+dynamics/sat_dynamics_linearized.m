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
<<<<<<< Updated upstream
=======

    q0 = Xr(7);
    q1 = Xr(8);
    q2 = Xr(9);
    q3 = Xr(10);
    wx = Xr(11);
    wy = Xr(12);
    wz = Xr(13);

    mu = 3.986e14;
    r = norm(X(1:3));
    m = 12;
    c = .3;
    b = .2; 
    a = .1;
    inertiaBody = 1 / 12 * m * [ (b^2 + c^2) .01 .01 ; .01 (a^2 + c^2) .01 ; .01 .01 (a^2 + b^2)];

>>>>>>> Stashed changes
    Xd(1:3) = X(4:6);

    rNorm = norm(X(1:3));
    Xd(4:6) = -mu * X(1:3) / rNorm^3;

    omegaBody = X(11:13);
    omegaRefBody = Xr(5:7);
    Xd(7:10) = quaternion_kinematics(X(7:10), omegaBody);

<<<<<<< Updated upstream
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
=======
    P = [.07 0 0; 0 .07 0; 0 0 .07];
    Kp = [.0015 0 0 ; 0 .0015 0 ; 0 0 .0015];
    delw = X(11:13);
    delx = Ref(t) - X(8:10);
    u = -Kp *delx;
    u = u- P * delw;

    bMatrix = [zeros(4, 3); inv(inertiaBody)];
    Xd(7:13) = aMatrix * X(7:13) + bMatrix * u;
>>>>>>> Stashed changes
end
