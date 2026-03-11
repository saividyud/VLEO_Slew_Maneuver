function Xd = sat_dynamics_controlled(t, X, targetPointingEci, params, varargin)
    parser = inputParser;
    addParameter(parser, 'useJ2', true, @islogical);
    addParameter(parser, 'useAtmDrag', true, @islogical);
    addParameter(parser, 'useControl', true, @islogical);
    addParameter(parser, 'useSRP', false, @islogical);
    parse(parser, varargin{:});
    flags = parser.Results;

    Xd = zeros(13, 1);
    mu = params.mu;
    inertiaBody = params.I_CB;
    kpAttitude = params.Kp_att;

    rEci = X(1:3);
    vEci = X(4:6);
    qEciToBody = X(7:10);
    omegaEci = X(11:13);

    qEciToBody = qEciToBody / norm(qEciToBody);
    rEciToBody = quat2dcm(qEciToBody');
    rBodyToEci = rEciToBody';
    omegaBody = rEciToBody * omegaEci;

    Xd(1:3) = vEci;
    rNorm = norm(rEci);
    a2body = -mu * rEci / rNorm^3;
    aTotal = a2body;
    if flags.useJ2 %#ok<NASGU>
    end
    if flags.useAtmDrag %#ok<NASGU>
    end
    if flags.useSRP %#ok<NASGU>
    end
    Xd(4:6) = aTotal;

    omegaMatrix = [0, -omegaBody(1), -omegaBody(2), -omegaBody(3); ...
        omegaBody(1), 0, omegaBody(3), -omegaBody(2); ...
        omegaBody(2), -omegaBody(3), 0, omegaBody(1); ...
        omegaBody(3), omegaBody(2), -omegaBody(1), 0];
    Xd(7:10) = 0.5 * omegaMatrix * qEciToBody;

    tauBody = zeros(3, 1);
    if flags.useControl
        cameraBody = [0; 0; 1];
        currentPointingEci = rBodyToEci * cameraBody;
        currentPointingEci = currentPointingEci / norm(currentPointingEci);
        targetPointingEci = targetPointingEci / norm(targetPointingEci);
        pointingErrorEci = cross(currentPointingEci, targetPointingEci);
        pointingErrorBody = rEciToBody * pointingErrorEci;
        tauBody = kpAttitude * pointingErrorBody;
        params.tau_eci = rBodyToEci * tauBody; %#ok<NASGU>
    end

    alphaBody = inertiaBody \ (tauBody - cross(omegaBody, inertiaBody * omegaBody));
    Xd(11:13) = rBodyToEci * alphaBody;
end
