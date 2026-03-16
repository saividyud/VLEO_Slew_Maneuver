function Xd = sat_dynamics_controlled(t, X, targetPointingEci, params, varargin)
    parser = inputParser;
    addParameter(parser, 'useJ2', true, @islogical);
    addParameter(parser, 'useAtmDrag', true, @islogical);
    addParameter(parser, 'useControl', true, @islogical);
    parse(parser, varargin{:});
    qNorm = norm(X(7:10));
    if qNorm < 1e-10 || ~isfinite(qNorm)
        Xd = zeros(13, 1);
        return;
    end

    X(7:10) = X(7:10) / qNorm;
    qEciToBody = X(7:10);
    rEciToBody = quat2dcm(qEciToBody');
    rBodyToEci = rEciToBody';

    tauBody = zeros(3, 1);
    if parser.Results.useControl
        cameraBody = [0; 0; 1];
        currentPointingEci = rBodyToEci * cameraBody;
        currentPointingEci = currentPointingEci / norm(currentPointingEci);
        targetPointingEci = targetPointingEci / norm(targetPointingEci);
        pointingErrorEci = cross(currentPointingEci, targetPointingEci);
        pointingErrorBody = rEciToBody * pointingErrorEci;
        tauBody = params.Kp_att * pointingErrorBody;
    end

    aeroConfig = [];
    if parser.Results.useAtmDrag && isfield(params, 'aero') && ~isempty(params.aero)
        aeroConfig = params.aero;
        aeroConfig.mass = params.mass;
    end

    Xd = vleo.dynamics.sat_dynamics_base(t, X, params.mu, aeroConfig, tauBody, params.I_CB, parser.Results.useJ2);
end
