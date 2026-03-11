function tau = evaluate_controller_torque(simParams, t, X)
    tau = [0; 0; 0];

    if ~vleo.util.get_mode_flag(simParams, 'enableControl', true)
        return;
    end

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Controller')
        return;
    end

    try
        [controllerFunc, ~] = vleo.control.resolve_controller_handle(simParams.initParams.Controller);
        tauCandidate = reshape(controllerFunc(t, X), [], 1);
        if numel(tauCandidate) == 3 && all(isfinite(tauCandidate))
            tau = tauCandidate;
        end
    catch
        tau = [0; 0; 0];
    end
end
