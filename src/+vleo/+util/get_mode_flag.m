function enabled = get_mode_flag(simParams, fieldName, defaultValue)
    enabled = defaultValue;

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Modes')
        return;
    end

    if isfield(simParams.initParams.Modes, fieldName)
        enabled = logical(simParams.initParams.Modes.(fieldName));
    end
end
