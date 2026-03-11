function value = get_env_value(simParams, fieldName, defaultValue)
    value = defaultValue;

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Environment')
        return;
    end

    if isfield(simParams.initParams.Environment, fieldName)
        value = simParams.initParams.Environment.(fieldName);
    end
end
