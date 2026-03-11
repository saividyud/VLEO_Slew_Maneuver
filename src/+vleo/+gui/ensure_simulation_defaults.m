function simParams = ensure_simulation_defaults(simParams)
    defaults = vleo.gui.default_simulation_params();

    simParams = fill_missing_fields(simParams, defaults, {'saved', 'named', 'fileName', 'filePath'});

    if ~isfield(simParams, 'initParams') || ~isstruct(simParams.initParams)
        simParams.initParams = defaults.initParams;
        return;
    end

    simParams.initParams = fill_missing_fields(simParams.initParams, defaults.initParams, ...
        {'Simulation', 'Orbit', 'Attitude', 'Environment', 'Controller', 'Modes'});
    simParams.initParams.Simulation = fill_missing_fields(simParams.initParams.Simulation, defaults.initParams.Simulation);
    simParams.initParams.Orbit = fill_missing_fields(simParams.initParams.Orbit, defaults.initParams.Orbit);
    simParams.initParams.Attitude = fill_missing_fields(simParams.initParams.Attitude, defaults.initParams.Attitude);
    simParams.initParams.Environment = fill_missing_fields(simParams.initParams.Environment, defaults.initParams.Environment);
    simParams.initParams.Modes = fill_missing_fields(simParams.initParams.Modes, defaults.initParams.Modes);
    simParams.initParams.Controller = fill_missing_fields(simParams.initParams.Controller, defaults.initParams.Controller);

    if ~isfield(simParams.initParams.Environment, 'gasSurfaceInteractionModel') && ...
            isfield(simParams.initParams.Environment, 'gasSurfaceModel')
        simParams.initParams.Environment.gasSurfaceInteractionModel = simParams.initParams.Environment.gasSurfaceModel;
    end

    simParams.initParams.Environment.gasSurfaceInteractionModel = ...
        vleo.util.normalize_gsi_model(simParams.initParams.Environment.gasSurfaceInteractionModel);

    [controllerFunc, controllerName] = vleo.control.resolve_controller_handle(simParams.initParams.Controller);
    simParams.initParams.Controller.Func = controllerFunc;
    simParams.initParams.Controller.functionName = controllerName;
    simParams.initParams.Controller.functionFile = controllerName;
end

function target = fill_missing_fields(target, defaults, fieldNames)
    if nargin < 3
        fieldNames = fieldnames(defaults);
    end

    if isempty(target)
        target = defaults;
        return;
    end

    for i = 1:numel(fieldNames)
        fieldName = fieldNames{i};
        if ~isfield(target, fieldName)
            target.(fieldName) = defaults.(fieldName);
        end
    end
end
