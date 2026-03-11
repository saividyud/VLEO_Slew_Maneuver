function [controllerFunc, controllerName] = resolve_controller_handle(controllerSpec)
    controllerFunc = [];
    controllerName = 'vleo.control.attitude_pd_controller';

    if isa(controllerSpec, 'function_handle')
        controllerFunc = controllerSpec;
        controllerName = func2str(controllerFunc);
        return;
    end

    if isstruct(controllerSpec)
        if isfield(controllerSpec, 'Func') && isa(controllerSpec.Func, 'function_handle')
            controllerFunc = controllerSpec.Func;
            controllerName = func2str(controllerFunc);
            return;
        end

        if isfield(controllerSpec, 'functionName') && ~isempty(controllerSpec.functionName)
            controllerName = controllerSpec.functionName;
        elseif isfield(controllerSpec, 'functionFile') && ~isempty(controllerSpec.functionFile)
            controllerName = controllerSpec.functionFile;
        end
    elseif ~(isempty(controllerSpec) || (isstring(controllerSpec) && strlength(controllerSpec) == 0))
        controllerName = controllerSpec;
    end

    controllerName = vleo.control.normalize_controller_name(controllerName);
    controllerFunc = str2func(controllerName);
end
