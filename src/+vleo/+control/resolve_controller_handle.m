function [controllerFunc, controllerName] = resolve_controller_handle(controllerSpec)
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

    controllerName = normalize_controller_name(controllerName);
    controllerFunc = str2func(controllerName);
end

function controllerName = normalize_controller_name(controllerName)
    controllerName = strtrim(char(controllerName));

    if isempty(controllerName)
        controllerName = 'vleo.control.attitude_pd_controller';
        return;
    end

    [~, fileName, extension] = fileparts(controllerName);
    aliasName = resolve_alias(controllerName, [fileName, extension]);
    if ~isempty(aliasName)
        controllerName = aliasName;
        return;
    end

    if contains(controllerName, '+') || contains(controllerName, filesep)
        parts = regexp(controllerName, '[\\/]+', 'split');
        parts = parts(~cellfun('isempty', parts));
        packageParts = {};
        for i = 1:numel(parts)
            part = parts{i};
            if startsWith(part, '+')
                packageParts{end + 1} = part(2:end); %#ok<AGROW>
            elseif i == numel(parts)
                [~, finalName] = fileparts(part);
                packageParts{end + 1} = finalName; %#ok<AGROW>
            end
        end
        if ~isempty(packageParts)
            controllerName = strjoin(packageParts, '.');
            return;
        end
    end

    if endsWith(controllerName, '.m')
        controllerName = controllerName(1:end - 2);
    end
end

function aliasName = resolve_alias(varargin)
    aliasName = '';
    aliasPairs = { ...
        'control_torques', 'vleo.control.attitude_pd_controller'; ...
        'control_torques.m', 'vleo.control.attitude_pd_controller'; ...
        'myController', 'vleo.control.zero_torque_controller'; ...
        'myController.m', 'vleo.control.zero_torque_controller'};

    for i = 1:nargin
        candidate = varargin{i};
        if isempty(candidate)
            continue;
        end
        for j = 1:size(aliasPairs, 1)
            if strcmp(candidate, aliasPairs{j, 1})
                aliasName = aliasPairs{j, 2};
                return;
            end
        end
    end
end
