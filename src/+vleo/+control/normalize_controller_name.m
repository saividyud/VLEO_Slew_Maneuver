function controllerName = normalize_controller_name(controllerName)
    controllerName = strtrim(char(controllerName));

    if isempty(controllerName)
        controllerName = 'vleo.control.attitude_pd_controller';
        return;
    end

    [~, fileName, extension] = fileparts(controllerName);
    aliases = vleo.control.controller_name_aliases();

    if isempty(extension) && aliases.isKey(controllerName)
        controllerName = aliases(controllerName);
        return;
    end

    keyWithExtension = [fileName, extension];
    if aliases.isKey(keyWithExtension)
        controllerName = aliases(keyWithExtension);
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
