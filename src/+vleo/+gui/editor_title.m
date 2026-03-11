function titleText = editor_title(simParams, editorMode)
    simName = simParams.initParams.Simulation.name;
    switch lower(char(editorMode))
        case 'load'
            prefix = 'Load Simulation - ';
        otherwise
            prefix = 'New Simulation - ';
    end

    titleText = [prefix, simName];
    if isfield(simParams, 'saved') && ~simParams.saved
        titleText = [titleText, ' *'];
    end
end
