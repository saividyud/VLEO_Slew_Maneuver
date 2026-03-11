function open_load_simulation(mainFigHandle)
    simulationsDirectory = vleo.util.simulations_dir();
    [fileName, pathName] = uigetfile(fullfile(simulationsDirectory, '*.mat'), 'Select a simulation file to load');
    if isequal(fileName, 0) || isequal(pathName, 0)
        return;
    end

    loadedData = load(fullfile(pathName, fileName), 'simParams');
    if ~isfield(loadedData, 'simParams')
        uialert(mainFigHandle, 'The selected file does not contain simParams.', 'Load Error', 'Icon', 'error');
        return;
    end

    simParams = vleo.gui.ensure_simulation_defaults(loadedData.simParams);
    simParams.fileName = fileName;
    simParams.filePath = fullfile(pathName, fileName);
    simParams.named = true;
    vleo.gui.open_simulation_editor(mainFigHandle, simParams, 'load');
end
