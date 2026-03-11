function simParams = mark_simulation_dirty(subFigHandle)
    simParams = guidata(subFigHandle);
    simParams.saved = false;
    guidata(subFigHandle, simParams);

    editorMode = 'new';
    if isprop(subFigHandle, 'UserData') && isstruct(subFigHandle.UserData) && ...
            isfield(subFigHandle.UserData, 'editorMode')
        editorMode = subFigHandle.UserData.editorMode;
    end

    subFigHandle.Name = vleo.gui.editor_title(simParams, editorMode);

    if isprop(subFigHandle, 'UserData') && isstruct(subFigHandle.UserData) && ...
            isfield(subFigHandle.UserData, 'titleLabel') && isvalid(subFigHandle.UserData.titleLabel)
        subFigHandle.UserData.titleLabel.Text = simParams.initParams.Simulation.name;
    end
end
