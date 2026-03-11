function go_back(currentFig, mainFig)
    simParams = guidata(currentFig);
    if isfield(simParams, 'saved') && ~simParams.saved
        selection = uiconfirm(currentFig, ...
            'You have unsaved changes. Are you sure you want to go back and lose these changes?', ...
            'Unsaved Changes', ...
            'Options', {'Yes', 'No'}, ...
            'DefaultOption', 2);
        if strcmp(selection, 'No')
            return;
        end
    end

    delete(currentFig);
    mainFig.Visible = 'on';
end
