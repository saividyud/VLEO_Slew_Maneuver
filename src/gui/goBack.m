%% Helper Function: The Back Logic
function goBack(currentFig, mainFig)
    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    delete(currentFig);      % Destroy the sub-window
    mainFig.Visible = 'on';  % Reveal the main menu
end
