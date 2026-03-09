%% Helper Function: Styling
% We use this to apply Times New Roman to buttons to avoid repeating code
function applyButtonStyle(btn)
    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    lightGrayColor = [211, 211, 211] / 256;
    
    btn.FontName = 'Times New Roman';
    btn.FontSize = 14;
    btn.FontWeight = 'bold';
    btn.BackgroundColor = lightGrayColor;
end
