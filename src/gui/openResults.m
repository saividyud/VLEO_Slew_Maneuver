%% Helper Function: Open Results Window
function openResults(mainFigHandle)
    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    % 1. Hide the Main Window
    mainFigHandle.Visible = 'off';

    % 2. Create the New Window
    subFig = uifigure('Name', 'Results Viewer', ...
                      'Position', [500, 300, 800, 450]);

    % Define the Salmon Color (RGB triplet for #FF7E70)
    salmonColor = [1, 0.4941, 0.4392];
    
    % 3. Create a Grid for the Sub-Window
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {'1x', 'fit'}; % Content area, then bottom row
    subGl.ColumnWidth = {'fit', '1x'}; % Left corner, rest of space

    % 4. Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', ['Current Mode: ' 'Results Viewer']);
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 18;
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1, 2];
    lbl.HorizontalAlignment = 'center';

    % 5. The Back Button (Lower Left)
    btnBack = uibutton(subGl, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 14;
    btnBack.Layout.Row = 2;     % Bottom Row
    btnBack.Layout.Column = 1;  % Left Column
    btnBack.BackgroundColor = salmonColor;
    
    % 6. Link Back Button to Logic
    % Logic: Delete this window, turn Main Window back on
    btnBack.ButtonPushedFcn = @(~,~) goBack(subFig, mainFigHandle);
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);
end
