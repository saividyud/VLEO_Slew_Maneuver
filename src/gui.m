function simulation_gui
    % 1. Create the Main Window (16:9 Aspect Ratio)
    % Position format: [left, bottom, width, height]
    % 800 width / 450 height = 1.77 (16:9)
    mainFig = uifigure('Name', 'Simulation Manager', ...
                       'Position', [500, 300, 800, 450]);

    % Define the Salmon Color (RGB triplet for #FF7E70)
    salmonColor = [1, 0.4941, 0.4392];

    % 2. Setup the Main Layout Manager
    % Columns:
    %   100  = Fixed space for Quit button
    %   '1x' = Flexible space (left side)
    %   200  = Fixed width for the center buttons (keeps them from stretching)
    %   '1x' = Flexible space (right side)
    %   100  = Fixed space
    % Rows: 
    %   'fit' = Title height
    %   '1x'  = Flexible space (top gap)
    %   'fit' = Button 1
    %   'fit' = Button 2
    %   'fit' = Button 3
    %   '1x'  = Flexible space (bottom gap)
    %   'fit' = Quit button row
    gl = uigridlayout(mainFig);
    gl.ColumnWidth = {100, '1x', 200, '1x', 100};
    gl.RowHeight = {'fit', '1x', 'fit', 'fit', 'fit', '1x', 'fit'};

    % 3. The Main Title
    lblTitle = uilabel(gl, 'Text', 'Simulation Dashboard');
    lblTitle.FontName = 'Times New Roman';
    lblTitle.FontSize = 24;
    lblTitle.FontWeight = 'bold';
    lblTitle.HorizontalAlignment = 'center';
    
    % Position: Row 1, Span all 5 columns
    lblTitle.Layout.Row = 1;
    lblTitle.Layout.Column = [1, 5]; 

    % 4. The Three Central Buttons
    % Button 1: New Simulation
    btnNew = uibutton(gl, 'Text', 'New Simulation');
    applyButtonStyle(btnNew);
    btnNew.Layout.Row = 3;
    btnNew.Layout.Column = 3; % Center column
    btnNew.ButtonPushedFcn = @(~,~) openNewSimulation(mainFig);

    % Button 2: Load Saved
    btnLoad = uibutton(gl, 'Text', 'Load Saved Simulation');
    applyButtonStyle(btnLoad);
    btnLoad.Layout.Row = 4;
    btnLoad.Layout.Column = 3; % Center column
    btnLoad.ButtonPushedFcn = @(~,~) openLoadSimulation(mainFig);

    % Button 3: Show Results
    btnRes = uibutton(gl, 'Text', 'Show Results');
    applyButtonStyle(btnRes);
    btnRes.Layout.Row = 5;
    btnRes.Layout.Column = 3; % Center column
    btnRes.ButtonPushedFcn = @(~,~) openResults(mainFig);

    % 5. The Quit Button (Lower Left)
    btnQuit = uibutton(gl, 'Text', 'Quit');
    applyButtonStyle(btnQuit);
    btnQuit.Layout.Row = 7;     % Bottom Row
    btnQuit.Layout.Column = 1;  % Left Column
    btnQuit.BackgroundColor = salmonColor;
    % Logic: Close the main figure entirely
    btnQuit.ButtonPushedFcn = @(~,~) delete(mainFig);
end

%% Helper Function: Styling
% We use this to apply Times New Roman to buttons to avoid repeating code
function applyButtonStyle(btn)
    btn.FontName = 'Times New Roman';
    btn.FontSize = 14;
    btn.FontWeight = 'bold';
end

%% Helper Function: Open New Simulation Window
function openNewSimulation(mainFigHandle)
    % 1. Hide the Main Window
    mainFigHandle.Visible = 'off';

    % 2. Create the New Window
    subFig = uifigure('Name', 'New Simulation', ...
                      'Position', [500, 300, 800, 450]);

    % Define the Salmon Color (RGB triplet for #FF7E70)
    salmonColor = [1, 0.4941, 0.4392];
    
    % 3. Create a Grid for the Sub-Window
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {'1x', 'fit'}; % Content area, then bottom row
    subGl.ColumnWidth = {'fit', '1x'}; % Left corner, rest of space

    % 4. Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', ['Current Mode: ', 'New Simulation']);
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 18;
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1, 2];
    lbl.HorizontalAlignment = 'center';

    % 5. The Back Button (Lower Left)
    btnBack = uibutton(subGl, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 12;
    btnBack.Layout.Row = 2;     % Bottom Row
    btnBack.Layout.Column = 1;  % Left Column
    btnBack.BackgroundColor = salmonColor;
    
    % 6. Link Back Button to Logic
    % Logic: Delete this window, turn Main Window back on
    btnBack.ButtonPushedFcn = @(~,~) goBack(subFig, mainFigHandle);
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);
end

%% Helper Function: Open Load Simulation Window
function openLoadSimulation(mainFigHandle)
    % 1. Hide the Main Window
    mainFigHandle.Visible = 'off';

    % 2. Create the New Window
    subFig = uifigure('Name', 'Load Simulation', ...
                      'Position', [500, 300, 800, 450]);

    % Define the Salmon Color (RGB triplet for #FF7E70)
    salmonColor = [1, 0.4941, 0.4392];
    
    % 3. Create a Grid for the Sub-Window
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {'1x', 'fit'}; % Content area, then bottom row
    subGl.ColumnWidth = {'fit', '1x'}; % Left corner, rest of space

    % 4. Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', ['Current Mode: ' 'Load Simulation']);
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 18;
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1, 2];
    lbl.HorizontalAlignment = 'center';

    % 5. The Back Button (Lower Left)
    btnBack = uibutton(subGl, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 12;
    btnBack.Layout.Row = 2;     % Bottom Row
    btnBack.Layout.Column = 1;  % Left Column
    btnBack.BackgroundColor = salmonColor;
    
    % 6. Link Back Button to Logic
    % Logic: Delete this window, turn Main Window back on
    btnBack.ButtonPushedFcn = @(~,~) goBack(subFig, mainFigHandle);
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);
end

%% Helper Function: Open Results Window
function openResults(mainFigHandle)
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
    btnBack.FontSize = 12;
    btnBack.Layout.Row = 2;     % Bottom Row
    btnBack.Layout.Column = 1;  % Left Column
    btnBack.BackgroundColor = salmonColor;
    
    % 6. Link Back Button to Logic
    % Logic: Delete this window, turn Main Window back on
    btnBack.ButtonPushedFcn = @(~,~) goBack(subFig, mainFigHandle);
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);
end

%% Helper Function: The Back Logic
function goBack(currentFig, mainFig)
    delete(currentFig);      % Destroy the sub-window
    mainFig.Visible = 'on';  % Reveal the main menu
end

simulation_gui