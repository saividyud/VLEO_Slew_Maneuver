function openSimulationGUI
    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    evalin('base', 'clear'); % Clear the base workspace to ensure a clean state for the GUI
    clc
    close all

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
