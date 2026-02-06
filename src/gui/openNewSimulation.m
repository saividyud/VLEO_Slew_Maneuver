%% Helper Function: Open New Simulation Window
function openNewSimulation(mainFigHandle)
    % Hide the Main Window
    mainFigHandle.Visible = 'off';

    % Create the New Window
    subFig = uifigure('Name', 'New Simulation', ...
                      'Position', [500, 300, 800, 450]);

    % Define colors
    salmonColor = [1, 0.4941, 0.4392];
    yellowColor = [1, 0.87, 0.13];
    greenColor = [0, 191, 0] / 256;
    lightblueColor = [173, 216, 230] / 256;
    
    % Create a Grid for the Sub-Window
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        30,...
        '1x',...
        50,...
        50,...
        50,...
        50,...
        '1x',...
        30,...
    };
    subGl.ColumnWidth = {
        'fit',...
        '1x',...
        300,...
        '1x',...
        'fit',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'New Simulation');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 24;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = 3;
    lbl.HorizontalAlignment = 'center';

    %% Adding buttons
    % The Back Button (Lower Left)
    btnBack = uibutton(subGl, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 14;
    btnBack.Layout.Row = 9;     % Bottom Row
    btnBack.Layout.Column = 1;  % Left Column
    btnBack.BackgroundColor = salmonColor;

    % Logic: Delete this window, turn Main Window back on
    btnBack.ButtonPushedFcn = @(~,~) goBack(subFig, mainFigHandle);

    % Save button (upper left)
    btnSave = uibutton(subGl, 'Text', 'Save');
    btnSave.FontName = 'Times New Roman';
    btnSave.FontSize = 14;
    btnSave.BackgroundColor = greenColor;
    btnSave.Layout.Row = 1; % Top row
    btnSave.Layout.Column = 1; % Left column

    % Run button (upper right)
    btnRun = uibutton(subGl, 'Text', 'Run');
    btnRun.FontName = 'Times New Roman';
    btnRun.FontSize = 14;
    btnRun.BackgroundColor = yellowColor;
    btnRun.Layout.Row = 1; % Top row
    btnRun.Layout.Column = 5; % Left column

    % Results button (upper right)
    btnResults = uibutton(subGl, 'Text', 'Results');
    btnResults.FontName = 'Times New Roman';
    btnResults.FontSize = 14;
    btnResults.BackgroundColor = lightblueColor;
    btnResults.Layout.Row = 2;
    btnResults.Layout.Column = 5;

    % Orbital parameter button
    btnOrbital = uibutton(subGl, 'Text', 'Set Orbital Parameters');
    applyButtonStyle(btnOrbital)
    btnOrbital.Layout.Row = 4;
    btnOrbital.Layout.Column = 3;

    % Attitude parameter button
    btnAttitude = uibutton(subGl, 'Text', 'Set Attitude Parameters');
    applyButtonStyle(btnAttitude)
    btnAttitude.Layout.Row = 5;
    btnAttitude.Layout.Column = 3;

    % Environmental parameter button
    btnEnvironmental = uibutton(subGl, 'Text', 'Set Environmental Parameters');
    applyButtonStyle(btnEnvironmental)
    btnEnvironmental.Layout.Row = 6;
    btnEnvironmental.Layout.Column = 3;

    % Controller parameter button
    btnController = uibutton(subGl, 'Text', 'Set Controller Parameters');
    applyButtonStyle(btnController)
    btnController.Layout.Row = 7;
    btnController.Layout.Column = 3;
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);

    %% Defining simulation structure in base workspace
    % This is where we initialize the structure that will hold all of our simulation parameters. We can update this structure as the user changes parameters in the GUI, and then use it to run the simulation when the user clicks "Run".
    simParams = struct();
    simParams.saved = false; % Flag to track if data has been saved
    simParams.named = false; % Flag to track if simulation has been named
    guidata(subFig, simParams); % Store in figure's guidata for access in callbacks

    %% Adding logic to buttons
    % Define button actions
    btnSave.ButtonPushedFcn = @(~,~) saveSimulationData();
    btnRun.ButtonPushedFcn = @(~,~) runSimulation(subFig);
    btnResults.ButtonPushedFcn = @(~,~) displayResults(subFig);
    btnOrbital.ButtonPushedFcn = @(~,~) setOrbitalParameters(subFig);
    btnAttitude.ButtonPushedFcn = @(~,~) setAttitudeParameters(subFig);
    btnEnvironmental.ButtonPushedFcn = @(~,~) setEnvironmentalParameters(subFig);
    btnController.ButtonPushedFcn = @(~,~) setControllerParameters(subFig);

    function saveSimulationData()
        disp('Saving simulation data...');

        simParams = guidata(subFig);

        if simParams.named && simParams.saved
            disp('Simulation data already saved.');

        elseif simParams.named && ~simParams.saved
            % If already named but not saved, overwrite existing file
            fullFileName = simParams.filePath; % Assuming simParams has a field for file path
            save(fullFileName, 'simParams');
            disp(['Data overwritten to: ', fullFileName]);
            simParams.saved = true; % Mark as saved
            guidata(subFig, simParams); % Update the guidata with the new saved status

        else
            % If not saved and not named, prompt for file name
            % Open the save file dialog, defaulting to a .mat file
            [file, path] = uiputfile('*.mat', 'Save Simulation Parameters');

            % Check if the user selected a file or clicked Cancel
            if isequal(file, 0) || isequal(path, 0)
                disp('User clicked Cancel');
            else
                % Construct the full file path and save the data
                fullFileName = fullfile(path, file);
                simParams.named = true; % Mark as named since we now have a file name
                simParams.fileName = file; % Store the file name in simParams
                simParams.saved = true; % Mark as saved
                simParams.filePath = fullFileName; % Store the file path in simParams
                save(fullFileName, 'simParams'); % Saves the variable 'simParams' to the .mat file
                disp(['Data saved to: ', fullFileName]);
                guidata(subFig, simParams); % Update the guidata with the new saved status

                subFig.Name = ['New Simulation - ', file]; % Update the window title to include the file name
                lbl.Text = ['New Simulation - ', file]; % Update title label text
            end
        end

        assignin('base', 'simParams', simParams);

    end
end

function setOrbitalParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];

    % Create the New Window
    subFig = uifigure('Name', 'Set Orbital Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        'fit',...
        100,...
        75,...
        'fit',...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Orbital Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 5];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining entries
    % Semimajor Axis
    lblSMA = uilabel(subGl, 'Text', 'Semi-Major Axis (a) [m]');
    lblSMA.FontName = 'Times New Roman';
    lblSMA.FontSize = 14;
    lblSMA.Layout.Row = 3;
    lblSMA.Layout.Column = 4;
    entrySMA = uieditfield(subGl, 'numeric');
    entrySMA.FontName = 'Times New Roman';
    entrySMA.FontSize = 14;
    entrySMA.Layout.Row = 3;
    entrySMA.Layout.Column = 3;
    entrySMA.Value = 400e3; % Default to 400 km

    % Eccentricity
    lblEcc = uilabel(subGl, 'Text', 'Eccentricity (e)');
    lblEcc.FontName = 'Times New Roman';
    lblEcc.FontSize = 14;
    lblEcc.Layout.Row = 4;
    lblEcc.Layout.Column = 4;
    entryEcc = uieditfield(subGl, 'numeric');
    entryEcc.FontName = 'Times New Roman';
    entryEcc.FontSize = 14;
    entryEcc.Layout.Row = 4;
    entryEcc.Layout.Column = 3;
    entryEcc.Value = 0; % Default to 0

    % Inclination
    lblInc = uilabel(subGl, 'Text', 'Inclination (i) [deg]');
    lblInc.FontName = 'Times New Roman';
    lblInc.FontSize = 14;
    lblInc.Layout.Row = 5;
    lblInc.Layout.Column = 4;
    entryInc = uieditfield(subGl, 'numeric');
    entryInc.FontName = 'Times New Roman';
    entryInc.FontSize = 14;
    entryInc.Layout.Row = 5;
    entryInc.Layout.Column = 3;
    entryInc.Value = 0; % Default to 0

    % Argument of periapse
    lblAOP = uilabel(subGl, 'Text', 'Argument of periapse (ω) [deg]');
    lblAOP.FontName = 'Times New Roman';
    lblAOP.FontSize = 14;
    lblAOP.Layout.Row = 6;
    lblAOP.Layout.Column = 4;
    entryAOP = uieditfield(subGl, 'numeric');
    entryAOP.FontName = 'Times New Roman';
    entryAOP.FontSize = 14;
    entryAOP.Layout.Row = 6;
    entryAOP.Layout.Column = 3;
    entryAOP.Value = 0; % Default to 0

    % Right Ascension of Ascending Node
    lblRAAN = uilabel(subGl, 'Text', 'Right ascension of ascending node (Ω) [deg]');
    lblRAAN.FontName = 'Times New Roman';
    lblRAAN.FontSize = 14;
    lblRAAN.Layout.Row = 7;
    lblRAAN.Layout.Column = 4;
    entryRAAN = uieditfield(subGl, 'numeric');
    entryRAAN.FontName = 'Times New Roman';
    entryRAAN.FontSize = 14;
    entryRAAN.Layout.Row = 7;
    entryRAAN.Layout.Column = 3;
    entryRAAN.Value = 0; % Default to 0

    % True anomaly
    lblTA = uilabel(subGl, 'Text', 'True anomaly (ν) [deg]');
    lblTA.FontName = 'Times New Roman';
    lblTA.FontSize = 14;
    lblTA.Layout.Row = 8;
    lblTA.Layout.Column = 4;
    entryTA = uieditfield(subGl, 'numeric');
    entryTA.FontName = 'Times New Roman';
    entryTA.FontSize = 14;
    entryTA.Layout.Row = 8;
    entryTA.Layout.Column = 3;
    entryTA.Value = 0; % Default to 0

    %% Saving logic
    entryStruct = struct(...
        'semiMajorAxis', entrySMA, ...
        'eccentricity', entryEcc, ...
        'inclination', entryInc, ...
        'argPeriapse', entryAOP, ...
        'RAAN', entryRAAN, ...
        'trueAnomaly', entryTA ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose();

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.params.Orbital.semiMajorAxis = entryStruct.semiMajorAxis.Value;
        simParams.params.Orbital.eccentricity = entryStruct.eccentricity.Value;
        simParams.params.Orbital.inclination = entryStruct.inclination.Value;
        simParams.params.Orbital.argPeriapse = entryStruct.argPeriapse.Value;
        simParams.params.Orbital.RAAN = entryStruct.RAAN.Value;
        simParams.params.Orbital.trueAnomaly = entryStruct.trueAnomaly.Value;

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end
end

function setAttitudeParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    
    % Create the New Window
    subFig = uifigure('Name', 'Set Attitude Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        'fit',...
        150,...
        75,...
        'fit',...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Attitude Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 5];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining entries
    % Roll
    lblSMA = uilabel(subGl, 'Text', 'Roll ($\phi$) [deg]', 'Interpreter', 'latex');
    lblSMA.FontName = 'Times New Roman';
    lblSMA.FontSize = 14;
    lblSMA.Layout.Row = 3;
    lblSMA.Layout.Column = 4;
    entrySMA = uieditfield(subGl, 'numeric');
    entrySMA.FontName = 'Times New Roman';
    entrySMA.FontSize = 14;
    entrySMA.Layout.Row = 3;
    entrySMA.Layout.Column = 3;
    entrySMA.Value = 0; % Default to 0

    % Pitch
    lblEcc = uilabel(subGl, 'Text', 'Pitch ($\theta$) [deg]', 'Interpreter', 'latex');
    lblEcc.FontName = 'Times New Roman';
    lblEcc.FontSize = 14;
    lblEcc.Layout.Row = 4;
    lblEcc.Layout.Column = 4;
    entryEcc = uieditfield(subGl, 'numeric');
    entryEcc.FontName = 'Times New Roman';
    entryEcc.FontSize = 14;
    entryEcc.Layout.Row = 4;
    entryEcc.Layout.Column = 3;
    entryEcc.Value = 0; % Default to 0

    % Yaw
    lblInc = uilabel(subGl, 'Text', 'Yaw ($\psi$) [deg]', 'Interpreter', 'latex');
    lblInc.FontName = 'Times New Roman';
    lblInc.FontSize = 14;
    lblInc.Layout.Row = 5;
    lblInc.Layout.Column = 4;
    entryInc = uieditfield(subGl, 'numeric');
    entryInc.FontName = 'Times New Roman';
    entryInc.FontSize = 14;
    entryInc.Layout.Row = 5;
    entryInc.Layout.Column = 3;
    entryInc.Value = 0; % Default to 0

    % Roll rate
    lblAOP = uilabel(subGl, 'Text', 'Roll rate ($\dot{\phi}$) [deg/s]', 'Interpreter', 'latex');
    lblAOP.FontName = 'Times New Roman';
    lblAOP.FontSize = 14;
    lblAOP.Layout.Row = 6;
    lblAOP.Layout.Column = 4;
    entryAOP = uieditfield(subGl, 'numeric');
    entryAOP.FontName = 'Times New Roman';
    entryAOP.FontSize = 14;
    entryAOP.Layout.Row = 6;
    entryAOP.Layout.Column = 3;
    entryAOP.Value = 0; % Default to 0

    % Pitch rate
    lblRAAN = uilabel(subGl, 'Text', 'Pitch rate ($\dot{\theta}$) [deg/s]', 'Interpreter', 'latex');
    lblRAAN.FontName = 'Times New Roman';
    lblRAAN.FontSize = 14;
    lblRAAN.Layout.Row = 7;
    lblRAAN.Layout.Column = 4;
    entryRAAN = uieditfield(subGl, 'numeric');
    entryRAAN.FontName = 'Times New Roman';
    entryRAAN.FontSize = 14;
    entryRAAN.Layout.Row = 7;
    entryRAAN.Layout.Column = 3;
    entryRAAN.Value = 0; % Default to 0

    % Yaw rate
    lblTA = uilabel(subGl, 'Text', 'Yaw rate ($\dot{\psi}$) [deg/s]', 'Interpreter', 'latex');
    lblTA.FontName = 'Times New Roman';
    lblTA.FontSize = 14;
    lblTA.Layout.Row = 8;
    lblTA.Layout.Column = 4;
    entryTA = uieditfield(subGl, 'numeric');
    entryTA.FontName = 'Times New Roman';
    entryTA.FontSize = 14;
    entryTA.Layout.Row = 8;
    entryTA.Layout.Column = 3;
    entryTA.Value = 0; % Default to 0

    %% Saving logic
    entryStruct = struct(...
        'roll', entrySMA, ...
        'pitch', entryEcc, ...
        'yaw', entryInc, ...
        'rollRate', entryAOP, ...
        'pitchRate', entryRAAN, ...
        'yawRate', entryTA ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose;

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.params.Attitude.roll = entryStruct.roll.Value;
        simParams.params.Attitude.pitch = entryStruct.pitch.Value;
        simParams.params.Attitude.yaw = entryStruct.yaw.Value;
        simParams.params.Attitude.rollRate = entryStruct.rollRate.Value;
        simParams.params.Attitude.pitchRate = entryStruct.pitchRate.Value;
        simParams.params.Attitude.yawRate = entryStruct.yawRate.Value;

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end
end