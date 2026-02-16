%% Helper Function: Open New Simulation Window
function openNewSimulation(mainFigHandle)
    % Hide the Main Window
    mainFigHandle.Visible = 'off';

    % Create the New Window
    subFig = uifigure('Name', 'New Simulation - Unsaved', ...
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
        40,...
        40,...
        40,...
        40,...
        40,...
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
    lbl = uilabel(subGl, 'Text', 'New Simulation - Unsaved');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 24;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1, 5];
    lbl.HorizontalAlignment = 'center';

    %% Adding buttons
    % The Back Button (Lower Left)
    btnBack = uibutton(subGl, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 14;
    btnBack.Layout.Row = 10;     % Bottom Row
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

    % Simulation parameter button
    btnSimParams = uibutton(subGl, 'Text', 'Set Simulation Parameters');
    applyButtonStyle(btnSimParams)
    btnSimParams.Layout.Row = 4;
    btnSimParams.Layout.Column = 3;

    % Orbital parameter button
    btnOrbital = uibutton(subGl, 'Text', 'Set Orbital Parameters');
    applyButtonStyle(btnOrbital)
    btnOrbital.Layout.Row = 5;
    btnOrbital.Layout.Column = 3;

    % Attitude parameter button
    btnAttitude = uibutton(subGl, 'Text', 'Set Attitude Parameters');
    applyButtonStyle(btnAttitude)
    btnAttitude.Layout.Row = 6;
    btnAttitude.Layout.Column = 3;

    % Environmental parameter button
    btnEnvironmental = uibutton(subGl, 'Text', 'Set Environmental Parameters');
    applyButtonStyle(btnEnvironmental)
    btnEnvironmental.Layout.Row = 7;
    btnEnvironmental.Layout.Column = 3;

    % Controller parameter button
    btnController = uibutton(subGl, 'Text', 'Set Controller Parameters');
    applyButtonStyle(btnController)
    btnController.Layout.Row = 8;
    btnController.Layout.Column = 3;
    
    % Optional: Handle if user clicks the 'X' on window frame instead of Back
    subFig.CloseRequestFcn = @(~,~) goBack(subFig, mainFigHandle);

    %% Defining simulation structure in base workspace
    % This is where we initialize the structure that will hold all of our simulation parameters. We can update this structure as the user changes parameters in the GUI, and then use it to run the simulation when the user clicks "Run".
    simParams = struct();

    % Initializing all defaults
    simParams.saved = false; % Flag to track if data has been saved
    simParams.named = false; % Flag to track if simulation has been named
    simParams.fileName = 'My Simulation.mat'; % To store the file name if saved
    simParams.filePath = './simulations/My Simulation.mat'; % To store the file path if saved
    simParams.initParams = struct(); % This will hold all the actual simulation parameters (like time, orbital params, etc.) in a nested structure
    simParams.initParams.Simulation = struct(...
        'name', 'My Simulation', ...
        'type', 'Nonlinear', ...
        'initialTime', 0, ...
        'finalTime', 1000, ...
        'timeStep', 1, ...
        'relTol', 1e-6, ...
        'absTol', 1e-9 ...
    );

    earthRadius = 6378.14e3; % Earth equatorial radius [m]
    simParams.params.Orbit = struct(...
        'altitude', 400, ...
        'semiMajorAxis', 400 + earthRadius/1000, ... % Convert altitude to semimajor axis [km]
        'eccentricity', 0, ...
        'inclination', 0, ...
        'argPeriapse', 0, ...
        'RAAN', 0, ...
        'trueAnomaly', 0 ...
    );

    simParams.initParams.Attitude = struct(...
        'roll', 0, ...
        'pitch', 0, ...
        'yaw', 0, ...
        'rollRate', 0, ...
        'pitchRate', 0, ...
        'yawRate', 0, ...
        'Quaternion', [1, 0, 0, 0] ...
    );

    simParams.initParams.Environment = struct(...
        'secondsOfDay', 0, ...
        'dayOfYear', 1, ...
        'F107Average', 65, ...
        'F107Daily', 65, ...
        'magneticIndices', [3, 3, 3, 3, 3, 3, 3], ...
        'gasSurfaceInteractionModel', 'cook', ...
        'accommodationCoefficient', 1, ...
        'wallTemperature', 300, ...
        'specularReflectivity', 0.15, ...
        'diffuseReflectivity', 0.25, ...
        'enableAnomalousOxygen', false, ...
        'enableShadowAnalysis', true, ...
        'enableSolarRadiationPressure', true, ...
        'mu', 3.986004418e14 ... % Earth's gravitational parameter in m^3/s^2
    );

    simParams.initParams.Controller = struct(...
        'functionFile', 'myController.m', ...
        'Func', @(t, X) [0; 0; 0] ... % Default to a dummy controller that applies zero torque
    );

    guidata(subFig, simParams); % Store in figure's guidata for access in callbacks

    %% Adding logic to buttons
    % Define button actions
    btnSave.ButtonPushedFcn = @(~,~) saveSimulationData();
    btnRun.ButtonPushedFcn = @(~,~) runSimulation(subFig);
    btnResults.ButtonPushedFcn = @(~,~) AttitudeAnimatorGUI(subFig);
    btnSimParams.ButtonPushedFcn = @(~,~) functionalButtons(subFig, 'Set Simulation Parameters');
    btnOrbital.ButtonPushedFcn = @(~,~) functionalButtons(subFig, 'Set Orbital Parameters');
    btnAttitude.ButtonPushedFcn = @(~,~) functionalButtons(subFig, 'Set Attitude Parameters');
    btnEnvironmental.ButtonPushedFcn = @(~,~) functionalButtons(subFig, 'Set Environmental Parameters');
    btnController.ButtonPushedFcn = @(~,~) functionalButtons(subFig, 'Set Controller Parameters');

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

            % Checking if simulation name is passed in and updating title if so
            simName = simParams.initParams.Simulation.name;
            subFig.Name = ['New Simulation - ', simName]; % Update the window title to include the simulation name
            lbl.Text = simName; % Update title label text

        else
            % If not saved and not named, prompt for file name
            % Open the save file dialog, defaulting to a .mat file
            % Defaulting to simulation name if it exists, otherwise default to "My Simulation.mat"
            simulationsDirectory = './simulations/';
            defaultFileName = [simParams.initParams.Simulation.name, '.mat'];

            [file, path] = uiputfile([simulationsDirectory, defaultFileName], 'Save Simulation Parameters');

            % Check if the user selected a file or clicked Cancel
            if isequal(file, 0) || isequal(path, 0)
                disp('User clicked Cancel');
            else
                % Construct the full file path and save the data
                fullFileName = fullfile(simulationsDirectory, file);
                simParams.named = true; % Mark as named since we now have a file name
                simParams.fileName = file; % Store the file name in simParams
                simParams.saved = true; % Mark as saved
                simParams.filePath = fullFileName; % Store the file path in simParams
                save(fullFileName, 'simParams'); % Saves the variable 'simParams' to the .mat file
                disp(['Data saved to: ', fullFileName]);
                guidata(subFig, simParams); % Update the guidata with the new saved status

                % Checking if simulation name is passed in and updating title if so
                simName = simParams.initParams.Simulation.name;
                subFig.Name = ['New Simulation - ', simName]; % Update the window title to include the simulation name
                lbl.Text = simName; % Update title label text
            end
        end

        assignin('base', 'simParams', simParams);

    end

    % Helper Function: The Back Logic
    function goBack(currentFig, mainFig)
        % If there are unsaved changes, prompt the user to confirm before going back
        simParams = guidata(currentFig);
        if isfield(simParams, 'saved') && ~simParams.saved
            selection = uiconfirm(currentFig, 'You have unsaved changes. Are you sure you want to go back and lose these changes?', ...
                'Unsaved Changes', ...
                'Options', {'Yes', 'No'}, ...
                'DefaultOption', 2);
            if strcmp(selection, 'No')
                return; % User chose not to go back, so exit the function
            end
        end
        delete(currentFig);      % Destroy the sub-window
        mainFig.Visible = 'on';  % Reveal the main menu
    end

end