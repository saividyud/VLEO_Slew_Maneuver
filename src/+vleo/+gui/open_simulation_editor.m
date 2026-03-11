function open_simulation_editor(mainFigHandle, simParams, editorMode)
    simParams = vleo.gui.ensure_simulation_defaults(simParams);
    mainFigHandle.Visible = 'off';

    subFig = uifigure('Name', vleo.gui.editor_title(simParams, editorMode), ...
        'Position', [500, 300, 800, 450]);

    salmonColor = [1, 0.4941, 0.4392];
    yellowColor = [1, 0.87, 0.13];
    greenColor = [0, 191, 0] / 256;
    lightblueColor = [173, 216, 230] / 256;

    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, 30, '1x', 40, 40, 40, 40, 40, 40, 30};
    subGrid.ColumnWidth = {'fit', '1x', 300, '1x', 'fit'};

    titleLabel = uilabel(subGrid, 'Text', simParams.initParams.Simulation.name);
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 24;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 5];
    titleLabel.HorizontalAlignment = 'center';
    subFig.UserData = struct('editorMode', editorMode, 'titleLabel', titleLabel);

    btnBack = uibutton(subGrid, 'Text', 'Back');
    btnBack.FontName = 'Times New Roman';
    btnBack.FontSize = 14;
    btnBack.Layout.Row = 10;
    btnBack.Layout.Column = 1;
    btnBack.BackgroundColor = salmonColor;
    btnBack.ButtonPushedFcn = @(~, ~) vleo.gui.go_back(subFig, mainFigHandle);

    btnSave = uibutton(subGrid, 'Text', 'Save');
    btnSave.FontName = 'Times New Roman';
    btnSave.FontSize = 14;
    btnSave.BackgroundColor = greenColor;
    btnSave.Layout.Row = 1;
    btnSave.Layout.Column = 1;

    btnRun = uibutton(subGrid, 'Text', 'Run');
    btnRun.FontName = 'Times New Roman';
    btnRun.FontSize = 14;
    btnRun.BackgroundColor = yellowColor;
    btnRun.Layout.Row = 1;
    btnRun.Layout.Column = 5;

    btnResults = uibutton(subGrid, 'Text', 'Results');
    btnResults.FontName = 'Times New Roman';
    btnResults.FontSize = 14;
    btnResults.BackgroundColor = lightblueColor;
    btnResults.Layout.Row = 2;
    btnResults.Layout.Column = 5;

    btnSimParams = uibutton(subGrid, 'Text', 'Set Simulation Parameters');
    vleo.gui.apply_button_style(btnSimParams);
    btnSimParams.Layout.Row = 4;
    btnSimParams.Layout.Column = 3;

    btnOrbital = uibutton(subGrid, 'Text', 'Set Orbital Parameters');
    vleo.gui.apply_button_style(btnOrbital);
    btnOrbital.Layout.Row = 5;
    btnOrbital.Layout.Column = 3;

    btnAttitude = uibutton(subGrid, 'Text', 'Set Attitude Parameters');
    vleo.gui.apply_button_style(btnAttitude);
    btnAttitude.Layout.Row = 6;
    btnAttitude.Layout.Column = 3;

    btnEnvironmental = uibutton(subGrid, 'Text', 'Set Environmental Parameters');
    vleo.gui.apply_button_style(btnEnvironmental);
    btnEnvironmental.Layout.Row = 7;
    btnEnvironmental.Layout.Column = 3;

    btnController = uibutton(subGrid, 'Text', 'Set Controller Parameters');
    vleo.gui.apply_button_style(btnController);
    btnController.Layout.Row = 8;
    btnController.Layout.Column = 3;

    btnModes = uibutton(subGrid, 'Text', 'Set Aero/Control Modes');
    vleo.gui.apply_button_style(btnModes);
    btnModes.Layout.Row = 9;
    btnModes.Layout.Column = 3;

    subFig.CloseRequestFcn = @(~, ~) vleo.gui.go_back(subFig, mainFigHandle);
    guidata(subFig, simParams);

    btnSave.ButtonPushedFcn = @(~, ~) save_simulation_data();
    btnRun.ButtonPushedFcn = @(~, ~) vleo.gui.run_simulation(subFig);
    btnResults.ButtonPushedFcn = @(~, ~) vleo.gui.show_results(subFig);
    btnSimParams.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Simulation Parameters');
    btnOrbital.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Orbital Parameters');
    btnAttitude.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Attitude Parameters');
    btnEnvironmental.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Environmental Parameters');
    btnController.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Controller Parameters');
    btnModes.ButtonPushedFcn = @(~, ~) vleo.gui.dialogs.open_settings_dialog(subFig, 'Set Aero/Control Modes');

    function save_simulation_data()
        currentParams = guidata(subFig);
        simulationsDirectory = vleo.util.simulations_dir();

        if currentParams.named && currentParams.saved
            return;
        end

        if currentParams.named && ~currentParams.saved
            fullFileName = currentParams.filePath;
        else
            defaultFileName = [currentParams.initParams.Simulation.name, '.mat'];
            [fileName, pathName] = uiputfile(fullfile(simulationsDirectory, defaultFileName), 'Save Simulation Parameters');
            if isequal(fileName, 0) || isequal(pathName, 0)
                return;
            end
            fullFileName = fullfile(pathName, fileName);
            currentParams.named = true;
            currentParams.fileName = fileName;
            currentParams.filePath = fullFileName;
        end

        simParams = currentParams;
        simParams.saved = true;
        save(fullFileName, 'simParams');
        guidata(subFig, simParams);
        subFig.Name = vleo.gui.editor_title(simParams, editorMode);
        titleLabel.Text = simParams.initParams.Simulation.name;
        assignin('base', 'simParams', simParams);
    end
end
