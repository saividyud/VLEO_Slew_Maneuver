function set_simulation_parameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    simParams = guidata(subFigHandle);

    subFig = uifigure('Name', 'Set Simulation Parameters', 'Position', [500, 300, 800, 450]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 30, 30, 30, 30, 30, 30, 30, '1x', 30};
    subGrid.ColumnWidth = {'fit', 50, 175, 'fit', '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Simulation Parameters');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 20;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 5];
    titleLabel.HorizontalAlignment = 'center';

    btnSaveClose = uibutton(subGrid, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 11;
    btnSaveClose.Layout.Column = 1;
    btnSaveClose.BackgroundColor = salmonColor;

    lblName = uilabel(subGrid, 'Text', 'Simulation name');
    lblName.FontName = 'Times New Roman';
    lblName.FontSize = 16;
    lblName.Layout.Row = 3;
    lblName.Layout.Column = 4;
    entryName = uieditfield(subGrid, 'text');
    entryName.FontName = 'Times New Roman';
    entryName.FontSize = 16;
    entryName.Layout.Row = 3;
    entryName.Layout.Column = 3;
    entryName.Value = simParams.initParams.Simulation.name;

    lblType = uilabel(subGrid, 'Text', 'Simulation type');
    lblType.FontName = 'Times New Roman';
    lblType.FontSize = 16;
    lblType.Layout.Row = 4;
    lblType.Layout.Column = 4;
    entryType = uidropdown(subGrid, 'Items', {'Nonlinear', 'Linear'});
    entryType.FontName = 'Times New Roman';
    entryType.FontSize = 16;
    entryType.Layout.Row = 4;
    entryType.Layout.Column = 3;
    entryType.Value = simParams.initParams.Simulation.type;

    entryTime = add_numeric_entry(subGrid, 5, 'Initial time ($t_0$) [s]', simParams.initParams.Simulation.initialTime);
    entryFinalTime = add_numeric_entry(subGrid, 6, 'Final time ($t_f$) [s]', simParams.initParams.Simulation.finalTime);
    entryTimeStep = add_numeric_entry(subGrid, 7, 'Time step ($\Delta t$) [s]', simParams.initParams.Simulation.timeStep);
    entryRelTol = add_numeric_entry(subGrid, 8, 'Relative tolerance ($\epsilon_{rel}$)', simParams.initParams.Simulation.relTol);
    entryAbsTol = add_numeric_entry(subGrid, 9, 'Absolute tolerance ($\epsilon_{abs}$)', simParams.initParams.Simulation.absTol);

    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function save_and_close()
        currentParams = guidata(subFigHandle);
        currentParams.initParams.Simulation.name = entryName.Value;
        currentParams.initParams.Simulation.type = entryType.Value;
        currentParams.initParams.Simulation.initialTime = entryTime.Value;
        currentParams.initParams.Simulation.finalTime = entryFinalTime.Value;
        currentParams.initParams.Simulation.timeStep = entryTimeStep.Value;
        currentParams.initParams.Simulation.relTol = entryRelTol.Value;
        currentParams.initParams.Simulation.absTol = entryAbsTol.Value;
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end

function entry = add_numeric_entry(parent, rowNumber, labelText, defaultValue)
    label = uilabel(parent, 'Text', labelText, 'Interpreter', 'latex');
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = rowNumber;
    label.Layout.Column = 4;

    entry = uieditfield(parent, 'numeric');
    entry.FontName = 'Times New Roman';
    entry.FontSize = 14;
    entry.Layout.Row = rowNumber;
    entry.Layout.Column = 3;
    entry.Value = defaultValue;
end
