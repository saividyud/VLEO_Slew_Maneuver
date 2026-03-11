function set_aero_control_modes(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    simParams = guidata(subFigHandle);
    simParams = vleo.gui.ensure_simulation_defaults(simParams);

    subFig = uifigure('Name', 'Set Aero/Control Modes', 'Position', [500, 300, 600, 320]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 40, 40, '1x', 30};
    subGrid.ColumnWidth = {'fit', 180, 220, '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Aero/Control Modes');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 20;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 4];
    titleLabel.HorizontalAlignment = 'center';

    lblAero = uilabel(subGrid, 'Text', 'Enable aerodynamic model');
    lblAero.FontName = 'Times New Roman';
    lblAero.FontSize = 14;
    lblAero.Layout.Row = 3;
    lblAero.Layout.Column = 3;

    entryAero = uidropdown(subGrid, 'Items', {'true', 'false'});
    entryAero.FontName = 'Times New Roman';
    entryAero.FontSize = 14;
    entryAero.Layout.Row = 3;
    entryAero.Layout.Column = 2;
    entryAero.Value = ternary_value(simParams.initParams.Modes.enableAero);

    lblControl = uilabel(subGrid, 'Text', 'Enable control torques');
    lblControl.FontName = 'Times New Roman';
    lblControl.FontSize = 14;
    lblControl.Layout.Row = 4;
    lblControl.Layout.Column = 3;

    entryControl = uidropdown(subGrid, 'Items', {'true', 'false'});
    entryControl.FontName = 'Times New Roman';
    entryControl.FontSize = 14;
    entryControl.Layout.Row = 4;
    entryControl.Layout.Column = 2;
    entryControl.Value = ternary_value(simParams.initParams.Modes.enableControl);

    btnSaveClose = uibutton(subGrid, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 6;
    btnSaveClose.Layout.Column = 1;
    btnSaveClose.BackgroundColor = salmonColor;
    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function save_and_close()
        currentParams = guidata(subFigHandle);
        currentParams.initParams.Modes.enableAero = strcmp(entryAero.Value, 'true');
        currentParams.initParams.Modes.enableControl = strcmp(entryControl.Value, 'true');
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end

function valueText = ternary_value(flag)
    if flag
        valueText = 'true';
    else
        valueText = 'false';
    end
end
