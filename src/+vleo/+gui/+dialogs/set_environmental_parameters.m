function set_environmental_parameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    simParams = guidata(subFigHandle);

    subFig = uifigure('Name', 'Set Environmental Parameters', 'Position', [500, 300, 800, 450]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 20, 20, 20, 20, 20, 20, '1x', 20, 20, 20, 20, '1x', 30};
    subGrid.ColumnWidth = {'fit', '1x', 100, 'fit', '1x', 100, 'fit', '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Environmental Parameters');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 20;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 8];
    titleLabel.HorizontalAlignment = 'center';

    btnSaveClose = uibutton(subGrid, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 15;
    btnSaveClose.Layout.Column = [1, 3];
    btnSaveClose.BackgroundColor = salmonColor;

    entrySOD = add_env_numeric(subGrid, 4, 3, 4, 'Seconds of day [s]', simParams.initParams.Environment.secondsOfDay);
    entryDOY = add_env_numeric(subGrid, 5, 3, 4, 'Day of year', simParams.initParams.Environment.dayOfYear);
    entryF107 = add_env_numeric(subGrid, 6, 3, 4, '81-day average F10.7 solar flux', simParams.initParams.Environment.F107Average);
    entryF107D = add_env_numeric(subGrid, 7, 3, 4, 'Daily F10.7 solar flux', simParams.initParams.Environment.F107Daily);

    magneticGrid = uigridlayout(subGrid, [1, 8]);
    magneticGrid.RowHeight = {20};
    magneticGrid.ColumnWidth = {30, 30, 30, 30, 30, 30, 30, 'fit'};
    magneticGrid.Layout.Row = 8;
    magneticGrid.Layout.Column = [3, 4];
    magneticGrid.Padding = [0, 0, 0, 0];
    magneticLabel = uilabel(magneticGrid, 'Text', 'Magnetic indices');
    magneticLabel.FontName = 'Times New Roman';
    magneticLabel.FontSize = 14;
    magneticLabel.Layout.Row = 1;
    magneticLabel.Layout.Column = 8;
    entryMags = gobjects(1, 7);
    for idx = 1:7
        entryMags(idx) = uieditfield(magneticGrid, 'numeric');
        entryMags(idx).FontName = 'Times New Roman';
        entryMags(idx).FontSize = 14;
        entryMags(idx).Layout.Row = 1;
        entryMags(idx).Layout.Column = idx;
        entryMags(idx).ValueDisplayFormat = '%.1f';
        entryMags(idx).HorizontalAlignment = 'center';
        entryMags(idx).Value = simParams.initParams.Environment.magneticIndices(idx);
    end

    entryGSI = add_env_dropdown(subGrid, 4, 6, 7, 'Gas-surface interaction model', {'cook', 'sentman'}, ...
        vleo.util.normalize_gsi_model(simParams.initParams.Environment.gasSurfaceInteractionModel));
    entryAccom = add_env_numeric(subGrid, 5, 6, 7, 'Accommodation coefficient', simParams.initParams.Environment.accommodationCoefficient);
    entryTemp = add_env_numeric(subGrid, 6, 6, 7, 'Wall temperature [K]', simParams.initParams.Environment.wallTemperature);
    entrySpec = add_env_numeric(subGrid, 7, 6, 7, 'Specular reflectivity', simParams.initParams.Environment.specularReflectivity);
    entryDiff = add_env_numeric(subGrid, 8, 6, 7, 'Diffuse reflectivity', simParams.initParams.Environment.diffuseReflectivity);

    entryAO = add_boolean_dropdown(subGrid, 11, 'Enable anomalous oxygen', simParams.initParams.Environment.enableAnomalousOxygen);
    entryShadow = add_boolean_dropdown(subGrid, 12, 'Enable shadow analysis', simParams.initParams.Environment.enableShadowAnalysis);
    entrySRP = add_boolean_dropdown(subGrid, 13, 'Enable solar radiation pressure', simParams.initParams.Environment.enableSolarRadiationPressure);

    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function save_and_close()
        currentParams = guidata(subFigHandle);
        currentParams.initParams.Environment.secondsOfDay = entrySOD.Value;
        currentParams.initParams.Environment.dayOfYear = entryDOY.Value;
        currentParams.initParams.Environment.F107Average = entryF107.Value;
        currentParams.initParams.Environment.F107Daily = entryF107D.Value;
        currentParams.initParams.Environment.magneticIndices = arrayfun(@(entry) entry.Value, entryMags);
        currentParams.initParams.Environment.gasSurfaceInteractionModel = entryGSI.Value;
        currentParams.initParams.Environment.accommodationCoefficient = entryAccom.Value;
        currentParams.initParams.Environment.wallTemperature = entryTemp.Value;
        currentParams.initParams.Environment.specularReflectivity = entrySpec.Value;
        currentParams.initParams.Environment.diffuseReflectivity = entryDiff.Value;
        currentParams.initParams.Environment.enableAnomalousOxygen = strcmp(entryAO.Value, 'true');
        currentParams.initParams.Environment.enableShadowAnalysis = strcmp(entryShadow.Value, 'true');
        currentParams.initParams.Environment.enableSolarRadiationPressure = strcmp(entrySRP.Value, 'true');
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end

function entry = add_env_numeric(parent, rowNumber, valueColumn, labelColumn, labelText, defaultValue)
    label = uilabel(parent, 'Text', labelText);
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = rowNumber;
    label.Layout.Column = labelColumn;
    entry = uieditfield(parent, 'numeric');
    entry.FontName = 'Times New Roman';
    entry.FontSize = 14;
    entry.Layout.Row = rowNumber;
    entry.Layout.Column = valueColumn;
    entry.Value = defaultValue;
end

function entry = add_env_dropdown(parent, rowNumber, valueColumn, labelColumn, labelText, items, defaultValue)
    label = uilabel(parent, 'Text', labelText);
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = rowNumber;
    label.Layout.Column = labelColumn;
    entry = uidropdown(parent, 'Items', items);
    entry.FontName = 'Times New Roman';
    entry.FontSize = 14;
    entry.Layout.Row = rowNumber;
    entry.Layout.Column = valueColumn;
    entry.Value = defaultValue;
end

function entry = add_boolean_dropdown(parent, rowNumber, labelText, defaultValue)
    label = uilabel(parent, 'Text', labelText);
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = rowNumber;
    label.Layout.Column = 4;
    entry = uidropdown(parent, 'Items', {'true', 'false'});
    entry.FontName = 'Times New Roman';
    entry.FontSize = 14;
    entry.Layout.Row = rowNumber;
    entry.Layout.Column = 3;
    if defaultValue
        entry.Value = 'true';
    else
        entry.Value = 'false';
    end
end
