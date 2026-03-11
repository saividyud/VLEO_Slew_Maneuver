function set_controller_parameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    simParams = guidata(subFigHandle);
    [~, controllerName] = vleo.control.resolve_controller_handle(simParams.initParams.Controller);

    subFig = uifigure('Name', 'Set Controller Parameters', 'Position', [500, 300, 800, 450]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 30, 30, 30, 30, 30, 30, '1x', 30};
    subGrid.ColumnWidth = {'fit', 100, 260, 'fit', '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Controller Parameters');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 20;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 5];
    titleLabel.HorizontalAlignment = 'center';

    btnSaveClose = uibutton(subGrid, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10;
    btnSaveClose.Layout.Column = 1;
    btnSaveClose.BackgroundColor = salmonColor;

    label = uilabel(subGrid, 'Text', 'Controller function');
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = 3;
    label.Layout.Column = 4;

    entryFunc = uieditfield(subGrid, 'text');
    entryFunc.FontName = 'Times New Roman';
    entryFunc.FontSize = 14;
    entryFunc.Layout.Row = 3;
    entryFunc.Layout.Column = 3;
    entryFunc.Value = controllerName;

    btnBrowse = uibutton(subGrid, 'Text', 'Browse');
    btnBrowse.FontName = 'Times New Roman';
    btnBrowse.FontSize = 14;
    btnBrowse.Layout.Row = 3;
    btnBrowse.Layout.Column = 2;
    btnBrowse.ButtonPushedFcn = @(~, ~) browse_function();

    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function browse_function()
        controllerPath = vleo.util.project_path('src', '+vleo', '+control');
        [fileName, pathName] = uigetfile(fullfile(controllerPath, '*.m'), 'Select Controller Function File');
        if isequal(fileName, 0)
            return;
        end
        entryFunc.Value = vleo.control.normalize_controller_name(fullfile(pathName, fileName));
    end

    function save_and_close()
        currentParams = guidata(subFigHandle);
        controllerSpec = entryFunc.Value;
        [controllerFunc, normalizedName] = vleo.control.resolve_controller_handle(controllerSpec);
        currentParams.initParams.Controller.functionName = normalizedName;
        currentParams.initParams.Controller.functionFile = normalizedName;
        currentParams.initParams.Controller.Func = controllerFunc;
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end
