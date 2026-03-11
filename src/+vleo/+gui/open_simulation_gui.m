function open_simulation_gui()
    evalin('base', 'clear');
    clc;
    close all;

    mainFig = uifigure('Name', 'Simulation Manager', 'Position', [500, 300, 800, 450]);
    salmonColor = [1, 0.4941, 0.4392];

    gridLayout = uigridlayout(mainFig);
    gridLayout.ColumnWidth = {100, '1x', 200, '1x', 100};
    gridLayout.RowHeight = {'fit', '1x', 'fit', 'fit', 'fit', '1x', 'fit'};

    titleLabel = uilabel(gridLayout, 'Text', 'Simulation Dashboard');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 24;
    titleLabel.FontWeight = 'bold';
    titleLabel.HorizontalAlignment = 'center';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 5];

    btnNew = uibutton(gridLayout, 'Text', 'New Simulation');
    vleo.gui.apply_button_style(btnNew);
    btnNew.Layout.Row = 3;
    btnNew.Layout.Column = 3;
    btnNew.ButtonPushedFcn = @(~, ~) vleo.gui.open_new_simulation(mainFig);

    btnLoad = uibutton(gridLayout, 'Text', 'Load Saved Simulation');
    vleo.gui.apply_button_style(btnLoad);
    btnLoad.Layout.Row = 4;
    btnLoad.Layout.Column = 3;
    btnLoad.ButtonPushedFcn = @(~, ~) vleo.gui.open_load_simulation(mainFig);

    btnResults = uibutton(gridLayout, 'Text', 'Show Results');
    vleo.gui.apply_button_style(btnResults);
    btnResults.Layout.Row = 5;
    btnResults.Layout.Column = 3;
    btnResults.ButtonPushedFcn = @(~, ~) vleo.gui.open_results_window(mainFig);

    btnQuit = uibutton(gridLayout, 'Text', 'Quit');
    vleo.gui.apply_button_style(btnQuit);
    btnQuit.Layout.Row = 7;
    btnQuit.Layout.Column = 1;
    btnQuit.BackgroundColor = salmonColor;
    btnQuit.ButtonPushedFcn = @(~, ~) delete(mainFig);
end
