function display_results(subFigHandle)
    simParams = guidata(subFigHandle);
    if ~isfield(simParams, 'results')
        uialert(subFigHandle, 'No simulation results found. Please run the simulation first.', 'No Results', 'Icon', 'warning');
        return;
    end

    vleo.viz.plot_results_summary(simParams.results, ...
        'FigureName', ['Results - ', simParams.initParams.Simulation.name], ...
        'TitleText', ['Simulation Results: ', simParams.initParams.Simulation.name]);
end
