function show_results(subFigHandle)
    simParams = guidata(subFigHandle);
    if ~isfield(simParams, 'results')
        uialert(subFigHandle, 'No simulation results found. Please run the simulation first.', 'No Results', 'Icon', 'warning');
        return;
    end

    vleo.gui.display_results(subFigHandle);
    vleo.viz.animate_results(subFigHandle);
end
