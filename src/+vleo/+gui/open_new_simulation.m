function open_new_simulation(mainFigHandle)
    simParams = vleo.gui.default_simulation_params();
    vleo.gui.open_simulation_editor(mainFigHandle, simParams, 'new');
end
