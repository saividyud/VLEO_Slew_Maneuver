function run_simulation(subFigHandle)
    simParams = guidata(subFigHandle);
    simParams = vleo.gui.ensure_simulation_defaults(simParams);
    guidata(subFigHandle, simParams);

    X0 = vleo.dynamics.initial_state_from_sim_params(simParams);
    t0 = simParams.initParams.Simulation.initialTime;
    tf = simParams.initParams.Simulation.finalTime;
    dt = simParams.initParams.Simulation.timeStep;
    ts = t0:dt:tf;
    opts = odeset('RelTol', simParams.initParams.Simulation.relTol, ...
        'AbsTol', simParams.initParams.Simulation.absTol);

    isNonlinear = strcmpi(simParams.initParams.Simulation.type, 'Nonlinear');
    if isNonlinear
        odeFun = @(t, X) vleo.dynamics.sat_dynamics_gui(t, X, subFigHandle);
    else
        Xr = [1; 0; 0; 0; 0; 0; 0];
        odeFun = @(t, X) vleo.dynamics.sat_dynamics_linearized(t, X, Xr);
    end

    progressDialog = uiprogressdlg(subFigHandle, ...
        'Title', 'Running Simulation', ...
        'Message', sprintf('Integrating %s equations of motion...', simParams.initParams.Simulation.type), ...
        'Indeterminate', 'on');

    try
        [tOut, XOut] = ode45(odeFun, ts, X0, opts);

        results = struct();
        results.t = tOut;
        results.X = XOut;
        results.rs = XOut(:, 1:3);
        results.vs = XOut(:, 4:6);
        results.betas = XOut(:, 7:10);
        results.omegas = XOut(:, 11:13);

        nSteps = length(tOut);
        torques = zeros(nSteps, 3);
        aeroTorques = zeros(nSteps, 3);
        totalTorques = zeros(nSteps, 3);
        aeroConfig = vleo.dynamics.build_aero_config_from_sim_params(simParams);

        for k = 1:nSteps
            tauControl = vleo.dynamics.evaluate_controller_torque(simParams, tOut(k), XOut(k, :)');
            if isNonlinear && vleo.util.get_mode_flag(simParams, 'enableAero', true)
                [~, tauAero] = vleo.dynamics.compute_aero_effects(tOut(k), XOut(k, :)', aeroConfig);
            else
                tauAero = [0; 0; 0];
            end
            torques(k, :) = tauControl';
            aeroTorques(k, :) = tauAero';
            totalTorques(k, :) = (tauControl + tauAero)';
        end

        results.torques = torques;
        results.aeroTorques = aeroTorques;
        results.totalTorques = totalTorques;
        results.maxAeroTorqueNorm = max(sqrt(sum(aeroTorques.^2, 2)));
        results.torqueCsvPath = save_torque_breakdown_csv(simParams, tOut, torques, aeroTorques, totalTorques);

        simParams.results = results;
        guidata(subFigHandle, simParams);

        close(progressDialog);
        if isempty(results.torqueCsvPath)
            successMessage = sprintf('Simulation completed successfully!\n%d time steps integrated.\nMax |tau_aero| = %.3e N m\nTorque breakdown CSV could not be saved.', ...
                nSteps, results.maxAeroTorqueNorm);
        else
            [~, csvName, csvExt] = fileparts(results.torqueCsvPath);
            successMessage = sprintf('Simulation completed successfully!\n%d time steps integrated.\nMax |tau_aero| = %.3e N m\nTorque breakdown saved to simulations/%s%s', ...
                nSteps, results.maxAeroTorqueNorm, csvName, csvExt);
        end

        uialert(subFigHandle, successMessage, 'Success', 'Icon', 'success');
    catch ME
        close(progressDialog);
        uialert(subFigHandle, sprintf('Simulation failed:\n%s', ME.message), 'Error', 'Icon', 'error');
    end
end

function csvPath = save_torque_breakdown_csv(simParams, tOut, controlTorques, aeroTorques, totalTorques)
    simName = simParams.initParams.Simulation.name;
    csvPath = fullfile(vleo.util.simulations_dir(), [vleo.util.sanitize_filename(simName), '_torque_breakdown.csv']);

    torqueTable = table( ...
        tOut, ...
        controlTorques(:, 1), controlTorques(:, 2), controlTorques(:, 3), ...
        aeroTorques(:, 1), aeroTorques(:, 2), aeroTorques(:, 3), ...
        totalTorques(:, 1), totalTorques(:, 2), totalTorques(:, 3), ...
        'VariableNames', { ...
            'time_s', ...
            'tau_control_x', 'tau_control_y', 'tau_control_z', ...
            'tau_aero_x', 'tau_aero_y', 'tau_aero_z', ...
            'tau_total_x', 'tau_total_y', 'tau_total_z'});

    try
        writetable(torqueTable, csvPath);
    catch
        csvPath = '';
    end
end
