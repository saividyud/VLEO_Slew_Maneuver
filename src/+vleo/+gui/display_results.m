function display_results(subFigHandle)
    simParams = guidata(subFigHandle);
    if ~isfield(simParams, 'results')
        uialert(subFigHandle, 'No simulation results found. Please run the simulation first.', 'No Results', 'Icon', 'warning');
        return;
    end

    results = simParams.results;
    tMinutes = results.t / 60;
    eulerDeg = vleo.util.quat_history_to_euler_deg(results.betas);
    omegasDeg = rad2deg(results.omegas);
    torques = results.torques;
    if isfield(results, 'aeroTorques')
        aeroTorques = results.aeroTorques;
    else
        aeroTorques = zeros(size(torques));
    end

    fig = figure('Name', ['Results - ', simParams.initParams.Simulation.name], 'NumberTitle', 'off');
    fig.WindowState = 'maximized';

    tile = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tile, ['Simulation Results: ', simParams.initParams.Simulation.name], ...
        'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

    ax1 = nexttile(tile);
    hold(ax1, 'on');
    plot(ax1, tMinutes, eulerDeg(:, 1), 'r', 'LineWidth', 1.2);
    plot(ax1, tMinutes, eulerDeg(:, 2), 'b', 'LineWidth', 1.2);
    plot(ax1, tMinutes, eulerDeg(:, 3), 'k', 'LineWidth', 1.2);
    hold(ax1, 'off');
    xlabel(ax1, 'Time [min]');
    ylabel(ax1, 'Angle [deg]');
    title(ax1, 'Euler Angles');
    legend(ax1, {'Roll \phi', 'Pitch \theta', 'Yaw \psi'}, 'Location', 'best');
    grid(ax1, 'on');

    ax2 = nexttile(tile);
    hold(ax2, 'on');
    plot(ax2, tMinutes, omegasDeg(:, 1), 'r', 'LineWidth', 1.2);
    plot(ax2, tMinutes, omegasDeg(:, 2), 'b', 'LineWidth', 1.2);
    plot(ax2, tMinutes, omegasDeg(:, 3), 'k', 'LineWidth', 1.2);
    hold(ax2, 'off');
    xlabel(ax2, 'Time [min]');
    ylabel(ax2, 'Angular Rate [deg/s]');
    title(ax2, 'Angular Velocity');
    legend(ax2, {'\omega_1', '\omega_2', '\omega_3'}, 'Location', 'best');
    grid(ax2, 'on');

    ax3 = nexttile(tile);
    hold(ax3, 'on');
    plot(ax3, tMinutes, torques(:, 1), 'r', 'LineWidth', 1.2);
    plot(ax3, tMinutes, torques(:, 2), 'b', 'LineWidth', 1.2);
    plot(ax3, tMinutes, torques(:, 3), 'k', 'LineWidth', 1.2);
    plot(ax3, tMinutes, aeroTorques(:, 1), 'r--', 'LineWidth', 1.0);
    plot(ax3, tMinutes, aeroTorques(:, 2), 'b--', 'LineWidth', 1.0);
    plot(ax3, tMinutes, aeroTorques(:, 3), 'k--', 'LineWidth', 1.0);
    hold(ax3, 'off');
    xlabel(ax3, 'Time [min]');
    ylabel(ax3, 'Torque [N m]');
    title(ax3, 'Control (solid) and Aero (dashed) Torques');
    legend(ax3, {'\tau_{c,1}', '\tau_{c,2}', '\tau_{c,3}', '\tau_{a,1}', '\tau_{a,2}', '\tau_{a,3}'}, 'Location', 'best');
    grid(ax3, 'on');

    ax4 = nexttile(tile);
    earthRadius = 6378.14e3;
    [xe, ye, ze] = sphere(40);
    surf(ax4, xe * earthRadius, ye * earthRadius, ze * earthRadius, ...
        'FaceColor', [0.3, 0.5, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold(ax4, 'on');
    plot3(ax4, results.rs(:, 1), results.rs(:, 2), results.rs(:, 3), 'k', 'LineWidth', 1.5);
    scatter3(ax4, results.rs(1, 1), results.rs(1, 2), results.rs(1, 3), 60, 'g', 'filled', 'DisplayName', 'Start');
    scatter3(ax4, results.rs(end, 1), results.rs(end, 2), results.rs(end, 3), 60, 'r', 'filled', 'DisplayName', 'End');
    hold(ax4, 'off');
    axis(ax4, 'equal');
    xlabel(ax4, 'X [m]');
    ylabel(ax4, 'Y [m]');
    zlabel(ax4, 'Z [m]');
    title(ax4, 'Orbital Trajectory');
    legend(ax4, {'Earth', 'Trajectory', 'Start', 'End'}, 'Location', 'best');
    grid(ax4, 'on');
    view(ax4, 30, 30);
end
