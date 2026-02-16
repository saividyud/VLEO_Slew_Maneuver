function displayResults(subFigHandle)
% DISPLAYRESULTS Plots the simulation results stored in the figure's
% guidata. Creates a multi-panel figure showing Euler angles, angular
% rates, control torques, and the orbital trajectory.
%
% Parameters
% ----------
% subFigHandle : matlab.ui.Figure
%   Handle to the simulation configuration figure (must contain
%   simParams.results in guidata after a successful runSimulation call)

    %% Check for results
    simParams = guidata(subFigHandle);

    if ~isfield(simParams, 'results')
        uialert(subFigHandle, ...
            'No simulation results found. Please run the simulation first.', ...
            'No Results', 'Icon', 'warning');
        return;
    end

    results = simParams.results;
    t       = results.t;
    rs      = results.rs;
    betas   = results.betas;
    omegas  = results.omegas;
    torques = results.torques;

    %% Convert time to minutes for plotting
    t_min = t / 60;

    %% Compute Euler angles from quaternions
    nSteps = length(t);
    rolls   = zeros(nSteps, 1);
    pitches = zeros(nSteps, 1);
    yaws    = zeros(nSteps, 1);

    for k = 1:nSteps
        R_BI = DCMfromQ(betas(k, :)');
        % ZYX Euler angle extraction from DCM
        pitches(k) = -asind(R_BI(1, 3));
        rolls(k)   = atan2d(R_BI(2, 3), R_BI(3, 3));
        yaws(k)    = atan2d(R_BI(1, 2), R_BI(1, 1));
    end

    % Convert angular rates to deg/s for display
    omegas_deg = rad2deg(omegas);

    %% Create results figure
    simName = simParams.params.Simulation.name;
    fig = figure('Name', ['Results - ' simName], ...
                 'NumberTitle', 'off');
    fig.WindowState = 'maximized';

    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['Simulation Results: ' simName], ...
          'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

    %% Panel 1: Euler Angles
    ax1 = nexttile(tl);
    hold(ax1, 'on');
    plot(ax1, t_min, rolls,   'r', 'LineWidth', 1.2);
    plot(ax1, t_min, pitches, 'b', 'LineWidth', 1.2);
    plot(ax1, t_min, yaws,    'k', 'LineWidth', 1.2);
    hold(ax1, 'off');
    xlabel(ax1, 'Time [min]');
    ylabel(ax1, 'Angle [deg]');
    title(ax1, 'Euler Angles');
    legend(ax1, {'Roll \phi', 'Pitch \theta', 'Yaw \psi'}, 'Location', 'best');
    grid(ax1, 'on');

    %% Panel 2: Angular Velocity
    ax2 = nexttile(tl);
    hold(ax2, 'on');
    plot(ax2, t_min, omegas_deg(:, 1), 'r', 'LineWidth', 1.2);
    plot(ax2, t_min, omegas_deg(:, 2), 'b', 'LineWidth', 1.2);
    plot(ax2, t_min, omegas_deg(:, 3), 'k', 'LineWidth', 1.2);
    hold(ax2, 'off');
    xlabel(ax2, 'Time [min]');
    ylabel(ax2, 'Angular Rate [deg/s]');
    title(ax2, 'Angular Velocity');
    legend(ax2, {'\omega_1', '\omega_2', '\omega_3'}, 'Location', 'best');
    grid(ax2, 'on');

    %% Panel 3: Control Torques
    ax3 = nexttile(tl);
    hold(ax3, 'on');
    plot(ax3, t_min, torques(:, 1), 'r', 'LineWidth', 1.2);
    plot(ax3, t_min, torques(:, 2), 'b', 'LineWidth', 1.2);
    plot(ax3, t_min, torques(:, 3), 'k', 'LineWidth', 1.2);
    hold(ax3, 'off');
    xlabel(ax3, 'Time [min]');
    ylabel(ax3, 'Torque [N m]');
    title(ax3, 'Control Torques');
    legend(ax3, {'\tau_1', '\tau_2', '\tau_3'}, 'Location', 'best');
    grid(ax3, 'on');

    %% Panel 4: Orbital Trajectory
    ax4 = nexttile(tl);

    % Draw Earth
    earthRadius = 6378.14e3; % [m]
    [xe, ye, ze] = sphere(40);
    xe = xe * earthRadius;
    ye = ye * earthRadius;
    ze = ze * earthRadius;
    surf(ax4, xe, ye, ze, ...
         'FaceColor', [0.3 0.5 0.8], 'FaceAlpha', 0.3, ...
         'EdgeColor', 'none');
    hold(ax4, 'on');

    % Plot trajectory
    plot3(ax4, rs(:, 1), rs(:, 2), rs(:, 3), 'k', 'LineWidth', 1.5);

    % Mark start and end positions
    scatter3(ax4, rs(1, 1), rs(1, 2), rs(1, 3), 60, 'g', 'filled', ...
             'DisplayName', 'Start');
    scatter3(ax4, rs(end, 1), rs(end, 2), rs(end, 3), 60, 'r', 'filled', ...
             'DisplayName', 'End');
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
