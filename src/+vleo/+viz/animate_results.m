function animate_results(source)
    simParams = resolve_simulation_source(source);
    results = simParams.results;

    earthRadius = 6378.14e3;
    rs = results.rs;
    vs = results.vs;
    betas = results.betas;
    omegas = results.omegas;
    ts = results.t;
    torques = results.torques;
    if isfield(results, 'aeroTorques')
        aeroTorques = results.aeroTorques;
    else
        aeroTorques = zeros(size(torques));
    end

    r0 = rs(1, :);
    v0 = vs(1, :);
    b1Initial = v0' / norm(v0);
    b3Initial = -r0' / norm(r0);
    b2Initial = cross(b3Initial, b1Initial);
    rBodyFromInertial0 = [b1Initial'; b2Initial'; b3Initial'];

    fig = figure('Name', 'CubeSat');
    set(fig, 'CloseRequestFcn', @close_figure);
    layout = tiledlayout(fig, 6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    set(fig, 'Position', [150, 150, 1280, 720]);

    ax1 = nexttile(layout, 1, [3, 1]);
    E = wgs84Ellipsoid;
    [x, y, z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
    surf(x, y, z, 'FaceAlpha', 'texturemap', 'FaceColor', 'texturemap', 'EdgeAlpha', 'texturemap', 'Parent', ax1);
    hold(ax1, 'on');
    plot3(ax1, rs(:, 1), rs(:, 2), rs(:, 3), 'k--', 'LineWidth', 2);
    posMarker = scatter3(ax1, rs(1, 1), rs(1, 2), rs(1, 3), 'r', 'filled');
    hold(ax1, 'off');
    axis(ax1, 'equal');
    view(ax1, 30, 30);
    xlabel(ax1, 'X Position (m)');
    ylabel(ax1, 'Y Position (m)');
    zlabel(ax1, 'Z Position (m)');
    title(ax1, 'Satellite Trajectory');
    bounds = 1.1 * earthRadius;
    xlim(ax1, [-bounds, bounds]);
    ylim(ax1, [-bounds, bounds]);
    zlim(ax1, [-bounds, bounds]);
    grid(ax1, 'on');

    ax2 = nexttile(layout, 7, [3, 1]);
    scatter3(ax2, [0, 0, 0], [0, 0, 0], [0, 0, 0], 40, 'black', '+');
    hold(ax2, 'on');
    i1 = [1; 0; 0];
    i2 = [0; 1; 0];
    i3 = [0; 0; 1];
    inertialArrows = draw_frame([i1, i2, i3], 'black', ax2);
    text(i1(1), i1(2), i1(3), 'i_1', 'Color', 'k', 'Parent', ax2);
    text(i2(1), i2(2), i2(3), 'i_2', 'Color', 'k', 'Parent', ax2);
    text(i3(1), i3(2), i3(3), 'i_3', 'Color', 'k', 'Parent', ax2);

    pos = rs(1, :);
    vel = vs(1, :);
    o1 = pos' / norm(pos);
    o2Guess = vel' / norm(vel);
    o3 = cross(o1, o2Guess);
    o2 = cross(o3, o1);
    orbitalArrows = draw_frame([o1, o2, o3], 'red', ax2);
    o1Text = text(o1(1), o1(2), o1(3), 'o_1', 'Color', 'r', 'Parent', ax2);
    o2Text = text(o2(1), o2(2), o2(3), 'o_2', 'Color', 'r', 'Parent', ax2);
    o3Text = text(o3(1), o3(2), o3(3), 'o_3', 'Color', 'r', 'Parent', ax2);

    vehicle = vleo.dynamics.default_vehicle_config();
    data = stlread(vehicle.stlFilePath);
    triangleMat = data.ConnectivityList;
    points = 2 * data.Points;
    points = points - mean(points, 1);
    rStlToInertial = [0, 0, 1; 0, 1, 0; -1, 0, 0];
    points = rStlToInertial * points';
    pointsInitial = rBodyFromInertial0' * points;
    cubesat = trimesh(triangleMat, pointsInitial(1, :), pointsInitial(2, :), pointsInitial(3, :), ...
        'FaceColor', 'blue', 'FaceAlpha', 0.1, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5, 'Parent', ax2);

    bodyArrows = draw_frame([b1Initial, b2Initial, b3Initial], 'blue', ax2);
    b1Text = text(b1Initial(1), b1Initial(2), b1Initial(3), 'b_1', 'Color', 'b', 'Parent', ax2);
    b2Text = text(b2Initial(1), b2Initial(2), b2Initial(3), 'b_2', 'Color', 'b', 'Parent', ax2);
    b3Text = text(b3Initial(1), b3Initial(2), b3Initial(3), 'b_3', 'Color', 'b', 'Parent', ax2);
    traj = plot3(ax2, rs(:, 1) - rs(1, 1), rs(:, 2) - rs(1, 2), rs(:, 3) - rs(1, 3), 'k--', 'LineWidth', 2);
    axis(ax2, 'equal');
    grid(ax2, 'on');
    view(ax2, 30, 30);
    xlim(ax2, [-1, 1]);
    ylim(ax2, [-1, 1]);
    zlim(ax2, [-1, 1]);
    xlabel(ax2, 'X');
    ylabel(ax2, 'Y');
    zlabel(ax2, 'Z');
    title(ax2, 'Inertial and Orbital Reference Frames');
    legend(ax2, [inertialArrows{1}, orbitalArrows{1}, bodyArrows{1}], {'Inertial Frame', 'Orbital Frame', 'Body Frame'}, 'Location', 'northwest');

    timeAnn = annotation(fig, 'textbox', [0.1, 0.9, 0.06, 0.04], 'String', {sprintf('Time = %.2f min', ts(1) / 60)}, 'FitBoxToText', 'off');
    timeAnn.HorizontalAlignment = 'center';
    timeAnn.VerticalAlignment = 'middle';
    timeAnn.BackgroundColor = 'white';

    eulerDeg = vleo.util.quat_history_to_euler_deg(betas);
    ax3 = nexttile(layout, 2, [2, 1]);
    plot(ax3, ts / 60, eulerDeg(:, 1), 'r', 'LineWidth', 1);
    hold(ax3, 'on');
    plot(ax3, ts / 60, eulerDeg(:, 2), 'b', 'LineWidth', 1);
    plot(ax3, ts / 60, eulerDeg(:, 3), 'k', 'LineWidth', 1);
    rollMarker = scatter(ax3, ts(1) / 60, eulerDeg(1, 1), 'r', 'filled');
    pitchMarker = scatter(ax3, ts(1) / 60, eulerDeg(1, 2), 'b', 'filled');
    yawMarker = scatter(ax3, ts(1) / 60, eulerDeg(1, 3), 'k', 'filled');
    ax3Line = xline(ax3, ts(1) / 60, 'k--', 'LineWidth', 1);
    hold(ax3, 'off');
    xlabel(ax3, 'Time [min]');
    ylabel(ax3, 'Angle [deg]');
    grid(ax3, 'on');
    set_limits(ax3, eulerDeg(:));
    legend(ax3, [rollMarker, pitchMarker, yawMarker], {'Roll \phi', 'Pitch \theta', 'Yaw \psi'}, 'Location', 'best');

    ax4 = nexttile(layout, 6, [2, 1]);
    omegasDeg = rad2deg(omegas);
    plot(ax4, ts / 60, omegasDeg(:, 1), 'r', 'LineWidth', 1);
    hold(ax4, 'on');
    plot(ax4, ts / 60, omegasDeg(:, 2), 'b', 'LineWidth', 1);
    plot(ax4, ts / 60, omegasDeg(:, 3), 'k', 'LineWidth', 1);
    rollRateMarker = scatter(ax4, ts(1) / 60, omegasDeg(1, 1), 'r', 'filled');
    pitchRateMarker = scatter(ax4, ts(1) / 60, omegasDeg(1, 2), 'b', 'filled');
    yawRateMarker = scatter(ax4, ts(1) / 60, omegasDeg(1, 3), 'k', 'filled');
    ax4Line = xline(ax4, ts(1) / 60, 'k--', 'LineWidth', 1);
    hold(ax4, 'off');
    xlabel(ax4, 'Time [min]');
    ylabel(ax4, 'Rotation Rate [deg/s]');
    grid(ax4, 'on');
    set_limits(ax4, omegasDeg(:));
    legend(ax4, [rollRateMarker, pitchRateMarker, yawRateMarker], {'\omega_1', '\omega_2', '\omega_3'}, 'Location', 'best');

    ax5 = nexttile(layout, 10, [2, 1]);
    plot(ax5, ts / 60, torques(:, 1), 'r', 'LineWidth', 1);
    hold(ax5, 'on');
    plot(ax5, ts / 60, torques(:, 2), 'b', 'LineWidth', 1);
    plot(ax5, ts / 60, torques(:, 3), 'k', 'LineWidth', 1);
    plot(ax5, ts / 60, aeroTorques(:, 1), 'r--', 'LineWidth', 1);
    plot(ax5, ts / 60, aeroTorques(:, 2), 'b--', 'LineWidth', 1);
    plot(ax5, ts / 60, aeroTorques(:, 3), 'k--', 'LineWidth', 1);
    tau1Marker = scatter(ax5, ts(1) / 60, torques(1, 1), 'r', 'filled');
    tau2Marker = scatter(ax5, ts(1) / 60, torques(1, 2), 'b', 'filled');
    tau3Marker = scatter(ax5, ts(1) / 60, torques(1, 3), 'k', 'filled');
    ax5Line = xline(ax5, ts(1) / 60, 'k--', 'LineWidth', 1);
    hold(ax5, 'off');
    xlabel(ax5, 'Time [min]');
    ylabel(ax5, 'Torque [Nm]');
    grid(ax5, 'on');
    set_limits(ax5, [torques(:); aeroTorques(:)]);
    legend(ax5, {'\tau_{c,1}', '\tau_{c,2}', '\tau_{c,3}', '\tau_{a,1}', '\tau_{a,2}', '\tau_{a,3}'}, 'Location', 'best');

    set(findall(fig, '-property', 'FontSize'), 'FontSize', 24);

    for idx = 1:5:length(rs)
        if ~isvalid(fig)
            return;
        end

        timeAnn.String = {sprintf('Time = %.2f min', ts(idx) / 60)};
        pos = rs(idx, :);
        vel = vs(idx, :);
        posMarker.XData = pos(1);
        posMarker.YData = pos(2);
        posMarker.ZData = pos(3);

        o1 = pos' / norm(pos);
        o2 = vel' / norm(vel);
        o3 = cross(o1, o2);
        set(orbitalArrows{1}, 'UData', o1(1), 'VData', o1(2), 'WData', o1(3));
        set(orbitalArrows{2}, 'UData', o2(1), 'VData', o2(2), 'WData', o2(3));
        set(orbitalArrows{3}, 'UData', o3(1), 'VData', o3(2), 'WData', o3(3));
        o1Text.Position = o1;
        o2Text.Position = o2;
        o3Text.Position = o3;

        rBodyFromInertial = quat2dcm(betas(idx, :));
        b1 = rBodyFromInertial(1, :)';
        b2 = rBodyFromInertial(2, :)';
        b3 = rBodyFromInertial(3, :)';
        set(bodyArrows{1}, 'UData', b1(1), 'VData', b1(2), 'WData', b1(3));
        set(bodyArrows{2}, 'UData', b2(1), 'VData', b2(2), 'WData', b2(3));
        set(bodyArrows{3}, 'UData', b3(1), 'VData', b3(2), 'WData', b3(3));
        b1Text.Position = b1;
        b2Text.Position = b2;
        b3Text.Position = b3;
        cubesat.Vertices = (rBodyFromInertial' * points)';

        followBounds = 3000e3;
        xlim(ax1, [pos(1) - followBounds, pos(1) + followBounds]);
        ylim(ax1, [pos(2) - followBounds, pos(2) + followBounds]);
        zlim(ax1, [pos(3) - followBounds, pos(3) + followBounds]);
        traj.XData = rs(:, 1) - pos(1);
        traj.YData = rs(:, 2) - pos(2);
        traj.ZData = rs(:, 3) - pos(3);

        rollMarker.XData = ts(idx) / 60;
        rollMarker.YData = eulerDeg(idx, 1);
        pitchMarker.XData = ts(idx) / 60;
        pitchMarker.YData = eulerDeg(idx, 2);
        yawMarker.XData = ts(idx) / 60;
        yawMarker.YData = eulerDeg(idx, 3);
        ax3Line.Value = ts(idx) / 60;

        rollRateMarker.XData = ts(idx) / 60;
        rollRateMarker.YData = omegasDeg(idx, 1);
        pitchRateMarker.XData = ts(idx) / 60;
        pitchRateMarker.YData = omegasDeg(idx, 2);
        yawRateMarker.XData = ts(idx) / 60;
        yawRateMarker.YData = omegasDeg(idx, 3);
        ax4Line.Value = ts(idx) / 60;

        tau1Marker.XData = ts(idx) / 60;
        tau1Marker.YData = torques(idx, 1);
        tau2Marker.XData = ts(idx) / 60;
        tau2Marker.YData = torques(idx, 2);
        tau3Marker.XData = ts(idx) / 60;
        tau3Marker.YData = torques(idx, 3);
        ax5Line.Value = ts(idx) / 60;
        drawnow;
    end
end

function simParams = resolve_simulation_source(source)
    if isstruct(source)
        simParams = source;
    else
        simParams = guidata(source);
    end

    if ~isfield(simParams, 'results')
        error('VLEO_Slew_Maneuver:NoResults', 'No simulation results were provided for animation.');
    end
end

function set_limits(ax, values)
    span = abs(max(values) - min(values));
    if span == 0
        span = 0.01;
    end
    ylim(ax, [min(values) - 0.1 * span, max(values) + 0.1 * span]);
end

function arrows = draw_frame(vectors, color, parent)
    arrows = cell(1, 3);
    for idx = 1:3
        arrows{idx} = quiver3(0, 0, 0, vectors(1, idx), vectors(2, idx), vectors(3, idx), color, 'LineWidth', 2, 'Parent', parent);
    end
end

function close_figure(src, ~)
    delete(src);
end
