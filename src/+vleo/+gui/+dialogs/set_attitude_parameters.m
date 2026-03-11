function set_attitude_parameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    simParams = guidata(subFigHandle);

    subFig = uifigure('Name', 'Set Attitude Parameters', 'Position', [500, 300, 800, 450]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 45, 45, 45, 45, 45, 45, '1x', 30};
    subGrid.ColumnWidth = {100, 160, '1x', 'fit', 'fit', 300, '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Attitude Parameters');
    titleLabel.FontName = 'Times New Roman';
    titleLabel.FontSize = 20;
    titleLabel.FontWeight = 'bold';
    titleLabel.Layout.Row = 1;
    titleLabel.Layout.Column = [1, 7];
    titleLabel.HorizontalAlignment = 'center';

    btnSaveClose = uibutton(subGrid, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10;
    btnSaveClose.Layout.Column = 1;
    btnSaveClose.BackgroundColor = salmonColor;

    ax = uiaxes(subGrid);
    ax.Layout.Row = [3, 8];
    ax.Layout.Column = [1, 2];
    ax.Title.String = 'Current Orientation';
    ax.Title.FontName = 'Times New Roman';
    ax.Title.FontSize = 14;
    ax.DataAspectRatio = [1, 1, 1];
    ax.XLim = [-1, 1];
    ax.YLim = [-1, 1];
    ax.ZLim = [-1, 1];
    ax.XColor = 'none';
    ax.YColor = 'none';
    ax.ZColor = 'none';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    hold(ax, 'on');

    inertialVectors = eye(3);
    [~, ~] = draw_frame(inertialVectors, 'k', ax, 'i');
    [rRefToInertial, ~, ~] = orbit_reference(simParams);
    orbitalVectors = [rRefToInertial(:, 1), rRefToInertial(:, 2), rRefToInertial(:, 3)];
    draw_frame(orbitalVectors, 'b', ax, 'o');

    entryRoll = add_attitude_entry(subGrid, 3, '$\phi$ [deg]', simParams.initParams.Attitude.roll, [-180, 180]);
    entryPitch = add_attitude_entry(subGrid, 4, '$\theta$ [deg]', simParams.initParams.Attitude.pitch, [-180, 180]);
    entryYaw = add_attitude_entry(subGrid, 5, '$\psi$ [deg]', simParams.initParams.Attitude.yaw, [-180, 180]);
    entryRollRate = add_attitude_entry(subGrid, 6, '$\dot{\phi}$ [deg/s]', simParams.initParams.Attitude.rollRate, [-180, 180]);
    entryPitchRate = add_attitude_entry(subGrid, 7, '$\dot{\theta}$ [deg/s]', simParams.initParams.Attitude.pitchRate, [-180, 180]);
    entryYawRate = add_attitude_entry(subGrid, 8, '$\dot{\psi}$ [deg/s]', simParams.initParams.Attitude.yawRate, [-180, 180]);

    sliderRoll = add_slider(subGrid, 3, simParams.initParams.Attitude.roll);
    sliderPitch = add_slider(subGrid, 4, simParams.initParams.Attitude.pitch);
    sliderYaw = add_slider(subGrid, 5, simParams.initParams.Attitude.yaw);
    sliderRollRate = add_slider(subGrid, 6, simParams.initParams.Attitude.rollRate);
    sliderPitchRate = add_slider(subGrid, 7, simParams.initParams.Attitude.pitchRate);
    sliderYawRate = add_slider(subGrid, 8, simParams.initParams.Attitude.yawRate);

    entryStruct = struct('roll', entryRoll, 'pitch', entryPitch, 'yaw', entryYaw, ...
        'rollRate', entryRollRate, 'pitchRate', entryPitchRate, 'yawRate', entryYawRate);

    bodyVectors = compute_body_vectors();
    [bodyArrows, bodyTexts] = draw_frame(bodyVectors, 'r', ax, 'b');
    hold(ax, 'off');
    view(ax, 30, 30);

    bind_slider(entryRoll, sliderRoll, true);
    bind_slider(entryPitch, sliderPitch, true);
    bind_slider(entryYaw, sliderYaw, true);
    bind_slider(entryRollRate, sliderRollRate, false);
    bind_slider(entryPitchRate, sliderPitchRate, false);
    bind_slider(entryYawRate, sliderYawRate, false);

    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function bind_slider(entry, slider, updatePlot)
        slider.ValueChangingFcn = @(~, event) update_value(event.Value, entry, slider, updatePlot);
        entry.ValueChangedFcn = @(~, event) update_value(event.Value, entry, slider, updatePlot);
    end

    function update_value(value, entry, slider, updatePlot)
        entry.Value = value;
        slider.Value = value;
        if updatePlot
            newVectors = compute_body_vectors();
            for idx = 1:3
                set(bodyArrows{idx}, 'UData', newVectors(1, idx), 'VData', newVectors(2, idx), 'WData', newVectors(3, idx));
                set(bodyTexts{idx}, 'Position', newVectors(:, idx)');
            end
        end
    end

    function vectors = compute_body_vectors()
        roll = deg2rad(entryStruct.roll.Value);
        pitch = deg2rad(entryStruct.pitch.Value);
        yaw = deg2rad(entryStruct.yaw.Value);
        rZ = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
        rY = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
        rX = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
        rBodyFromReference = rZ * rY * rX;
        vectors = rRefToInertial * rBodyFromReference;
    end

    function save_and_close()
        currentParams = guidata(subFigHandle);
        currentParams.initParams.Attitude.roll = entryStruct.roll.Value;
        currentParams.initParams.Attitude.pitch = entryStruct.pitch.Value;
        currentParams.initParams.Attitude.yaw = entryStruct.yaw.Value;
        currentParams.initParams.Attitude.rollRate = entryStruct.rollRate.Value;
        currentParams.initParams.Attitude.pitchRate = entryStruct.pitchRate.Value;
        currentParams.initParams.Attitude.yawRate = entryStruct.yawRate.Value;
        rZ = [cosd(currentParams.initParams.Attitude.yaw), -sind(currentParams.initParams.Attitude.yaw), 0; ...
            sind(currentParams.initParams.Attitude.yaw), cosd(currentParams.initParams.Attitude.yaw), 0; ...
            0, 0, 1];
        rY = [cosd(currentParams.initParams.Attitude.pitch), 0, sind(currentParams.initParams.Attitude.pitch); ...
            0, 1, 0; ...
            -sind(currentParams.initParams.Attitude.pitch), 0, cosd(currentParams.initParams.Attitude.pitch)];
        rX = [1, 0, 0; 0, cosd(currentParams.initParams.Attitude.roll), -sind(currentParams.initParams.Attitude.roll); ...
            0, sind(currentParams.initParams.Attitude.roll), cosd(currentParams.initParams.Attitude.roll)];
        currentParams.initParams.Attitude.Quaternion = dcm2quat(rRefToInertial * rZ * rY * rX);
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end

function entry = add_attitude_entry(parent, rowNumber, labelText, defaultValue, limits)
    label = uilabel(parent, 'Text', labelText, 'Interpreter', 'latex');
    label.FontName = 'Times New Roman';
    label.FontSize = 14;
    label.Layout.Row = rowNumber;
    label.Layout.Column = 5;
    entry = uieditfield(parent, 'numeric');
    entry.FontName = 'Times New Roman';
    entry.FontSize = 14;
    entry.Layout.Row = rowNumber;
    entry.Layout.Column = 4;
    entry.Limits = limits;
    entry.Value = defaultValue;
end

function slider = add_slider(parent, rowNumber, defaultValue)
    slider = uislider(parent);
    slider.Layout.Row = rowNumber;
    slider.Layout.Column = 6;
    slider.Limits = [-180, 180];
    slider.Value = defaultValue;
    slider.MajorTicks = -180:90:180;
end

function [rRefToInertial, orbitalVectors, inertialVectors] = orbit_reference(simParams)
    inertialVectors = eye(3);
    inc = deg2rad(simParams.initParams.Orbit.inclination);
    raan = deg2rad(simParams.initParams.Orbit.RAAN);
    aop = deg2rad(simParams.initParams.Orbit.argPeriapse);
    nu = deg2rad(simParams.initParams.Orbit.trueAnomaly);
    r3Omega = [cos(raan), -sin(raan), 0; sin(raan), cos(raan), 0; 0, 0, 1];
    r1Inc = [1, 0, 0; 0, cos(inc), -sin(inc); 0, sin(inc), cos(inc)];
    r3OmegaSmall = [cos(aop), -sin(aop), 0; sin(aop), cos(aop), 0; 0, 0, 1];
    rOrbitToInertial = r3Omega * r1Inc * r3OmegaSmall;
    o3 = rOrbitToInertial(:, 3);
    o1 = rOrbitToInertial(:, 1) * cos(nu) + rOrbitToInertial(:, 2) * sin(nu);
    o2 = cross(o3, o1);
    orbitalVectors = [o1, o2, o3];
    rRefToInertial = [o1 / norm(o1), o2 / norm(o2), o3 / norm(o3)];
end

function [arrows, texts] = draw_frame(vectors, color, parent, symbol)
    arrows = cell(1, 3);
    texts = cell(1, 3);
    for idx = 1:3
        arrows{idx} = quiver3(0, 0, 0, vectors(1, idx), vectors(2, idx), vectors(3, idx), color, 'LineWidth', 2, 'Parent', parent);
        arrows{idx}.AutoScale = 'off';
        texts{idx} = text(vectors(1, idx), vectors(2, idx), vectors(3, idx), [' ', symbol, '_', num2str(idx)], ...
            'FontSize', 16, 'FontName', 'Times New Roman', 'Parent', parent, 'Color', color);
    end
end
