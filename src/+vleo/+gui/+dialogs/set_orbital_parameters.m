function set_orbital_parameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    c = vleo.util.constants();
    earthRadius = c.R_earth;
    simParams = guidata(subFigHandle);

    subFig = uifigure('Name', 'Set Orbital Parameters', 'Position', [500, 300, 800, 450]);
    subGrid = uigridlayout(subFig);
    subGrid.RowHeight = {30, '1x', 45, 45, 45, 45, 45, 45, '1x', 30};
    subGrid.ColumnWidth = {100, 160, '1x', 'fit', 'fit', 300, '1x'};

    titleLabel = uilabel(subGrid, 'Text', 'Set Orbital Parameters');
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
    ax.XLabel.String = 'X';
    ax.YLabel.String = 'Y';
    ax.ZLabel.String = 'Z';
    ax.Title.String = 'Current Orbit';
    ax.Title.FontName = 'Times New Roman';
    ax.Title.FontSize = 14;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    ax.DataAspectRatio = [1, 1, 1];
    E = wgs84Ellipsoid;
    [x, y, z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis, 50);
    surf(ax, x, y, z, 'FaceAlpha', 'texturemap', 'FaceColor', 'texturemap', 'EdgeAlpha', 'texturemap');
    hold(ax, 'on');

    entryAltitude = add_orbit_numeric_entry(subGrid, 3, '$h$ [m]', simParams.initParams.Orbit.altitude, [100e3, 1000e3]);
    entryEcc = add_orbit_numeric_entry(subGrid, 4, '$e$', simParams.initParams.Orbit.eccentricity, [0, 0.9]);
    entryInc = add_orbit_numeric_entry(subGrid, 5, '$i$ [deg]', simParams.initParams.Orbit.inclination, [0, 180]);
    entryAOP = add_orbit_numeric_entry(subGrid, 6, '$\omega$ [deg]', simParams.initParams.Orbit.argPeriapse, [0, 360]);
    entryRAAN = add_orbit_numeric_entry(subGrid, 7, '$\Omega$ [deg]', simParams.initParams.Orbit.RAAN, [0, 360]);
    entryTA = add_orbit_numeric_entry(subGrid, 8, '$\nu$ [deg]', simParams.initParams.Orbit.trueAnomaly, [0, 360]);

    sliderAltitude = add_slider(subGrid, 3, [100e3, 1000e3], simParams.initParams.Orbit.altitude, {'100 km', '400 km', '700 km', '1,000 km'});
    sliderAltitude.MajorTicks = 100e3:300e3:1000e3;
    sliderEcc = add_slider(subGrid, 4, [0, 0.9], simParams.initParams.Orbit.eccentricity, {'0', '0.3', '0.6', '0.9'});
    sliderEcc.MajorTicks = 0:0.3:0.9;
    sliderInc = add_slider(subGrid, 5, [0, 180], simParams.initParams.Orbit.inclination, {'0°', '60°', '120°', '180°'});
    sliderInc.MajorTicks = 0:60:180;
    sliderAOP = add_slider(subGrid, 6, [0, 360], simParams.initParams.Orbit.argPeriapse, {'0°', '90°', '180°', '270°', '360°'});
    sliderAOP.MajorTicks = 0:90:360;
    sliderRAAN = add_slider(subGrid, 7, [0, 360], simParams.initParams.Orbit.RAAN, {'0°', '90°', '180°', '270°', '360°'});
    sliderRAAN.MajorTicks = 0:90:360;
    sliderTA = add_slider(subGrid, 8, [0, 360], simParams.initParams.Orbit.trueAnomaly, {'0°', '90°', '180°', '270°', '360°'});
    sliderTA.MajorTicks = 0:90:360;

    entryStruct = struct('altitude', entryAltitude, 'eccentricity', entryEcc, 'inclination', entryInc, ...
        'argPeriapse', entryAOP, 'RAAN', entryRAAN, 'trueAnomaly', entryTA);

    orbitLine = plot3(ax, 0, 0, 0, 'r', 'LineWidth', 2);
    orbitPosition = plot3(ax, 0, 0, 0, 'o', 'MarkerFaceColor', 'k');
    refresh_plot();
    hold(ax, 'off');

    bind_slider(entryAltitude, sliderAltitude);
    bind_slider(entryEcc, sliderEcc);
    bind_slider(entryInc, sliderInc);
    bind_slider(entryAOP, sliderAOP);
    bind_slider(entryRAAN, sliderRAAN);
    bind_slider(entryTA, sliderTA);

    btnSaveClose.ButtonPushedFcn = @(~, ~) save_and_close();

    function bind_slider(entry, slider)
        slider.ValueChangingFcn = @(~, event) update_from_slider(event.Value, entry, slider);
        entry.ValueChangedFcn = @(~, event) update_from_entry(event.Value, entry, slider);
    end

    function update_from_slider(value, entry, slider)
        entry.Value = value;
        slider.Value = value;
        refresh_plot();
    end

    function update_from_entry(value, entry, slider)
        entry.Value = value;
        slider.Value = value;
        refresh_plot();
    end

    function refresh_plot()
        a = entryStruct.altitude.Value + earthRadius;
        e = entryStruct.eccentricity.Value;
        inc = deg2rad(entryStruct.inclination.Value);
        raan = deg2rad(entryStruct.RAAN.Value);
        aop = deg2rad(entryStruct.argPeriapse.Value);
        nu = deg2rad(linspace(entryStruct.trueAnomaly.Value, entryStruct.trueAnomaly.Value + 360, 500));
        r = a * (1 - e^2) ./ (1 + e * cos(nu));
        rOrbital = [r .* cos(nu); r .* sin(nu); zeros(size(nu))];
        r3Omega = [cos(raan), -sin(raan), 0; sin(raan), cos(raan), 0; 0, 0, 1];
        r1Inc = [1, 0, 0; 0, cos(inc), -sin(inc); 0, sin(inc), cos(inc)];
        r3OmegaSmall = [cos(aop), -sin(aop), 0; sin(aop), cos(aop), 0; 0, 0, 1];
        rInertial = r3Omega * r1Inc * r3OmegaSmall * rOrbital;
        set(orbitLine, 'XData', rInertial(1, :), 'YData', rInertial(2, :), 'ZData', rInertial(3, :));
        set(orbitPosition, 'XData', rInertial(1, 1), 'YData', rInertial(2, 1), 'ZData', rInertial(3, 1));
        ax.XLim = [-1.1, 1.1] * (a * (1 + e));
        ax.YLim = [-1.1, 1.1] * (a * (1 + e));
        ax.ZLim = [-1.1, 1.1] * (a * (1 + e));
        view(ax, 30, 30);
    end

    function save_and_close()
        currentParams = guidata(subFigHandle);
        currentParams.initParams.Orbit.altitude = entryStruct.altitude.Value;
        currentParams.initParams.Orbit.semiMajorAxis = entryStruct.altitude.Value + earthRadius;
        currentParams.initParams.Orbit.eccentricity = entryStruct.eccentricity.Value;
        currentParams.initParams.Orbit.inclination = entryStruct.inclination.Value;
        currentParams.initParams.Orbit.argPeriapse = entryStruct.argPeriapse.Value;
        currentParams.initParams.Orbit.RAAN = entryStruct.RAAN.Value;
        currentParams.initParams.Orbit.trueAnomaly = entryStruct.trueAnomaly.Value;
        guidata(subFigHandle, currentParams);
        vleo.gui.mark_simulation_dirty(subFigHandle);
        delete(subFig);
    end
end

function entry = add_orbit_numeric_entry(parent, rowNumber, labelText, defaultValue, limits)
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

function slider = add_slider(parent, rowNumber, limits, defaultValue, tickLabels)
    slider = uislider(parent);
    slider.Layout.Row = rowNumber;
    slider.Layout.Column = 6;
    slider.Limits = limits;
    slider.Value = defaultValue;
    slider.MajorTickLabels = tickLabels;
end
