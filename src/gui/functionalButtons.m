function functionalButtons(subFigHandle, buttonType)
    switch buttonType
        case 'Set Simulation Parameters'
            setSimulationParameters(subFigHandle);
        case 'Set Orbital Parameters'
            setOrbitalParameters(subFigHandle);
        case 'Set Attitude Parameters'
            setAttitudeParameters(subFigHandle);
        case 'Set Environmental Parameters'
            setEnvironmentalParameters(subFigHandle);
        case 'Set Controller Parameters'
            setControllerParameters(subFigHandle);
        otherwise
            error('Unknown button type: %s', buttonType);
    end
end

function setSimulationParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Create the New Window
    subFig = uifigure('Name', 'Set Simulation Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    % subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        'fit',...
        50,...
        175,...
        'fit',...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Simulation Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 5];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 11; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining entries
    % Simulation name
    lblName = uilabel(subGl, 'Text', 'Simulation name');
    lblName.FontName = 'Times New Roman';
    lblName.FontSize = 16;
    lblName.Layout.Row = 3;
    lblName.Layout.Column = 4;
    entryName = uieditfield(subGl, 'text');
    entryName.FontName = 'Times New Roman';
    entryName.FontSize = 16;
    entryName.Layout.Row = 3;
    entryName.Layout.Column = 3;
    entryName.Value = simParams.initParams.Simulation.name; % Default name

    % Simulation type
    lblType = uilabel(subGl, 'Text', 'Simulation type');
    lblType.FontName = 'Times New Roman';
    lblType.FontSize = 16;
    lblType.Layout.Row = 4;
    lblType.Layout.Column = 4;
    entryType = uidropdown(subGl);
    entryType.Items = {'Nonlinear', 'Linear'};
    entryType.FontName = 'Times New Roman';
    entryType.FontSize = 16;
    entryType.Layout.Row = 4;
    entryType.Layout.Column = 3;
    entryType.Value = simParams.initParams.Simulation.type; % Default to Nonlinear

    % Initial time
    lblTime = uilabel(subGl, 'Text', 'Initial time ($t_0$) [s]', 'Interpreter', 'latex');
    lblTime.FontName = 'Times New Roman';
    lblTime.FontSize = 14;
    lblTime.Layout.Row = 5;
    lblTime.Layout.Column = 4;
    entryTime = uieditfield(subGl, 'numeric');
    entryTime.FontName = 'Times New Roman';
    entryTime.FontSize = 14;
    entryTime.Layout.Row = 5;
    entryTime.Layout.Column = 3;
    entryTime.Value = simParams.initParams.Simulation.initialTime; % Default to 0

    % Final time
    lblFinalTime = uilabel(subGl, 'Text', 'Final time ($t_f$) [s]', 'Interpreter', 'latex');
    lblFinalTime.FontName = 'Times New Roman';
    lblFinalTime.FontSize = 14;
    lblFinalTime.Layout.Row = 6;
    lblFinalTime.Layout.Column = 4;
    entryFinalTime = uieditfield(subGl, 'numeric');
    entryFinalTime.FontName = 'Times New Roman';
    entryFinalTime.FontSize = 14;
    entryFinalTime.Layout.Row = 6;
    entryFinalTime.Layout.Column = 3;
    entryFinalTime.Value = simParams.initParams.Simulation.finalTime; % Default to 1000 seconds

    % Time step
    lblTimeStep = uilabel(subGl, 'Text', 'Time step ($\Delta t$) [s]', 'Interpreter', 'latex');
    lblTimeStep.FontName = 'Times New Roman';
    lblTimeStep.FontSize = 14;
    lblTimeStep.Layout.Row = 7;
    lblTimeStep.Layout.Column = 4;
    entryTimeStep = uieditfield(subGl, 'numeric');
    entryTimeStep.FontName = 'Times New Roman';
    entryTimeStep.FontSize = 14;
    entryTimeStep.Layout.Row = 7;
    entryTimeStep.Layout.Column = 3;
    entryTimeStep.Value = simParams.initParams.Simulation.timeStep; % Default to 1 second

    % Relative tolerance
    lblRelTol = uilabel(subGl, 'Text', 'Relative tolerance ($\epsilon_{rel}$)', 'Interpreter', 'latex');
    lblRelTol.FontName = 'Times New Roman';
    lblRelTol.FontSize = 14;
    lblRelTol.Layout.Row = 8;
    lblRelTol.Layout.Column = 4;
    entryRelTol = uieditfield(subGl, 'numeric');
    entryRelTol.FontName = 'Times New Roman';
    entryRelTol.FontSize = 14;
    entryRelTol.Layout.Row = 8;
    entryRelTol.Layout.Column = 3;
    entryRelTol.Value = simParams.initParams.Simulation.relTol; % Default to 1e-6

    % Absolute tolerance
    lblAbsTol = uilabel(subGl, 'Text', 'Absolute tolerance ($\epsilon_{abs}$)', 'Interpreter', 'latex');
    lblAbsTol.FontName = 'Times New Roman';
    lblAbsTol.FontSize = 14;
    lblAbsTol.Layout.Row = 9;
    lblAbsTol.Layout.Column = 4;
    entryAbsTol = uieditfield(subGl, 'numeric');
    entryAbsTol.FontName = 'Times New Roman';
    entryAbsTol.FontSize = 14;
    entryAbsTol.Layout.Row = 9;
    entryAbsTol.Layout.Column = 3;
    entryAbsTol.Value = simParams.initParams.Simulation.absTol; % Default to 1e-9

    %% Saving logic
    entryStruct = struct(...
        'name', entryName, ...
        'type', entryType, ...
        'initialTime', entryTime, ...
        'finalTime', entryFinalTime, ...
        'timeStep', entryTimeStep, ...
        'relTol', entryRelTol, ...
        'absTol', entryAbsTol ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose();

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.initParams.Simulation.name = entryStruct.name.Value;
        simParams.initParams.Simulation.type = entryStruct.type.Value;
        simParams.initParams.Simulation.initialTime = entryStruct.initialTime.Value;
        simParams.initParams.Simulation.finalTime = entryStruct.finalTime.Value;
        simParams.initParams.Simulation.timeStep = entryStruct.timeStep.Value;
        simParams.initParams.Simulation.relTol = entryStruct.relTol.Value;
        simParams.initParams.Simulation.absTol = entryStruct.absTol.Value;

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update main window title with asterisk to indicate unsaved changes
        simName = simParams.initParams.Simulation.name;
        subFigHandle.Name = [simName, ' *']; % Update the window title to include the simulation name and asterisk for unsaved changes

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end

end

function setOrbitalParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];
    earthRadius = 6378.14e3; % Earth equatorial radius [m]

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Create the New Window
    subFig = uifigure('Name', 'Set Orbital Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    % subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        45,...
        45,...
        45,...
        45,...
        45,...
        45,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        100,...
        160,...
        '1x',...
        'fit',...
        'fit',...
        300,...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Orbital Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 7];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining figure to show current orbit
    % Axes for orbit plot
    ax = uiaxes(subGl);
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
    ax.DataAspectRatio = [1, 1, 1]; % Equal aspect ratio for all axes

    % Plotting Earth for reference
    % Defining Earth using WGS84 ellipsoid
    E = wgs84Ellipsoid;
    [x,y,z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis, 50); % Scale to m and use 50 points for smoothness
    surf(x, y, z, FaceAlpha="texturemap", FaceColor="texturemap", EdgeAlpha="texturemap", Parent=ax); % Plot the Earth
    hold(ax, "on");

    % Plotting current orbit
    a = simParams.initParams.Orbit.semiMajorAxis; % m
    e = simParams.initParams.Orbit.eccentricity;
    i = deg2rad(simParams.initParams.Orbit.inclination); % Convert to radians
    RAAN = deg2rad(simParams.initParams.Orbit.RAAN); % Convert to radians
    AOP = deg2rad(simParams.initParams.Orbit.argPeriapse); % Convert to radians
    nu = deg2rad(linspace(simParams.initParams.Orbit.trueAnomaly, simParams.initParams.Orbit.trueAnomaly + 360, 500)); % True anomaly from 0 to 360 degrees

    % Construct orbit in orbital frame
    r = a * (1 - e^2) ./ (1 + e * cos(nu)); % Orbital radius
    r_orbital = [r .* cos(nu); r .* sin(nu); zeros(size(nu))]; % Orbit in orbital plane

    % Creating transformation matrix from orbital frame to inertial frame
    R3_Omega = [cos(RAAN), -sin(RAAN), 0; 
                sin(RAAN), cos(RAAN), 0; 
                0, 0, 1]; % Rotation about Z-axis by RAAN

    R1_i = [1, 0, 0;
            0, cos(i), -sin(i);
            0, sin(i), cos(i)]; % Rotation about X-axis by inclination

    R3_omega = [cos(AOP), -sin(AOP), 0; 
                sin(AOP), cos(AOP), 0; 
                0, 0, 1]; % Rotation about Z-axis by argument of periapse

    R = R3_Omega * R1_i * R3_omega; % Combined rotation matrix

    % Rotate orbit from orbital frame to inertial frame
    r_inertial = R * r_orbital;
    x_orbit = r_inertial(1, :);
    y_orbit = r_inertial(2, :);
    z_orbit = r_inertial(3, :);
    
    % Plotting
    orbitLine = plot3(ax, x_orbit, y_orbit, z_orbit, 'r', 'LineWidth', 2); % Plot the orbit in red
    orbitPosition = plot3(ax, x_orbit(1), y_orbit(1), z_orbit(1), 'o', 'MarkerFaceColor', 'k'); % Mark the current position on the orbit

    ax.XLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
    ax.YLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
    ax.ZLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit

    view(ax, 30, 30)

    hold(ax, "off");

    %% Defining entries
    % Altitude
    lblAltitude = uilabel(subGl, 'Text', '$h$ [m]', 'Interpreter', 'latex');
    lblAltitude.FontName = 'Times New Roman';
    lblAltitude.FontSize = 14;
    lblAltitude.Layout.Row = 3;
    lblAltitude.Layout.Column = 5;
    entryAltitude = uieditfield(subGl, 'numeric');
    entryAltitude.FontName = 'Times New Roman';
    entryAltitude.FontSize = 14;
    entryAltitude.Layout.Row = 3;
    entryAltitude.Layout.Column = 4;
    entryAltitude.Tooltip = 'Altitude of the orbit (m)';
    entryAltitude.Limits = [100e3, 1000e3]; % Limit altitude entry to 100,000 m to 1,000,000 m
    entryAltitude.Value = simParams.initParams.Orbit.altitude; % Default to 400,000 m altitude

    % Eccentricity
    lblEcc = uilabel(subGl, 'Text', '$e$', 'Interpreter', 'latex');
    lblEcc.FontName = 'Times New Roman';
    lblEcc.FontSize = 14;
    lblEcc.Layout.Row = 4;
    lblEcc.Layout.Column = 5;
    entryEcc = uieditfield(subGl, 'numeric');
    entryEcc.FontName = 'Times New Roman';
    entryEcc.FontSize = 14;
    entryEcc.Layout.Row = 4;
    entryEcc.Layout.Column = 4;
    entryEcc.Tooltip = 'Eccentricity';
    entryEcc.Limits = [0, 0.9]; % Limit eccentricity entry to 0 to 0.9
    entryEcc.Value = simParams.initParams.Orbit.eccentricity; % Default to 0

    % Inclination
    lblInc = uilabel(subGl, 'Text', '$i$ [deg]', 'Interpreter', 'latex');
    lblInc.FontName = 'Times New Roman';
    lblInc.FontSize = 14;
    lblInc.Layout.Row = 5;
    lblInc.Layout.Column = 5;
    entryInc = uieditfield(subGl, 'numeric');
    entryInc.FontName = 'Times New Roman';
    entryInc.FontSize = 14;
    entryInc.Layout.Row = 5;
    entryInc.Layout.Column = 4;
    entryInc.Tooltip = 'Inclination';
    entryInc.Limits = [0, 180]; % Limit inclination entry to 0 to 180 degrees
    entryInc.Value = simParams.initParams.Orbit.inclination; % Default to 0

    % Argument of periapse
    lblAOP = uilabel(subGl, 'Text', '$\omega$ [deg]', 'Interpreter', 'latex');
    lblAOP.FontName = 'Times New Roman';
    lblAOP.FontSize = 14;
    lblAOP.Layout.Row = 6;
    lblAOP.Layout.Column = 5;
    entryAOP = uieditfield(subGl, 'numeric');
    entryAOP.FontName = 'Times New Roman';
    entryAOP.FontSize = 14;
    entryAOP.Layout.Row = 6;
    entryAOP.Layout.Column = 4;
    entryAOP.Tooltip = 'Argument of periapse';
    entryAOP.Limits = [0, 360]; % Limit argument of periapse entry to 0 to 360 degrees
    entryAOP.Value = simParams.initParams.Orbit.argPeriapse; % Default to 0

    % Right ascension of ascending node
    lblRAAN = uilabel(subGl, 'Text', '$\Omega$ [deg]', 'Interpreter', 'latex');
    lblRAAN.FontName = 'Times New Roman';
    lblRAAN.FontSize = 14;
    lblRAAN.Layout.Row = 7;
    lblRAAN.Layout.Column = 5;
    entryRAAN = uieditfield(subGl, 'numeric');
    entryRAAN.FontName = 'Times New Roman';
    entryRAAN.FontSize = 14;
    entryRAAN.Layout.Row = 7;
    entryRAAN.Layout.Column = 4;
    entryRAAN.Tooltip = 'Right ascension of the ascending node';
    entryRAAN.Limits = [0, 360]; % Limit RAAN entry to 0 to 360 degrees
    entryRAAN.Value = simParams.initParams.Orbit.RAAN; % Default to 0

    % True anomaly
    lblTA = uilabel(subGl, 'Text', '$\nu$ [deg]', 'Interpreter', 'latex');
    lblTA.FontName = 'Times New Roman';
    lblTA.FontSize = 14;
    lblTA.Layout.Row = 8;
    lblTA.Layout.Column = 5;
    entryTA = uieditfield(subGl, 'numeric');
    entryTA.FontName = 'Times New Roman';
    entryTA.FontSize = 14;
    entryTA.Layout.Row = 8;
    entryTA.Layout.Column = 4;
    entryTA.Tooltip = 'True anomaly';
    entryTA.Limits = [0, 360]; % Limit true anomaly entry to 0 to 360 degrees
    entryTA.Value = simParams.initParams.Orbit.trueAnomaly; % Default to 0

    %% Defining sliders
    % Altitude slider
    sliderAltitude = uislider(subGl);
    sliderAltitude.Layout.Row = 3;
    sliderAltitude.Layout.Column = 6;
    sliderAltitude.Limits = [100e3, 1000e3]; % 100 km to 1,000 km
    sliderAltitude.Value = simParams.initParams.Orbit.altitude; % Default to 400 km altitude
    sliderAltitude.MajorTicks = 100e3:300e3:1000e3; % Major ticks every 300 km
    sliderAltitude.MajorTickLabels = {'100 km', '400 km', '700 km', '1,000 km'};
    sliderAltitude.ValueChangingFcn = @(src, event) updateEntry(src, event, entryAltitude);
    entryAltitude.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderAltitude);

    % Eccentricity slider
    sliderEcc = uislider(subGl);
    sliderEcc.Layout.Row = 4;
    sliderEcc.Layout.Column = 6;
    sliderEcc.Limits = [0, 0.9]; % 0 to 0.9
    sliderEcc.Value = simParams.initParams.Orbit.eccentricity; % Default to 0
    sliderEcc.MajorTicks = 0:0.3:0.9; % Major ticks every 0.3
    sliderEcc.MajorTickLabels = {'0', '0.3', '0.6', '0.9'};
    sliderEcc.ValueChangingFcn = @(src, event) updateEntry(src, event, entryEcc);
    entryEcc.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderEcc);

    % Inclination slider
    sliderInc = uislider(subGl);
    sliderInc.Layout.Row = 5;
    sliderInc.Layout.Column = 6;
    sliderInc.Limits = [0, 180]; % 0 to 180 degrees
    sliderInc.Value = simParams.initParams.Orbit.inclination; % Default to 0
    sliderInc.MajorTicks = 0:60:180; % Major ticks every 60 degrees
    sliderInc.MajorTickLabels = {'0°', '60°', '120°', '180°'};
    sliderInc.ValueChangingFcn = @(src, event) updateEntry(src, event, entryInc);
    entryInc.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderInc);

    % Argument of periapse slider
    sliderAOP = uislider(subGl);
    sliderAOP.Layout.Row = 6;
    sliderAOP.Layout.Column = 6;
    sliderAOP.Limits = [0, 360]; % 0 to 360 degrees
    sliderAOP.Value = simParams.initParams.Orbit.argPeriapse; % Default to 0
    sliderAOP.MajorTicks = 0:90:360; % Major ticks every 90 degrees
    sliderAOP.MajorTickLabels = {'0°', '90°', '180°', '270°', '360°'};
    sliderAOP.ValueChangingFcn = @(src, event) updateEntry(src, event, entryAOP);
    entryAOP.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderAOP);

    % RAAN slider
    sliderRAAN = uislider(subGl);
    sliderRAAN.Layout.Row = 7;
    sliderRAAN.Layout.Column = 6;
    sliderRAAN.Limits = [0, 360]; % 0 to 360 degrees
    sliderRAAN.Value = simParams.initParams.Orbit.RAAN; % Default to 0
    sliderRAAN.MajorTicks = 0:90:360; % Major ticks every 90 degrees
    sliderRAAN.MajorTickLabels = {'0°', '90°', '180°', '270°', '360°'};
    sliderRAAN.ValueChangingFcn = @(src, event) updateEntry(src, event, entryRAAN);
    entryRAAN.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderRAAN);

    % True anomaly slider
    sliderTA = uislider(subGl);
    sliderTA.Layout.Row = 8;
    sliderTA.Layout.Column = 6;
    sliderTA.Limits = [0, 360]; % 0 to 360 degrees
    sliderTA.Value = simParams.initParams.Orbit.trueAnomaly; % Default to 0
    sliderTA.MajorTicks = 0:90:360; % Major ticks every 90 degrees
    sliderTA.MajorTickLabels = {'0°', '90°', '180°', '270°', '360°'};
    sliderTA.ValueChangingFcn = @(src, event) updateEntry(src, event, entryTA);
    entryTA.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderTA);

    function updateEntry(~, event, entry)
        value = event.Value;
        entry.Value = value;

        % Update the orbit plot in real-time as parameters change
        a = entryStruct.altitude.Value + earthRadius; % Convert altitude to semi-major axis
        e = entryStruct.eccentricity.Value;
        i = deg2rad(entryStruct.inclination.Value); % Convert to radians
        RAAN = deg2rad(entryStruct.RAAN.Value); % Convert to radians
        AOP = deg2rad(entryStruct.argPeriapse.Value); % Convert to radians
        nu = deg2rad(linspace(entryStruct.trueAnomaly.Value, entryStruct.trueAnomaly.Value + 360, 500)); % True anomaly from 0 to 360 degrees
        
        % Construct orbit in orbital frame
        r = a * (1 - e^2) ./ (1 + e * cos(nu)); % Orbital radius
        r_orbital = [r .* cos(nu); r .* sin(nu); zeros(size(nu))]; % Orbit in orbital plane

        % Creating transformation matrix from orbital frame to inertial frame
        R3_Omega = [cos(RAAN), -sin(RAAN), 0; 
                    sin(RAAN), cos(RAAN), 0; 
                    0, 0, 1]; % Rotation about Z-axis by RAAN

        R1_i = [1, 0, 0;
                 0, cos(i), -sin(i);
                 0, sin(i), cos(i)]; % Rotation about X-axis by inclination

        R3_omega = [cos(AOP), -sin(AOP), 0; 
                    sin(AOP), cos(AOP), 0; 
                    0, 0, 1]; % Rotation about Z-axis by argument of periapse

        R = R3_Omega * R1_i * R3_omega; % Combined rotation matrix

        % Rotate orbit from orbital frame to inertial frame
        r_inertial = R * r_orbital;
        x_orbit = r_inertial(1, :);
        y_orbit = r_inertial(2, :);
        z_orbit = r_inertial(3, :);
        
        % Plotting
        set(orbitLine, 'XData', x_orbit, 'YData', y_orbit, 'ZData', z_orbit); % Update the orbit plot with new coordinates
        set(orbitPosition, 'XData', x_orbit(1), 'YData', y_orbit(1), 'ZData', z_orbit(1)); % Update the position marker

        ax.XLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
        ax.YLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
        ax.ZLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit

    end

    function updateSlider(~, event, slider)
        value = event.Value;
        slider.Value = value;

        % Update the orbit plot in real-time as parameters change
        a = entryStruct.altitude.Value + earthRadius; % Convert altitude to semi-major axis
        e = entryStruct.eccentricity.Value;
        i = deg2rad(entryStruct.inclination.Value); % Convert to radians
        RAAN = deg2rad(entryStruct.RAAN.Value); % Convert to radians
        AOP = deg2rad(entryStruct.argPeriapse.Value); % Convert to radians
        nu = deg2rad(linspace(entryStruct.trueAnomaly.Value, entryStruct.trueAnomaly.Value + 360, 500)); % True anomaly from 0 to 360 degrees
        
        % Construct orbit in orbital frame
        r = a * (1 - e^2) ./ (1 + e * cos(nu)); % Orbital radius
        r_orbital = [r .* cos(nu); r .* sin(nu); zeros(size(nu))]; % Orbit in orbital plane

        % Creating transformation matrix from orbital frame to inertial frame
        R3_Omega = [cos(RAAN), -sin(RAAN), 0; 
                    sin(RAAN), cos(RAAN), 0; 
                    0, 0, 1]; % Rotation about Z-axis by RAAN

        R1_i = [1, 0, 0;
                 0, cos(i), -sin(i);
                 0, sin(i), cos(i)]; % Rotation about X-axis by inclination

        R3_omega = [cos(AOP), -sin(AOP), 0; 
                    sin(AOP), cos(AOP), 0; 
                    0, 0, 1]; % Rotation about Z-axis by argument of periapse

        R = R3_Omega * R1_i * R3_omega; % Combined rotation matrix

        % Rotate orbit from orbital frame to inertial frame
        r_inertial = R * r_orbital;
        x_orbit = r_inertial(1, :);
        y_orbit = r_inertial(2, :);
        z_orbit = r_inertial(3, :);
        
        % Plotting
        set(orbitLine, 'XData', x_orbit, 'YData', y_orbit, 'ZData', z_orbit); % Update the orbit plot with new coordinates
        set(orbitPosition, 'XData', x_orbit(1), 'YData', y_orbit(1), 'ZData', z_orbit(1)); % Update the position marker

        ax.XLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
        ax.YLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
        ax.ZLim = [-1.1, 1.1] * (a * (1 + e)); % Set limits to show the entire orbit
    end

    %% Saving logic
    entryStruct = struct(...
        'altitude', entryAltitude, ...
        'eccentricity', entryEcc, ...
        'inclination', entryInc, ...
        'argPeriapse', entryAOP, ...
        'RAAN', entryRAAN, ...
        'trueAnomaly', entryTA ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose();

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.initParams.Orbit.altitude = entryStruct.altitude.Value;
        simParams.initParams.Orbit.semiMajorAxis = entryStruct.altitude.Value + earthRadius/1000; % Convert altitude to semi-major axis [km]
        simParams.initParams.Orbit.eccentricity = entryStruct.eccentricity.Value;
        simParams.initParams.Orbit.inclination = entryStruct.inclination.Value;
        simParams.initParams.Orbit.argPeriapse = entryStruct.argPeriapse.Value;
        simParams.initParams.Orbit.RAAN = entryStruct.RAAN.Value;
        simParams.initParams.Orbit.trueAnomaly = entryStruct.trueAnomaly.Value;

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update main window title with asterisk to indicate unsaved changes
        simName = simParams.initParams.Simulation.name;
        subFigHandle.Name = [simName, ' *']; % Update the window title to include the simulation name and asterisk for unsaved changes

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end
end

function setAttitudeParameters(subFigHandle)
   salmonColor = [1, 0.4941, 0.4392];

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Create the New Window
    subFig = uifigure('Name', 'Set Orbital Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    % subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        45,...
        45,...
        45,...
        45,...
        45,...
        45,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        100,...
        160,...
        '1x',...
        'fit',...
        'fit',...
        300,...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Orbital Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 7];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining figure to show current orientation
    % Axes for orientation plot
    ax = uiaxes(subGl);
    ax.Layout.Row = [3, 8];
    ax.Layout.Column = [1, 2];
    ax.Title.String = 'Current Orientation';
    ax.Title.FontName = 'Times New Roman';
    ax.Title.FontSize = 14;
    ax.DataAspectRatio = [1, 1, 1]; % Equal aspect ratio for all axes
    ax.XLim = [-1.0, 1.0];
    ax.YLim = [-1.0, 1.0];
    ax.ZLim = [-1.0, 1.0];
    % Turning axes off for cleaner look since we're just showing orientation
    ax.XColor = 'none';
    ax.YColor = 'none';
    ax.ZColor = 'none';
    % Adding a faint grid for reference
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.ZGrid = 'on';
    hold(ax, "on");

    % Plotting inertial axes
    i_1 = [1, 0, 0]'; % Inertial X-axis
    i_2 = [0, 1, 0]'; % Inertial Y-axis
    i_3 = [0, 0, 1]'; % Inertial Z-axis
    inertialVectors = [i_1, i_2, i_3];
    [inertialArrows, inertialTexts] = draw_frame(inertialVectors, 'k', ax, 'i'); % Plot inertial frame in black

    % Plotting orbital axes
    a = simParams.initParams.Orbit.semiMajorAxis; % m
    e = simParams.initParams.Orbit.eccentricity;
    i = deg2rad(simParams.initParams.Orbit.inclination); % Convert to radians
    RAAN = deg2rad(simParams.initParams.Orbit.RAAN); % Convert to radians
    AOP = deg2rad(simParams.initParams.Orbit.argPeriapse); % Convert to radians
    nu = deg2rad(linspace(simParams.initParams.Orbit.trueAnomaly, simParams.initParams.Orbit.trueAnomaly + 360, 500)); % True anomaly from 0 to 360 degrees
    
    % Creating transformation matrix from orbital frame to inertial frame (R_IO)
    R3_Omega = [cos(RAAN), -sin(RAAN), 0; 
                sin(RAAN), cos(RAAN), 0; 
                0, 0, 1]; % Rotation about Z-axis by RAAN

    R1_i = [1, 0, 0;
            0, cos(i), -sin(i);
            0, sin(i), cos(i)]; % Rotation about X-axis by inclination

    R3_omega = [cos(AOP), -sin(AOP), 0; 
                sin(AOP), cos(AOP), 0; 
                0, 0, 1]; % Rotation about Z-axis by argument of periapse

    R = R3_Omega * R1_i * R3_omega; % Combined rotation matrix

    o_3 = R(:, 3); % Orbital Z-axis (angular momentum vector)
    o_1 = R(:, 1)*cos(nu(1)) + R(:, 2)*sin(nu(1)); % Orbital X-axis (points radially outward towards current position)
    o_2 = cross(o_3, o_1); % Orbital Y-axis (completes right-handed system)
    orbitalVectors = [o_1, o_2, o_3];
    [orbitalArrows, orbitalTexts] = draw_frame(orbitalVectors, 'b', ax, 'o'); % Plot orbital frame in blue

    % Plotting body axes defined with respect to orbital frame
    v_I = R * [-sin(nu(1)); e + cos(nu(1)); 0]; % Velocity vector in inertial frame
    ref1 = v_I / norm(v_I); % Reference X-axis points along velocity vector (in direction of motion)
    ref3 = -o_3; % Reference Z-axis aligned with orbital angular momentum vector
    ref2 = cross(ref3, ref1); % Reference Y-axis completes right-handed system
    R_Ref2I = [ref1, ref2, ref3];

    % Rotating into body frame using attitude angles (roll, pitch, yaw)
    roll = deg2rad(simParams.initParams.Attitude.roll);
    pitch = deg2rad(simParams.initParams.Attitude.pitch);
    yaw = deg2rad(simParams.initParams.Attitude.yaw);

    % Rotation about Z (yaw)
    Rz_yaw = [cos(yaw), -sin(yaw), 0;
              sin(yaw), cos(yaw), 0;
              0, 0, 1];

    % Rotation about Y (pitch)
    Ry_pitch = [cos(pitch), 0, sin(pitch);
                0, 1, 0;
                -sin(pitch), 0, cos(pitch)];

    % Rotation about X (roll)
    Rx_roll = [1, 0, 0;
               0, cos(roll), -sin(roll);
               0, sin(roll), cos(roll)];

    % Combined rotation from reference frame to body frame
    R_B2Ref = Rz_yaw * Ry_pitch * Rx_roll;

    % Final rotation from inertial frame to body frame
    R_IB = R_Ref2I * R_B2Ref;

    b_1 = R_IB(:, 1); % Body X-axis
    b_2 = R_IB(:, 2); % Body Y-axis
    b_3 = R_IB(:, 3); % Body Z-axis
    bodyVectors = [b_1, b_2, b_3];
    [bodyArrows, bodyTexts] = draw_frame(bodyVectors, 'r', ax, 'b'); % Plot body frame in red

    view(ax, 30, 30)

    hold(ax, "off");

    function [arrows, texts] = draw_frame(vectors, color, parent, symbol)
        arrows = cell(1, 3);
        texts = cell(1, 3);
        for i = 1:3
            arrows{i} = quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), color, 'LineWidth', 2, Parent=parent);
            arrows{i}.AutoScale = 'off'; % Disable auto-scaling to keep arrows at correct length
            if nargin == 4
                texts{i} = text(vectors(1, i), vectors(2, i), vectors(3, i), [' ', symbol, '_', num2str(i)], 'FontSize', 16, 'FontName', 'Times New Roman', Parent=parent, Color=color);
            end
        end
    end

    %% Defining entries
    % Roll
    lblRoll = uilabel(subGl, 'Text', '$\phi$ [deg]', 'Interpreter', 'latex');
    lblRoll.FontName = 'Times New Roman';
    lblRoll.FontSize = 14;
    lblRoll.Layout.Row = 3;
    lblRoll.Layout.Column = 5;
    entryRoll = uieditfield(subGl, 'numeric');
    entryRoll.FontName = 'Times New Roman';
    entryRoll.FontSize = 14;
    entryRoll.Layout.Row = 3;
    entryRoll.Layout.Column = 4;
    entryRoll.Tooltip = 'Roll angle';
    entryRoll.Limits = [-180, 180]; % Limit roll entry to -180 to 180 degrees
    entryRoll.Value = simParams.initParams.Attitude.roll; % Default to 0

    % Pitch
    lblPitch = uilabel(subGl, 'Text', '$\theta$ [deg]', 'Interpreter', 'latex');
    lblPitch.FontName = 'Times New Roman';
    lblPitch.FontSize = 14;
    lblPitch.Layout.Row = 4;
    lblPitch.Layout.Column = 5;
    entryPitch = uieditfield(subGl, 'numeric');
    entryPitch.FontName = 'Times New Roman';
    entryPitch.FontSize = 14;
    entryPitch.Layout.Row = 4;
    entryPitch.Layout.Column = 4;
    entryPitch.Tooltip = 'Pitch angle';
    entryPitch.Limits = [-180, 180]; % Limit pitch entry to -180 to 180 degrees
    entryPitch.Value = simParams.initParams.Attitude.pitch; % Default to 0

    % Yaw
    lblYaw = uilabel(subGl, 'Text', '$\psi$ [deg]', 'Interpreter', 'latex');
    lblYaw.FontName = 'Times New Roman';
    lblYaw.FontSize = 14;
    lblYaw.Layout.Row = 5;
    lblYaw.Layout.Column = 5;
    entryYaw = uieditfield(subGl, 'numeric');
    entryYaw.FontName = 'Times New Roman';
    entryYaw.FontSize = 14;
    entryYaw.Layout.Row = 5;
    entryYaw.Layout.Column = 4;
    entryYaw.Tooltip = 'Yaw angle';
    entryYaw.Limits = [-180, 180]; % Limit yaw entry to -180 to 180 degrees
    entryYaw.Value = simParams.initParams.Attitude.yaw; % Default to 0

    % Roll rate
    lblRollRate = uilabel(subGl, 'Text', '$\dot{\phi}$ [deg/s]', 'Interpreter', 'latex');
    lblRollRate.FontName = 'Times New Roman';
    lblRollRate.FontSize = 14;
    lblRollRate.Layout.Row = 6;
    lblRollRate.Layout.Column = 5;
    entryRollRate = uieditfield(subGl, 'numeric');
    entryRollRate.FontName = 'Times New Roman';
    entryRollRate.FontSize = 14;
    entryRollRate.Layout.Row = 6;
    entryRollRate.Layout.Column = 4;
    entryRollRate.Tooltip = 'Roll rate';
    entryRollRate.Limits = [-180, 180]; % Limit roll rate entry to -180 to 180 deg/s
    entryRollRate.Value = simParams.initParams.Attitude.rollRate; % Default to 0

    % Pitch rate
    lblPitchRate = uilabel(subGl, 'Text', '$\dot{\theta}$ [deg/s]', 'Interpreter', 'latex');
    lblPitchRate.FontName = 'Times New Roman';
    lblPitchRate.FontSize = 14;
    lblPitchRate.Layout.Row = 7;
    lblPitchRate.Layout.Column = 5;
    entryPitchRate = uieditfield(subGl, 'numeric');
    entryPitchRate.FontName = 'Times New Roman';
    entryPitchRate.FontSize = 14;
    entryPitchRate.Layout.Row = 7;
    entryPitchRate.Layout.Column = 4;
    entryPitchRate.Tooltip = 'Pitch rate';
    entryPitchRate.Limits = [-180, 180]; % Limit pitch rate entry to -180 to 180 deg/s
    entryPitchRate.Value = simParams.initParams.Attitude.pitchRate; % Default to 0

    % Yaw rate
    lblYawRate = uilabel(subGl, 'Text', '$\dot{\psi}$ [deg/s]', 'Interpreter', 'latex');
    lblYawRate.FontName = 'Times New Roman';
    lblYawRate.FontSize = 14;
    lblYawRate.Layout.Row = 8;
    lblYawRate.Layout.Column = 5;
    entryYawRate = uieditfield(subGl, 'numeric');
    entryYawRate.FontName = 'Times New Roman';
    entryYawRate.FontSize = 14;
    entryYawRate.Layout.Row = 8;
    entryYawRate.Layout.Column = 4;
    entryYawRate.Tooltip = 'Yaw rate';
    entryYawRate.Limits = [-180, 180]; % Limit yaw rate entry to -180 to 180 deg/s
    entryYawRate.Value = simParams.initParams.Attitude.yawRate; % Default to 0

    %% Defining sliders
    % Roll slider
    sliderRoll = uislider(subGl);
    sliderRoll.Layout.Row = 3;
    sliderRoll.Layout.Column = 6;
    sliderRoll.Limits = [-180, 180]; % -180 to 180 degrees
    sliderRoll.Value = simParams.initParams.Attitude.roll; % Default to 0
    sliderRoll.MajorTicks = -180:90:180; % Major ticks every 90 degrees
    sliderRoll.MajorTickLabels = {'-180°', '-90°', '0°', '90°', '180°'};
    sliderRoll.ValueChangingFcn = @(src, event) updateEntry(src, event, entryRoll, 'Roll');
    entryRoll.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderRoll, 'Roll');

    % Pitch slider
    sliderPitch = uislider(subGl);
    sliderPitch.Layout.Row = 4;
    sliderPitch.Layout.Column = 6;
    sliderPitch.Limits = [-180, 180]; % -180 to 180 degrees
    sliderPitch.Value = simParams.initParams.Attitude.pitch; % Default to 0
    sliderPitch.MajorTicks = -180:90:180; % Major ticks every 90 degrees
    sliderPitch.MajorTickLabels = {'-180°', '-90°', '0°', '90°', '180°'};
    sliderPitch.ValueChangingFcn = @(src, event) updateEntry(src, event, entryPitch, 'Pitch');
    entryPitch.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderPitch, 'Pitch');

    % Yaw slider
    sliderYaw = uislider(subGl);
    sliderYaw.Layout.Row = 5;
    sliderYaw.Layout.Column = 6;
    sliderYaw.Limits = [-180, 180]; % -180 to 180 degrees
    sliderYaw.Value = simParams.initParams.Attitude.yaw; % Default to 0
    sliderYaw.MajorTicks = -180:90:180; % Major ticks every 90 degrees
    sliderYaw.MajorTickLabels = {'-180°', '-90°', '0°', '90°', '180°'};
    sliderYaw.ValueChangingFcn = @(src, event) updateEntry(src, event, entryYaw, 'Yaw');
    entryYaw.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderYaw, 'Yaw');

    % Roll rate slider
    sliderRollRate = uislider(subGl);
    sliderRollRate.Layout.Row = 6;
    sliderRollRate.Layout.Column = 6;
    sliderRollRate.Limits = [-180, 180]; % -180 to 180 deg/s
    sliderRollRate.Value = simParams.initParams.Attitude.rollRate; % Default to 0
    sliderRollRate.MajorTicks = -180:90:180; % Major ticks every 90 deg/s
    sliderRollRate.MajorTickLabels = {'-180°/s', '-90°/s', '0°/s', '90°/s', '180°/s'};
    sliderRollRate.ValueChangingFcn = @(src, event) updateEntry(src, event, entryRollRate, 'RollRate');
    entryRollRate.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderRollRate, 'RollRate');

    % Pitch rate slider
    sliderPitchRate = uislider(subGl);
    sliderPitchRate.Layout.Row = 7;
    sliderPitchRate.Layout.Column = 6;
    sliderPitchRate.Limits = [-180, 180]; % -180 to 180 deg/s
    sliderPitchRate.Value = simParams.initParams.Attitude.pitchRate; % Default to 0
    sliderPitchRate.MajorTicks = -180:90:180; % Major ticks every 90 deg/s
    sliderPitchRate.MajorTickLabels = {'-180°/s', '-90°/s', '0°/s', '90°/s', '180°/s'};
    sliderPitchRate.ValueChangingFcn = @(src, event) updateEntry(src, event, entryPitchRate, 'PitchRate');
    entryPitchRate.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderPitchRate, 'PitchRate');

    % Yaw rate slider
    sliderYawRate = uislider(subGl);
    sliderYawRate.Layout.Row = 8;
    sliderYawRate.Layout.Column = 6;
    sliderYawRate.Limits = [-180, 180]; % -180 to 180 deg/s
    sliderYawRate.Value = simParams.initParams.Attitude.yawRate; % Default to 0
    sliderYawRate.MajorTicks = -180:90:180; % Major ticks every 90 deg/s
    sliderYawRate.MajorTickLabels = {'-180°/s', '-90°/s', '0°/s', '90°/s', '180°/s'};
    sliderYawRate.ValueChangingFcn = @(src, event) updateEntry(src, event, entryYawRate, 'YawRate');
    entryYawRate.ValueChangedFcn = @(src, event) updateSlider(src, event, sliderYawRate, 'YawRate');

    function updateEntry(~, event, entry, rotationType)
        value = event.Value;
        entry.Value = value;

        if rotationType == "RollRate" || rotationType == "PitchRate" || rotationType == "YawRate"
            % If we're updating a rate, we don't need to update the orientation plot since it only depends on the angles
            return;
        end

        % Update the orientation plot in real-time as parameters change
        roll = deg2rad(entryStruct.roll.Value);
        pitch = deg2rad(entryStruct.pitch.Value);
        yaw = deg2rad(entryStruct.yaw.Value);

        % Rotation about Z (yaw)
        Rz_yaw = [cos(yaw), -sin(yaw), 0;
                  sin(yaw), cos(yaw), 0;
                  0, 0, 1];

        % Rotation about Y (pitch)
        Ry_pitch = [cos(pitch), 0, sin(pitch);
                    0, 1, 0;
                    -sin(pitch), 0, cos(pitch)];

        % Rotation about X (roll)
        Rx_roll = [1, 0, 0;
                   0, cos(roll), -sin(roll);
                   0, sin(roll), cos(roll)];

        % Combined rotation from reference frame to body frame
        R_B2Ref = Rz_yaw * Ry_pitch * Rx_roll;

        % Final rotation from inertial frame to body frame
        R_IB = R_Ref2I * R_B2Ref;

        b_1 = R_IB(:, 1); % Body X-axis
        b_2 = R_IB(:, 2); % Body Y-axis
        b_3 = R_IB(:, 3); % Body Z-axis
        bodyVectors = [b_1, b_2, b_3];
        
        % Update body frame arrows
        for i = 1:3
            set(bodyArrows{i}, 'UData', bodyVectors(1, i), 'VData', bodyVectors(2, i), 'WData', bodyVectors(3, i));
            set(bodyTexts{i}, 'Position', bodyVectors(:, i)');
        end

    end

    function updateSlider(~, event, slider, rotationType)
        value = event.Value;
        slider.Value = value;

        if rotationType == "RollRate" || rotationType == "PitchRate" || rotationType == "YawRate"
            % If we're updating a rate, we don't need to update the orientation plot since it only depends on the angles
            return;
        end

        % Update the orientation plot in real-time as parameters change
        roll = deg2rad(entryStruct.roll.Value);
        pitch = deg2rad(entryStruct.pitch.Value);
        yaw = deg2rad(entryStruct.yaw.Value);

        % Rotation about Z (yaw)
        Rz_yaw = [cos(yaw), -sin(yaw), 0;
                  sin(yaw), cos(yaw), 0;
                  0, 0, 1];

        % Rotation about Y (pitch)
        Ry_pitch = [cos(pitch), 0, sin(pitch);
                    0, 1, 0;
                    -sin(pitch), 0, cos(pitch)];

        % Rotation about X (roll)
        Rx_roll = [1, 0, 0;
                   0, cos(roll), -sin(roll);
                   0, sin(roll), cos(roll)];

        % Combined rotation from reference frame to body frame
        R_B2Ref = Rz_yaw * Ry_pitch * Rx_roll;

        % Final rotation from inertial frame to body frame
        R_IB = R_Ref2I * R_B2Ref;

        b_1 = R_IB(:, 1); % Body X-axis
        b_2 = R_IB(:, 2); % Body Y-axis
        b_3 = R_IB(:, 3); % Body Z-axis
        bodyVectors = [b_1, b_2, b_3];
        
        % Update body frame arrows
        for i = 1:3
            set(bodyArrows{i}, 'UData', bodyVectors(1, i), 'VData', bodyVectors(2, i), 'WData', bodyVectors(3, i));
            set(bodyTexts{i}, 'Position', bodyVectors(:, i)');
        end

    end

    %% Saving logic
    entryStruct = struct(...
        'roll', entryRoll, ...
        'pitch', entryPitch, ...
        'yaw', entryYaw, ...
        'rollRate', entryRollRate, ...
        'pitchRate', entryPitchRate, ...
        'yawRate', entryYawRate ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose;

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.initParams.Attitude.roll = entryStruct.roll.Value;
        simParams.initParams.Attitude.pitch = entryStruct.pitch.Value;
        simParams.initParams.Attitude.yaw = entryStruct.yaw.Value;
        simParams.initParams.Attitude.rollRate = entryStruct.rollRate.Value;
        simParams.initParams.Attitude.pitchRate = entryStruct.pitchRate.Value;
        simParams.initParams.Attitude.yawRate = entryStruct.yawRate.Value;

        % Converting Euler angles to quaternions for internal use in dynamics calculations
        % Euler angles encode the orientation between orbital frame and body frame, so must get DCM from body to inertial
        Rz_yaw = [cosd(simParams.initParams.Attitude.yaw), -sind(simParams.initParams.Attitude.yaw), 0;
                  sind(simParams.initParams.Attitude.yaw), cosd(simParams.initParams.Attitude.yaw), 0;
                  0, 0, 1];

        % Rotation about Y (pitch)
        Ry_pitch = [cosd(simParams.initParams.Attitude.pitch), 0, sind(simParams.initParams.Attitude.pitch);
                    0, 1, 0;
                    -sind(simParams.initParams.Attitude.pitch), 0, cosd(simParams.initParams.Attitude.pitch)];

        % Rotation about X (roll)
        Rx_roll = [1, 0, 0;
                   0, cosd(simParams.initParams.Attitude.roll), -sind(simParams.initParams.Attitude.roll);
                   0, sind(simParams.initParams.Attitude.roll), cosd(simParams.initParams.Attitude.roll)];

        R_B2Ref = Rz_yaw * Ry_pitch * Rx_roll; % Rotation from reference frame to body frame
        R_IB = R_Ref2I * R_B2Ref; % Rotation from inertial frame to body frame

        simParams.initParams.Attitude.Quaternion = dcm2quat(R_IB); % Convert DCM to quaternions

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update main window title with asterisk to indicate unsaved changes
        simName = simParams.initParams.Simulation.name;
        subFigHandle.Name = [simName, ' *']; % Update the window title to include the simulation name and asterisk for unsaved changes

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end
end

function setEnvironmentalParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Create the New Window
    subFig = uifigure('Name', 'Set Environmental Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    % subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        20,...
        20,...
        20,...
        20,...
        20,...
        20,...
        '1x',...
        20,...
        20,...
        20,...
        20,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        'fit',...
        '1x',...
        100,...
        'fit',...
        '1x',...
        100,...
        'fit',...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Environmental Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 8];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 15; % Bottom row
    btnSaveClose.Layout.Column = [1, 3]; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Defining entries
    % Parameters label
    lblParams = uilabel(subGl, 'Text', 'Parameters:');
    lblParams.FontName = 'Times New Roman';
    lblParams.FontSize = 16;
    lblParams.FontWeight = 'bold';
    lblParams.Layout.Row = 3;
    lblParams.Layout.Column = [3, 4];

    % Seconds of day
    lblSOD = uilabel(subGl, 'Text', 'Seconds of day [s]');
    lblSOD.FontName = 'Times New Roman';
    lblSOD.FontSize = 14;
    lblSOD.Layout.Row = 4;
    lblSOD.Layout.Column = 4;
    entrySOD = uieditfield(subGl, 'numeric');
    entrySOD.FontName = 'Times New Roman';
    entrySOD.FontSize = 14;
    entrySOD.Layout.Row = 4;
    entrySOD.Layout.Column = 3;
    entrySOD.Value = simParams.initParams.Environment.secondsOfDay; % Default to 0

    % Day of year
    lblDOY = uilabel(subGl, 'Text', 'Day of year');
    lblDOY.FontName = 'Times New Roman';
    lblDOY.FontSize = 14;
    lblDOY.Layout.Row = 5;
    lblDOY.Layout.Column = 4;
    entryDOY = uieditfield(subGl, 'numeric');
    entryDOY.FontName = 'Times New Roman';
    entryDOY.FontSize = 14;
    entryDOY.Layout.Row = 5;
    entryDOY.Layout.Column = 3;
    entryDOY.Value = simParams.initParams.Environment.dayOfYear; % Default to 1

    % 81-day average F10.7 solar flux
    lblF107 = uilabel(subGl, 'Text', '81-day average F10.7 solar flux');
    lblF107.FontName = 'Times New Roman';
    lblF107.FontSize = 14;
    lblF107.Layout.Row = 6;
    lblF107.Layout.Column = 4;
    entryF107 = uieditfield(subGl, 'numeric');
    entryF107.FontName = 'Times New Roman';
    entryF107.FontSize = 14;
    entryF107.Layout.Row = 6;
    entryF107.Layout.Column = 3;
    entryF107.Value = simParams.initParams.Environment.F107Average; % Default to 65

    % Daily F10.7 solar flux
    lblF107D = uilabel(subGl, 'Text', 'Daily F10.7 solar flux');
    lblF107D.FontName = 'Times New Roman';
    lblF107D.FontSize = 14;
    lblF107D.Layout.Row = 7;
    lblF107D.Layout.Column = 4;
    entryF107D = uieditfield(subGl, 'numeric');
    entryF107D.FontName = 'Times New Roman';
    entryF107D.FontSize = 14;
    entryF107D.Layout.Row = 7;
    entryF107D.Layout.Column = 3;
    entryF107D.Value = simParams.initParams.Environment.F107Daily; % Default to 65

    % Magnetic indices
    % Dividing grid section into 7 columns 
    magneticGl = uigridlayout(subGl, [1, 8]);
    magneticGl.RowHeight = 20;
    magneticGl.ColumnWidth = {
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        'fit',...
    };
    magneticGl.Layout.Row = 8;
    magneticGl.Layout.Column = [3, 4];
    magneticGl.Padding = [0, 0, 0, 0]; % Remove padding for tighter layout

    lblMag = uilabel(magneticGl, 'Text', 'Magnetic indices');
    lblMag.FontName = 'Times New Roman';
    lblMag.FontSize = 14;
    lblMag.Layout.Row = 1;
    lblMag.Layout.Column = 8;

    entryMags = gobjects(1, 7); % Preallocate array for magnetic index entries
    for i = 1 : 1 : 7
        entry = uieditfield(magneticGl, 'numeric');
        entry.FontName = 'Times New Roman';
        entry.FontSize = 14;
        entry.Layout.Row = 1;
        entry.Layout.Column = i;
        entry.Tooltip = sprintf('Magnetic index %d', i); % Add tooltip for clarity
        entry.ValueDisplayFormat = '%.1f'; % Format to 1 decimal place
        entry.HorizontalAlignment = 'center'; % Center align the text for better aesthetics
        entry.Value = simParams.initParams.Environment.magneticIndices(i); % Default to 0
        entryMags(i) = entry; % Store handle for later access
    end

    % Gas-surface interaction model
    lblGSI = uilabel(subGl, 'Text', 'Gas-surface interaction model');
    lblGSI.FontName = 'Times New Roman';
    lblGSI.FontSize = 14;
    lblGSI.Layout.Row = 4;
    lblGSI.Layout.Column = 7;
    entryGSI = uidropdown(subGl);
    entryGSI.Items = {'cook'};
    entryGSI.FontName = 'Times New Roman';
    entryGSI.FontSize = 14;
    entryGSI.Layout.Row = 4;
    entryGSI.Layout.Column = 6;
    entryGSI.Value = simParams.initParams.Environment.gasSurfaceInteractionModel; % Default to cook

    % Accommodation coefficient
    lblAccom = uilabel(subGl, 'Text', 'Accommodation coefficient');
    lblAccom.FontName = 'Times New Roman';
    lblAccom.FontSize = 14;
    lblAccom.Layout.Row = 5;
    lblAccom.Layout.Column = 7;
    entryAccom = uieditfield(subGl, 'numeric');
    entryAccom.FontName = 'Times New Roman';
    entryAccom.FontSize = 14;
    entryAccom.Layout.Row = 5;
    entryAccom.Layout.Column = 6;
    entryAccom.Value = simParams.initParams.Environment.accommodationCoefficient; % Default to 1.0

    % Wall temperature
    lblTemp = uilabel(subGl, 'Text', 'Wall temperature [K]');
    lblTemp.FontName = 'Times New Roman';
    lblTemp.FontSize = 14;
    lblTemp.Layout.Row = 6;
    lblTemp.Layout.Column = 7;
    entryTemp = uieditfield(subGl, 'numeric');
    entryTemp.FontName = 'Times New Roman';
    entryTemp.FontSize = 14;
    entryTemp.Layout.Row = 6;
    entryTemp.Layout.Column = 6;
    entryTemp.Value = simParams.initParams.Environment.wallTemperature; % Default to 300 K

    % Specular reflectivity
    lblSpec = uilabel(subGl, 'Text', 'Specular reflectivity');
    lblSpec.FontName = 'Times New Roman';
    lblSpec.FontSize = 14;
    lblSpec.Layout.Row = 7;
    lblSpec.Layout.Column = 7;
    entrySpec = uieditfield(subGl, 'numeric');
    entrySpec.FontName = 'Times New Roman';
    entrySpec.FontSize = 14;
    entrySpec.Layout.Row = 7;
    entrySpec.Layout.Column = 6;
    entrySpec.Value = simParams.initParams.Environment.specularReflectivity; % Default to 0.15

    % Diffuse reflectivity
    lblDiff = uilabel(subGl, 'Text', 'Diffuse reflectivity');
    lblDiff.FontName = 'Times New Roman';
    lblDiff.FontSize = 14;
    lblDiff.Layout.Row = 8;
    lblDiff.Layout.Column = 7;
    entryDiff = uieditfield(subGl, 'numeric');
    entryDiff.FontName = 'Times New Roman';
    entryDiff.FontSize = 14;
    entryDiff.Layout.Row = 8;
    entryDiff.Layout.Column = 6;
    entryDiff.Value = simParams.initParams.Environment.diffuseReflectivity; % Default to 0.25

    % Options
    lblOptions = uilabel(subGl, 'Text', 'Options:');
    lblOptions.FontName = 'Times New Roman';
    lblOptions.FontSize = 16;
    lblOptions.FontWeight = 'bold';
    lblOptions.Layout.Row = 10;
    lblOptions.Layout.Column = [3, 4];

    % Enable anomalous oxygen
    lblAO = uilabel(subGl, 'Text', 'Enable anomalous oxygen');
    lblAO.FontName = 'Times New Roman';
    lblAO.FontSize = 14;
    lblAO.Layout.Row = 11;
    lblAO.Layout.Column = 4;
    entryAO = uidropdown(subGl);
    entryAO.Items = {'true', 'false'};
    entryAO.FontName = 'Times New Roman';
    entryAO.FontSize = 14;
    entryAO.Layout.Row = 11;
    entryAO.Layout.Column = 3;
    if simParams.initParams.Environment.enableAnomalousOxygen == true
        entryAO.Value = 'true'; % Default to false
    else
        entryAO.Value = 'false';
    end

    % Enable shadow analysis
    lblShadow = uilabel(subGl, 'Text', 'Enable shadow analysis');
    lblShadow.FontName = 'Times New Roman';
    lblShadow.FontSize = 14;
    lblShadow.Layout.Row = 12;
    lblShadow.Layout.Column = 4;
    entryShadow = uidropdown(subGl);
    entryShadow.Items = {'true', 'false'};
    entryShadow.FontName = 'Times New Roman';
    entryShadow.FontSize = 14;
    entryShadow.Layout.Row = 12;
    entryShadow.Layout.Column = 3;
    if simParams.initParams.Environment.enableShadowAnalysis == true
        entryShadow.Value = 'true'; % Default to true
    else
        entryShadow.Value = 'false';
    end

    % Enable solar radiation pressure
    lblSRP = uilabel(subGl, 'Text', 'Enable solar radiation pressure');
    lblSRP.FontName = 'Times New Roman';
    lblSRP.FontSize = 14;
    lblSRP.Layout.Row = 13;
    lblSRP.Layout.Column = 4;
    entrySRP = uidropdown(subGl);
    entrySRP.Items = {'true', 'false'};
    entrySRP.FontName = 'Times New Roman';
    entrySRP.FontSize = 14;
    entrySRP.Layout.Row = 13;
    entrySRP.Layout.Column = 3;
    if simParams.initParams.Environment.enableSolarRadiationPressure == true
        entrySRP.Value = 'true'; % Default to true
    else
        entrySRP.Value = 'false';
    end

    %% Saving logic
    entryStruct = struct(...
        'secondsOfDay', entrySOD, ...
        'dayOfYear', entryDOY, ...
        'F107', entryF107, ...
        'F107D', entryF107D, ...
        'magneticIndices', entryMags, ...
        'gasSurfaceModel', entryGSI, ...
        'accommodationCoefficient', entryAccom, ...
        'wallTemperature', entryTemp, ...
        'specularReflectivity', entrySpec, ...
        'diffuseReflectivity', entryDiff, ...
        'enableAnomalousOxygen', entryAO, ...
        'enableShadowAnalysis', entryShadow, ...
        'enableSolarRadiationPressure', entrySRP ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose();

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.initParams.Environment.secondsOfDay = entryStruct.secondsOfDay.Value;
        simParams.initParams.Environment.dayOfYear = entryStruct.dayOfYear.Value;
        simParams.initParams.Environment.F107Average = entryStruct.F107.Value;
        simParams.initParams.Environment.F107Daily = entryStruct.F107D.Value;
        simParams.initParams.Environment.magneticIndices = arrayfun(@(e) e.Value, entryStruct.magneticIndices);
        simParams.initParams.Environment.gasSurfaceModel = entryStruct.gasSurfaceModel.Value;
        simParams.initParams.Environment.accommodationCoefficient = entryStruct.accommodationCoefficient.Value;
        simParams.initParams.Environment.wallTemperature = entryStruct.wallTemperature.Value;
        simParams.initParams.Environment.specularReflectivity = entryStruct.specularReflectivity.Value;
        simParams.initParams.Environment.diffuseReflectivity = entryStruct.diffuseReflectivity.Value;
        simParams.initParams.Environment.enableAnomalousOxygen = strcmp(entryStruct.enableAnomalousOxygen.Value, 'true');
        simParams.initParams.Environment.enableShadowAnalysis = strcmp(entryStruct.enableShadowAnalysis.Value, 'true');
        simParams.initParams.Environment.enableSolarRadiationPressure = strcmp(entryStruct.enableSolarRadiationPressure.Value, 'true');

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        % Update main window title with asterisk to indicate unsaved changes
        simName = simParams.initParams.Simulation.name;
        subFigHandle.Name = ['New Simulation - ', simName, ' *']; % Update the window title to include the simulation name and asterisk for unsaved changes

        delete(subFig);
    end

end

function setControllerParameters(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Create the New Window
    subFig = uifigure('Name', 'Set Controller Parameters', ...
                      'Position', [500, 300, 800, 450]);

    % Prevent users from clicking out of figure
    % subFig.WindowStyle = 'modal';

    % Initializing grid layout
    subGl = uigridlayout(subFig);
    subGl.RowHeight = {
        30,... 
        '1x',...
        30,...
        30,...
        30,...
        30,...
        30,...
        30,...
        '1x',...
        30,...
    };

    subGl.ColumnWidth = {
        'fit',...
        100,...
        200,...
        'fit',...
        '1x',...
    };

    % Add a Title to the Sub-Window (Center of content area)
    lbl = uilabel(subGl, 'Text', 'Set Controller Parameters');
    lbl.FontName = 'Times New Roman';
    lbl.FontSize = 20;
    lbl.FontWeight = 'bold';
    lbl.Layout.Row = 1;
    lbl.Layout.Column = [1 5];
    lbl.HorizontalAlignment = 'center';

    %% Defining buttons
    % Save and close button
    btnSaveClose = uibutton(subGl, 'Text', 'Save and Close');
    btnSaveClose.FontName = 'Times New Roman';
    btnSaveClose.FontSize = 14;
    btnSaveClose.Layout.Row = 10; % Bottom row
    btnSaveClose.Layout.Column = 1; % Left column
    btnSaveClose.BackgroundColor = salmonColor;

    %% Choosing a function file
    lblFunc = uilabel(subGl, 'Text', 'Controller function file');
    lblFunc.FontName = 'Times New Roman';
    lblFunc.FontSize = 14;
    lblFunc.Layout.Row = 3;
    lblFunc.Layout.Column = 4;
    entryFunc = uieditfield(subGl, 'text');
    entryFunc.FontName = 'Times New Roman';
    entryFunc.FontSize = 14;
    entryFunc.Layout.Row = 3;
    entryFunc.Layout.Column = 3;
    entryFunc.Value = simParams.initParams.Controller.functionFile; % Default to a function file
    btnBrowse = uibutton(subGl, 'Text', 'Browse');
    btnBrowse.FontName = 'Times New Roman';
    btnBrowse.FontSize = 14;
    btnBrowse.Layout.Row = 3;
    btnBrowse.Layout.Column = 2;
    btnBrowse.ButtonPushedFcn = @(~,~) BrowseFunction();

    function BrowseFunction()
        controllerFunctionPath = './src/controllers/'; % Default to controllers directory
        [file, ~] = uigetfile(fullfile(controllerFunctionPath, '*.m'), 'Select Controller Function File');
        if isequal(file, 0)
            return; % User canceled file selection
        end
        entryFunc.Value = fullfile(controllerFunctionPath, file); % Set the selected file path in the entry field
    end

     %% Saving logic
     entryStruct = struct(...
        'functionFile', entryFunc ...
    );

    % Logic: Save all entries and close this window
    btnSaveClose.ButtonPushedFcn = @(~,~) SaveClose();

    function SaveClose()
        % Retrieve simParams from figure's guidata
        simParams = guidata(subFigHandle);

        % Saving parameters to simParams structure
        simParams.initParams.Controller.functionFile = entryStruct.functionFile.Value;

        % Creating function handle for controller
        funcFile = simParams.initParams.Controller.functionFile(1:end-2); % Removing .m extension
        arr = split(funcFile, '/');
        funcFile = arr{end}; % Get just the file name without path
        simParams.initParams.Controller.Func = str2func(funcFile);

        simParams.saved = false; % Mark as unsaved since parameters changed

        % Update main window title with asterisk to indicate unsaved changes
        simName = simParams.initParams.Simulation.name;
        subFigHandle.Name = ['New Simulation - ', simName, ' *']; % Update the window title to include the simulation name and asterisk for unsaved changes

        % Update the guidata with the modified simParams
        guidata(subFigHandle, simParams);

        delete(subFig);
    end

end