function runSimulation(subFigHandle)
% RUNSIMULATION Runs the satellite simulation using parameters stored in
% the figure's guidata. Integrates the equations of motion using ode45 and
% stores results back in guidata for later visualization.
%
% Parameters
% ----------
% subFigHandle : matlab.ui.Figure
%   Handle to the simulation configuration figure (contains simParams in
%   guidata)

    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    %% Extract simulation parameters
    simParams = guidata(subFigHandle);
    simParams = ensureSimulationDefaults(simParams);
    guidata(subFigHandle, simParams);

    %% Build initial orbital state
    % simParams stores altitude and semiMajorAxis in m
    a   = simParams.initParams.Orbit.semiMajorAxis; % [m]
    e   = simParams.initParams.Orbit.eccentricity;
    inc = deg2rad(simParams.initParams.Orbit.inclination);
    raan = deg2rad(simParams.initParams.Orbit.RAAN);
    aop  = deg2rad(simParams.initParams.Orbit.argPeriapse);
    ta   = deg2rad(simParams.initParams.Orbit.trueAnomaly);

    orbit = [a, e, inc, raan, aop, ta];
    muOrbit = getEnvValue(simParams, 'mu', 3.986004e14);

    % Compute initial position and velocity in ECI frame
    [r_i_eci, v_i_eci] = keplerian_to_ijk_safe(orbit(1), orbit(2), ...
        rad2deg(orbit(3)), rad2deg(orbit(4)), rad2deg(orbit(5)), rad2deg(orbit(6)), ...
        'GravitationalParameter', muOrbit, 'Action', 'None');
    r_i = r_i_eci'; % [m]     1x3
    v_i = v_i_eci'; % [m/s]   1x3

    %% Compute initial quaternion
    % Reference body frame aligned with orbit:
    %   b1 = velocity direction, b3 = nadir (-r), b2 = completes triad
    b1 = v_i' / norm(v_i);
    b3 = -r_i' / norm(r_i);
    b2 = cross(b3, b1);
    b2 = b2 / norm(b2); % Ensure unit vector

    R_BI_ref = [b1'; b2'; b3']; % Reference DCM  (body <- inertial)

    % Apply user-specified Euler-angle offsets (ZYX: yaw-pitch-roll)
    roll  = deg2rad(simParams.initParams.Attitude.roll);
    pitch = deg2rad(simParams.initParams.Attitude.pitch);
    yaw   = deg2rad(simParams.initParams.Attitude.yaw);

    cr = cos(roll);  sr = sin(roll);
    cp = cos(pitch); sp = sin(pitch);
    cy = cos(yaw);   sy = sin(yaw);

    R_offset = [cy*cp,  cy*sp*sr - sy*cr,  cy*sp*cr + sy*sr;
                sy*cp,  sy*sp*sr + cy*cr,  sy*sp*cr - cy*sr;
                  -sp,           cp*sr,              cp*cr];

    R_BI_i = R_offset * R_BI_ref;
    beta_i = dcm2quat(R_BI_i);
    if beta_i(1) < 0
        beta_i = -beta_i;
    end

    % Initial angular velocity [rad/s]
    omega_i = deg2rad([simParams.initParams.Attitude.rollRate, ...
                       simParams.initParams.Attitude.pitchRate, ...
                       simParams.initParams.Attitude.yawRate]);

    %% Assemble initial state vector (13 x 1)
    X_i = [r_i, v_i, beta_i, omega_i]';

    %% Simulation time vector
    t0 = simParams.initParams.Simulation.initialTime;
    tf = simParams.initParams.Simulation.finalTime;
    dt = simParams.initParams.Simulation.timeStep;

    ts = t0 : dt : tf;

    %% ODE options
    opts = odeset('RelTol', simParams.initParams.Simulation.relTol, ...
                  'AbsTol', simParams.initParams.Simulation.absTol);

    %% Select dynamics model
    isNonlinear = strcmpi(simParams.initParams.Simulation.type, 'Nonlinear');
    if isNonlinear
        odeFun = @(t, X) sat_template_gui(t, X, subFigHandle);
    else
        Xr = [1; 0; 0; 0; 0; 0; 0];
        odeFun = @(t, X) Sat_template2_linear(t, X, Xr, X_i);
    end

    %% Run the integration with a progress dialog
    d = uiprogressdlg(subFigHandle, ...
        'Title', 'Running Simulation', ...
        'Message', sprintf('Integrating %s equations of motion...', ...
                           simParams.initParams.Simulation.type), ...
        'Indeterminate', 'on');

    try
        [t_out, X_out] = ode45(odeFun, ts, X_i, opts);

        %% Package results
        results        = struct();
        results.t      = t_out;
        results.X      = X_out;
        results.rs     = X_out(:, 1:3);
        results.vs     = X_out(:, 4:6);
        results.betas  = X_out(:, 7:10);
        results.omegas = X_out(:, 11:13);

        % Compute torque breakdown at each time step
        nSteps = length(t_out);
        torques = zeros(nSteps, 3);
        aeroTorques = zeros(nSteps, 3);
        totalTorques = zeros(nSteps, 3);

        for k = 1:nSteps
            tau = evaluateControlTorque(simParams, t_out(k), X_out(k, :)');
            if isNonlinear
                tauAero = evaluateAeroTorque(simParams, t_out(k), X_out(k, :)');
            else
                tauAero = [0; 0; 0];
            end

            torques(k, :) = tau';
            aeroTorques(k, :) = tauAero';
            totalTorques(k, :) = (tau + tauAero)';
        end

        results.torques = torques;
        results.aeroTorques = aeroTorques;
        results.totalTorques = totalTorques;
        results.maxAeroTorqueNorm = max(sqrt(sum(aeroTorques.^2, 2)));

        csvPath = saveTorqueBreakdownCsv(simParams, t_out, torques, aeroTorques, totalTorques);
        results.torqueCsvPath = csvPath;

        % Store results in guidata
        simParams.results = results;
        guidata(subFigHandle, simParams);

        close(d);
        if isempty(csvPath)
            successMsg = sprintf('Simulation completed successfully!\n%d time steps integrated.\nMax |tau_aero| = %.3e N m\nTorque breakdown CSV could not be saved.', nSteps, results.maxAeroTorqueNorm);
        else
            [~, csvName, csvExt] = fileparts(csvPath);
            successMsg = sprintf('Simulation completed successfully!\n%d time steps integrated.\nMax |tau_aero| = %.3e N m\nTorque breakdown saved to simulations/%s%s', nSteps, results.maxAeroTorqueNorm, csvName, csvExt);
        end

        uialert(subFigHandle, successMsg, 'Success', 'Icon', 'success');

    catch ME
        close(d);
        uialert(subFigHandle, ...
            sprintf('Simulation failed:\n%s', ME.message), ...
            'Error', 'Icon', 'error');
    end
end

function simParams = ensureSimulationDefaults(simParams)
    if ~isfield(simParams, 'initParams')
        simParams.initParams = struct();
    end

    if ~isfield(simParams.initParams, 'Modes')
        simParams.initParams.Modes = struct();
    end

    if ~isfield(simParams.initParams.Modes, 'enableAero')
        simParams.initParams.Modes.enableAero = true;
    end

    if ~isfield(simParams.initParams.Modes, 'enableControl')
        simParams.initParams.Modes.enableControl = true;
    end

    if ~isfield(simParams.initParams, 'Environment')
        simParams.initParams.Environment = struct();
    end

    if ~isfield(simParams.initParams.Environment, 'gasSurfaceInteractionModel') && ...
            isfield(simParams.initParams.Environment, 'gasSurfaceModel')
        simParams.initParams.Environment.gasSurfaceInteractionModel = ...
            simParams.initParams.Environment.gasSurfaceModel;
    end

    if ~isfield(simParams.initParams.Environment, 'gasSurfaceInteractionModel')
        simParams.initParams.Environment.gasSurfaceInteractionModel = 'cook';
    end

    simParams.initParams.Environment.gasSurfaceInteractionModel = ...
        normalizeGsiModelValue(simParams.initParams.Environment.gasSurfaceInteractionModel);

    if ~isfield(simParams.initParams.Environment, 'year')
        simParams.initParams.Environment.year = 2002;
    end

    if ~isfield(simParams.initParams, 'Controller')
        simParams.initParams.Controller = struct();
    end

    if ~isfield(simParams.initParams.Controller, 'functionFile') || ...
            isempty(simParams.initParams.Controller.functionFile)
        simParams.initParams.Controller.functionFile = 'control_torques.m';
    end

    if ~isfield(simParams.initParams.Controller, 'Func') || ...
            ~isa(simParams.initParams.Controller.Func, 'function_handle')
        [~, funcName, ~] = fileparts(simParams.initParams.Controller.functionFile);
        simParams.initParams.Controller.Func = str2func(funcName);
    end
end

function tau = evaluateControlTorque(simParams, t, X)
    tau = [0; 0; 0];

    if ~getModeFlag(simParams, 'enableControl', true)
        return;
    end

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Controller')
        return;
    end

    controller = simParams.initParams.Controller;
    controllerFunc = [];

    if isfield(controller, 'Func') && isa(controller.Func, 'function_handle')
        controllerFunc = controller.Func;
    elseif isfield(controller, 'functionFile') && ~isempty(controller.functionFile)
        [~, funcName, ~] = fileparts(controller.functionFile);
        controllerFunc = str2func(funcName);
    end

    if isempty(controllerFunc)
        return;
    end

    try
        tauCandidate = controllerFunc(t, X);
        tauCandidate = reshape(tauCandidate, [], 1);
        if numel(tauCandidate) == 3 && all(isfinite(tauCandidate))
            tau = tauCandidate;
        end
    catch
        tau = [0; 0; 0];
    end
end

function tauAero = evaluateAeroTorque(simParams, t, X)
    tauAero = [0; 0; 0];

    if ~getModeFlag(simParams, 'enableAero', true)
        return;
    end

    r = norm(X(1:3));
    if ~(isfinite(r) && r > 0)
        return;
    end

    thisFilePath = fileparts(mfilename('fullpath'));
    objFilePath = fullfile(thisFilePath, '..', 'dynamics', '6U CubeSat.obj');

    location.positionECI = X(1:3);

    q = X(7:10);
    v_eci = X(4:6);

    R_eci2body = quat2dcm(q');

    v_body = R_eci2body * v_eci;
    v_mag = norm(v_body);
    if v_mag > 1
        attitude.aoa = atan2d(v_body(3), v_body(1));
        attitude.aos = atan2d(v_body(2), v_body(1));
    else
        attitude.aoa = 0;
        attitude.aos = 0;
    end

    day_start = getEnvValue(simParams, 'dayOfYear', 1);
    seconds_start = getEnvValue(simParams, 'secondsOfDay', 0);
    time.year = getEnvValue(simParams, 'year', 2002);
    time.dayOfYear = day_start + floor((seconds_start + t) / 86400);
    time.UTseconds = mod(seconds_start + t, 86400);

    aeroOptions.f107Average = getEnvValue(simParams, 'F107Average', 150);
    aeroOptions.f107Daily = getEnvValue(simParams, 'F107Daily', 150);
    aeroOptions.magneticIndex = getEnvValue(simParams, 'magneticIndices', ones(1, 7) * 3);
    aeroOptions.anomalousOxygen = getEnvValue(simParams, 'enableAnomalousOxygen', false);
    aeroOptions.gsi_model = normalizeGsiModelValue(getEnvValue(simParams, 'gasSurfaceInteractionModel', 'cook'));
    aeroOptions.alpha = getEnvValue(simParams, 'accommodationCoefficient', 1);
    aeroOptions.Tw = getEnvValue(simParams, 'wallTemperature', 300);
    aeroOptions.enableShadow = getEnvValue(simParams, 'enableShadowAnalysis', true);
    aeroOptions.enableSolar = getEnvValue(simParams, 'enableSolarRadiationPressure', true);
    aeroOptions.sol_cR = getEnvValue(simParams, 'specularReflectivity', 0.15);
    aeroOptions.sol_cD = getEnvValue(simParams, 'diffuseReflectivity', 0.25);

    try
        aeroResults = computeAeroForces(objFilePath, location, attitude, time, aeroOptions);
        tauCandidate = reshape(aeroResults.M_aero, [], 1);

        if isfield(aeroResults, 'M_solar')
            tauCandidate = tauCandidate + reshape(aeroResults.M_solar, [], 1);
        end

        if numel(tauCandidate) == 3 && all(isfinite(tauCandidate))
            tauAero = tauCandidate;
        end
    catch
        tauAero = [0; 0; 0];
    end
end

function csvPath = saveTorqueBreakdownCsv(simParams, tOut, controlTorques, aeroTorques, totalTorques)
    csvPath = '';

    thisFilePath = fileparts(mfilename('fullpath'));
    simDirectory = fullfile(thisFilePath, '..', '..', 'simulations');
    if ~exist(simDirectory, 'dir')
        mkdir(simDirectory);
    end

    simName = 'simulation';
    if isfield(simParams, 'initParams') && isfield(simParams.initParams, 'Simulation') && ...
            isfield(simParams.initParams.Simulation, 'name') && ~isempty(simParams.initParams.Simulation.name)
        simName = simParams.initParams.Simulation.name;
    end

    simName = regexprep(simName, '[^a-zA-Z0-9_\- ]', '_');
    simName = strtrim(simName);
    simName = regexprep(simName, '\s+', '_');
    if isempty(simName)
        simName = 'simulation';
    end

    csvPath = fullfile(simDirectory, [simName, '_torque_breakdown.csv']);

    torqueTable = table( ...
        tOut, ...
        controlTorques(:, 1), controlTorques(:, 2), controlTorques(:, 3), ...
        aeroTorques(:, 1), aeroTorques(:, 2), aeroTorques(:, 3), ...
        totalTorques(:, 1), totalTorques(:, 2), totalTorques(:, 3), ...
        'VariableNames', { ...
            'time_s', ...
            'tau_control_x', 'tau_control_y', 'tau_control_z', ...
            'tau_aero_x', 'tau_aero_y', 'tau_aero_z', ...
            'tau_total_x', 'tau_total_y', 'tau_total_z' ...
        } ...
    );

    try
        writetable(torqueTable, csvPath);
    catch
        csvPath = '';
    end
end

function value = getEnvValue(simParams, fieldName, defaultValue)
    value = defaultValue;

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Environment')
        return;
    end

    if isfield(simParams.initParams.Environment, fieldName)
        value = simParams.initParams.Environment.(fieldName);
    end
end

function enabled = getModeFlag(simParams, fieldName, defaultValue)
    enabled = defaultValue;

    if ~isfield(simParams, 'initParams') || ~isfield(simParams.initParams, 'Modes')
        return;
    end

    if isfield(simParams.initParams.Modes, fieldName)
        enabled = logical(simParams.initParams.Modes.(fieldName));
    end
end

function modelName = normalizeGsiModelValue(modelName)
    modelName = lower(char(modelName));
    if ~any(strcmp(modelName, {'cook', 'sentman'}))
        modelName = 'cook';
    end
end
