function runSimulation(subFigHandle)
    salmonColor = [1, 0.4941, 0.4392];

    % Reading in simParams from guidata to pass into this function so we can update it with new values
    simParams = guidata(subFigHandle);

    % Extracting orbital parameters
    a = simParams.initParams.Orbit.semiMajorAxis;
    e = simParams.initParams.Orbit.eccentricity;
    i = simParams.initParams.Orbit.inclination;
    argPeriapse = simParams.initParams.Orbit.argPeriapse;
    RAAN = simParams.initParams.Orbit.RAAN;
    trueAnomaly = simParams.initParams.Orbit.trueAnomaly;

    orbit = [a, e, i, argPeriapse, RAAN, trueAnomaly];

    % Converting to state
    RV = RVfromOE(orbit, simParams.initParams.Environment.mu);
    r_i = RV(:, 1)'; % [km]
    v_i = RV(:, 2)'; % [km/s]

    % Extracting attitude parameters
    beta_i = simParams.initParams.Attitude.Quaternion;
    omega_i = [simParams.initParams.Attitude.rollRate, ...
               simParams.initParams.Attitude.pitchRate, ...
               simParams.initParams.Attitude.yawRate];

    X_i = [r_i, v_i, beta_i, omega_i]';

    % Extracting simulation parameters
    t_i = simParams.initParams.Simulation.initialTime;
    t_f = simParams.initParams.Simulation.finalTime;
    dt = simParams.initParams.Simulation.timeStep;

    % simParams.finalParams will be populated after simulation runs
    simParams.finalParams = struct();

    ts = t_i : dt : t_f;

    % Running the simulation
    opts = odeset('RelTol', simParams.initParams.Simulation.relTol, 'AbsTol', simParams.initParams.Simulation.absTol);
    [t, X] = ode45(@(t, X) sat_template_gui(t, X, subFigHandle), ts, X_i, opts);

    simParams.finalParams.X = X;
    simParams.finalParams.t = t;

    % Extracting torques from controller
    torques = zeros(length(t), 3);
    for idx = 1:length(t)
        torques(idx, :) = simParams.initParams.Controller.Func(t(idx), X(idx, :)');
    end

    simParams.finalParams.rs = X(:, 1:3);
    simParams.finalParams.vs = X(:, 4:6);
    simParams.finalParams.betas = X(:, 7:10);
    simParams.finalParams.omegas = X(:, 11:13);
    simParams.finalParams.ControlTorques = torques;

    guidata(subFigHandle, simParams); % Updating simParams in guidata to pass into sat_template_gui
    assignin("base", "simParams", simParams); % Assigning simParams to base workspace for plotting

end