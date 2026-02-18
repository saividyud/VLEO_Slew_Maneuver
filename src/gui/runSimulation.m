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

    %% Extract simulation parameters
    simParams = guidata(subFigHandle);

    % Constants
    earthRadius = 6378.14e3; % Earth equatorial radius [m]

    %% Build initial orbital state
    % simParams stores altitude and semiMajorAxis in m
    a   = simParams.initParams.Orbit.semiMajorAxis; % [m]
    e   = simParams.initParams.Orbit.eccentricity;
    inc = deg2rad(simParams.initParams.Orbit.inclination);
    raan = deg2rad(simParams.initParams.Orbit.RAAN);
    aop  = deg2rad(simParams.initParams.Orbit.argPeriapse);
    ta   = deg2rad(simParams.initParams.Orbit.trueAnomaly);

    orbit = [a, e, inc, raan, aop, ta];

    % Compute initial position and velocity in ECI frame
    RV  = RVfromOE(orbit);
    r_i = RV(:, 1)'; % [m]     1x3
    v_i = RV(:, 2)'; % [m/s]   1x3

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
    beta_i = QfromDCM(R_BI_i)'; % 1x4 quaternion [q0 q1 q2 q3]

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
    if strcmpi(simParams.initParams.Simulation.type, 'Nonlinear')
        odeFun = @Sat_template;
    else
        odeFun = @Sat_template2_linear;
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

        % Compute control torques at each time step
        nSteps = length(t_out);
        torques = zeros(nSteps, 3);
        for k = 1:nSteps
            torques(k, :) = control_torques(t_out(k), X_out(k, :)')';
        end
        results.torques = torques;

        % Store results in guidata
        simParams.results = results;
        guidata(subFigHandle, simParams);

        close(d);
        uialert(subFigHandle, ...
            sprintf('Simulation completed successfully!\n%d time steps integrated.', nSteps), ...
            'Success', 'Icon', 'success');

    catch ME
        close(d);
        uialert(subFigHandle, ...
            sprintf('Simulation failed:\n%s', ME.message), ...
            'Error', 'Icon', 'error');
    end
end
