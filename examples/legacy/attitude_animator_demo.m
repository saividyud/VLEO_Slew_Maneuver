projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(projectRoot);
setup_project();

simParams = vleo.gui.default_simulation_params();
simParams.initParams.Simulation.finalTime = 5 * 60;
simParams.initParams.Orbit.altitude = 250e3;
simParams.initParams.Orbit.semiMajorAxis = 6378.14e3 + 250e3;

X0 = vleo.dynamics.initial_state_from_sim_params(simParams);
ts = 0:1:simParams.initParams.Simulation.finalTime;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[tOut, XOut] = ode45(@vleo.dynamics.sat_dynamics_nonlinear, ts, X0, opts);

results = struct();
results.t = tOut;
results.rs = XOut(:, 1:3);
results.vs = XOut(:, 4:6);
results.betas = XOut(:, 7:10);
results.omegas = XOut(:, 11:13);
results.torques = zeros(length(tOut), 3);
for idx = 1:length(tOut)
    results.torques(idx, :) = vleo.control.attitude_pd_controller(tOut(idx), XOut(idx, :)')';
end

simParams.results = results;
vleo.viz.animate_results(simParams);
