function Xd = sat_dynamics_gui(t, X, subFigHandle)
    simParams = guidata(subFigHandle);
    c = vleo.util.constants();
    mu = vleo.util.get_env_value(simParams, 'mu', c.mu_earth);

    aeroConfig = [];
    if vleo.util.get_mode_flag(simParams, 'enableAero', true)
        aeroConfig = vleo.dynamics.build_aero_config_from_sim_params(simParams);
    end

    tauControl = vleo.dynamics.evaluate_controller_torque(simParams, t, X);
    Xd = vleo.dynamics.sat_dynamics_base(t, X, mu, aeroConfig, tauControl);
end
