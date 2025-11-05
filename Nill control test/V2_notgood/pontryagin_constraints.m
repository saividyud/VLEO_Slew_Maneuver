function [c, ceq] = pontryagin_constraints(x, X0, r_volcano_0, params, ...
    ra_target, dec_target, ra_dot_target, dec_dot_target, n_intervals)
    
    t_f = x(1);
    tau_all = reshape(x(2:end), n_intervals, 3);
    t_intervals = linspace(0, t_f, n_intervals+1);
    X_current = X0;
    
    for i = 1:n_intervals
        t_start = t_intervals(i);
        t_end = t_intervals(i+1);
        tau_i = tau_all(i, :)';
        
        % Fixed: Use Sat_template instead of Sat_template_optimal
        [~, X_seg] = ode45(@(t,X) Sat_template(t, X, tau_i, params), ...
            [t_start, t_end], X_current, ...
            odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
        X_current = X_seg(end, :)';
    end
    
    % Terminal constraints
    obs_final = state_to_observation(X_current, params);
    ra_final = deg2rad(obs_final.ra);
    dec_final = deg2rad(obs_final.dec);
    
    dt_check = 0.5;
    [~, X_plus] = ode45(@(t,X) Sat_template(t, X, zeros(3,1), params), ...
        [t_f, t_f+dt_check], X_current, ...
        odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    X_plus_state = X_plus(end, :)';
    obs_plus = state_to_observation(X_plus_state, params);
    
    ra_dot_final = (deg2rad(obs_plus.ra) - ra_final) / dt_check;
    dec_dot_final = (deg2rad(obs_plus.dec) - dec_final) / dt_check;
    
    ceq = [ra_final - ra_target; dec_final - dec_target; ...
           ra_dot_final - ra_dot_target; dec_dot_final - dec_dot_target];
    c = [];
end
