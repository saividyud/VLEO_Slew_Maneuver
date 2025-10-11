% =========================================================================
% plot_results.m
%
% Description:
% This function visualizes the results of the orbit propagation comparison.
% It generates a multi-panel figure showing the 3D orbit trajectories,
% altitude profiles, position errors relative to SGP4, and computation
% performance. It also prints a summary of the key metrics to the console.
%
% Inputs:
%   tsince_vec (array): Vector of time steps since epoch [minutes].
%   r_sgp4     (3xN): Position vectors from SGP4 [m].
%   r_custom   (3xN): Position vectors from the custom propagator [m].
%   r_hpop     (3xN): Position vectors from HPOP [m].
%   t_sgp4     (scalar): Computation time for SGP4 [s].
%   t_custom   (scalar): Computation time for the custom propagator [s].
%   t_hpop     (scalar): Computation time for HPOP [s].
% =========================================================================

function plot_results(tsince_vec, r_sgp4, r_custom, r_hpop, t_sgp4, t_custom, t_hpop)
    
    hpop_success = ~any(isnan(r_hpop(:)));
    n_steps = length(tsince_vec);
    propagation_hours = tsince_vec(end)/60;
    R_e = 6378137; % Earth radius in meters

    %% Calculate Errors and Altitudes
    err_custom = vecnorm(r_custom - r_sgp4);
    alt_sgp4   = (vecnorm(r_sgp4) - R_e) / 1e3;
    alt_custom = (vecnorm(r_custom) - R_e) / 1e3;
    
    if hpop_success
        err_hpop           = vecnorm(r_hpop - r_sgp4);
        err_custom_vs_hpop = vecnorm(r_custom - r_hpop);
        alt_hpop           = (vecnorm(r_hpop) - R_e) / 1e3;
    end
    
    %% Print Results Summary
    fprintf('========================================\n');
    fprintf(' RESULTS\n');
    fprintf('========================================\n\n');

    fprintf('Computation Time:\n');
    fprintf('  SGP4:            %.4f s (baseline)\n', t_sgp4);
    fprintf('  Custom Numerical: %.4f s (%.1fx)\n', t_custom, t_custom/t_sgp4);
    if hpop_success
        fprintf('  HPOP (20x20):     %.4f s (%.1fx)\n\n', t_hpop, t_hpop/t_sgp4);
    end

    fprintf('Position Error vs SGP4:\n');
    fprintf('  Custom Numerical:\n');
    fprintf('    Mean:  %.3f km\n', mean(err_custom)/1e3);
    fprintf('    Max:   %.3f km\n', max(err_custom)/1e3);
    fprintf('    Final: %.3f km\n', err_custom(end)/1e3);

    if hpop_success
        fprintf('  HPOP:\n');
        fprintf('    Mean:  %.3f km\n', mean(err_hpop)/1e3);
        fprintf('    Max:   %.3f km\n', max(err_hpop)/1e3);
        fprintf('    Final: %.3f km\n\n', err_hpop(end)/1e3);

        fprintf('Position Difference (Custom vs HPOP "Truth"):\n');
        fprintf('    Mean:  %.3f km\n', mean(err_custom_vs_hpop)/1e3);
        fprintf('    Max:   %.3f km\n', max(err_custom_vs_hpop)/1e3);
    end
    
    fprintf('\nAltitude Decay Over %.0f Hours:\n', propagation_hours);
    decay_sgp4 = alt_sgp4(1) - alt_sgp4(end);
    decay_custom = alt_custom(1) - alt_custom(end);
    fprintf('  SGP4:            %.3f km\n', decay_sgp4);
    fprintf('  Custom Numerical: %.3f km\n', decay_custom);
    if hpop_success
        decay_hpop = alt_hpop(1) - alt_hpop(end);
        fprintf('  HPOP:            %.3f km\n', decay_hpop);
    end

    %% Create Plots
    time_hours = tsince_vec / 60;
    figure('Name', 'Orbit Propagator Comparison', 'Position', [50, 50, 1600, 800]);

    % 3D Orbit Plot
    subplot(2, 3, 1);
    hold on; grid on; axis equal;
    plot3(r_sgp4(1,:)/1e3, r_sgp4(2,:)/1e3, r_sgp4(3,:)/1e3, 'b', 'LineWidth', 1.5, 'DisplayName', 'SGP4');
    plot3(r_custom(1,:)/1e3, r_custom(2,:)/1e3, r_custom(3,:)/1e3, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom');
    if hpop_success
        plot3(r_hpop(1,:)/1e3, r_hpop(2,:)/1e3, r_hpop(3,:)/1e3, 'g:', 'LineWidth', 2, 'DisplayName', 'HPOP');
    end
    [X,Y,Z] = sphere(50);
    surf(X*R_e/1e3, Y*R_e/1e3, Z*R_e/1e3, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    title('3D Orbit Trajectories'); xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    legend; view(3);

    % Altitude vs Time
    subplot(2, 3, 2);
    hold on; grid on;
    plot(time_hours, alt_sgp4, 'b', 'LineWidth', 1.5, 'DisplayName', 'SGP4');
    plot(time_hours, alt_custom, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom');
    if hpop_success
        plot(time_hours, alt_hpop, 'g:', 'LineWidth', 2, 'DisplayName', 'HPOP');
    end
    title('Altitude vs. Time'); xlabel('Time (hours)'); ylabel('Altitude (km)');
    legend('Location', 'southwest');
    
    % Position Error vs SGP4
    subplot(2, 3, 3);
    hold on; grid on;
    plot(time_hours, err_custom/1e3, 'r', 'LineWidth', 1.5, 'DisplayName', 'Custom vs SGP4');
    if hpop_success
        plot(time_hours, err_hpop/1e3, 'g', 'LineWidth', 1.5, 'DisplayName', 'HPOP vs SGP4');
    end
    title('Position Error relative to SGP4'); xlabel('Time (hours)'); ylabel('Error (km)');
    legend;
    
    % Altitude Decay
    subplot(2, 3, 4);
    hold on; grid on;
    plot(time_hours, (alt_sgp4 - alt_sgp4(1)), 'b', 'LineWidth', 1.5, 'DisplayName', 'SGP4');
    plot(time_hours, (alt_custom - alt_custom(1)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Custom');
    if hpop_success
        plot(time_hours, (alt_hpop - alt_hpop(1)), 'g:', 'LineWidth', 2, 'DisplayName', 'HPOP');
    end
    title('Altitude Decay from Initial'); xlabel('Time (hours)'); ylabel('Altitude Change (km)');
    legend('Location', 'southwest');

    % Difference between Custom and HPOP
    subplot(2, 3, 5);
    if hpop_success
        plot(time_hours, err_custom_vs_hpop/1e3, 'k', 'LineWidth', 1.5);
        grid on;
        title('Custom vs. HPOP ("Truth")'); xlabel('Time (hours)'); ylabel('Difference (km)');
    else
        axis off;
        text(0.5, 0.5, 'HPOP results not available', 'HorizontalAlignment', 'center');
    end

    % Computation Time
    subplot(2, 3, 6);
    times = [t_sgp4, t_custom];
    labels = {'SGP4', 'Custom'};
    if hpop_success
        times(3) = t_hpop;
        labels{3} = 'HPOP';
    end
    bar(times);
    set(gca, 'xticklabel', labels);
    grid on;
    title('Computation Time'); ylabel('Time (s)');

end
