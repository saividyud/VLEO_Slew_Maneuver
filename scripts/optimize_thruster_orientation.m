% optimize_thruster_orientation.m
% Compare the standard and optimized 8-corner thruster cant directions.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc;

params = vleo.dynamics.default_control_test_params(true);
standardLayout = vleo.control.corner_thruster_layout(params, 'standard');
optimizedLayout = vleo.control.corner_thruster_layout(params, 'optimized');

fprintf('=== OPTIMIZATION RESULTS ===\n');
fprintf('Objective: Maximize minimum angular acceleration (minimax agility)\n\n');

print_layout_summary('Standard [1,1,1] diagonal configuration', standardLayout);
print_layout_summary('Optimal asymmetric configuration', optimizedLayout);

improvementPct = 100 * (optimizedLayout.minimax_score / standardLayout.minimax_score - 1);
normalizedDirection = optimizedLayout.exhaust_direction_unit / max(optimizedLayout.exhaust_direction_unit);

fprintf('  Normalized optimized direction: [%.4f, %.4f, %.4f]\n', normalizedDirection);
fprintf('  Bottleneck-axis improvement: %.1f%%\n', improvementPct);

function print_layout_summary(titleText, layout)
    fprintf('%s:\n', titleText);
    fprintf('  Exhaust dir: [%.4f, %.4f, %.4f]\n', layout.exhaust_direction_unit);
    fprintf('  Translation [N/N]: X = %.4f, Y = %.4f, Z = %.4f\n', ...
        layout.translation_net_per_unit_thrust);
    fprintf('  Alpha [rad/s^2 per N]: X = %.4f, Y = %.4f, Z = %.4f\n', ...
        layout.angular_accel_per_unit_thrust);
    fprintf('  Min alpha: %.4f\n\n', layout.minimax_score);
end
