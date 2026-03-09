% test_Sat_template.m
% Simple test to verify Sat_template runs with aerodynamic forces
%
% Tests a 6U CubeSat at 200 km altitude for a short integration

run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));

clear; clc;

fprintf('=== Testing Sat_template with Aerodynamics ===\n\n');

%% Initial conditions

% Orbital parameters (200 km circular orbit)
R_E = 6.37813649e6;     % Earth radius [m]
mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
alt = 200e3;            % Altitude [m]
r0 = R_E + alt;         % Orbital radius [m]
v0 = sqrt(mu/r0);       % Circular velocity [m/s]

% Initial state vector (13 elements)
% Position: on x-axis
% Velocity: in y-direction (circular orbit in xy-plane)
% Quaternion: identity (no rotation)
% Angular velocity: zero

X0 = zeros(13, 1);
X0(1) = r0;             % x position [m]
X0(5) = v0;             % y velocity [m/s]
X0(7) = 1;              % q0 (scalar part of quaternion)
% X0(8:10) = 0          % q1, q2, q3 (vector part)
% X0(11:13) = 0         % angular velocity [rad/s]

fprintf('Initial conditions:\n');
fprintf('  Altitude: %.0f km\n', alt/1e3);
fprintf('  Orbital velocity: %.2f m/s\n', v0);
fprintf('  Position: [%.2f, %.2f, %.2f] km\n', X0(1:3)/1e3);
fprintf('  Velocity: [%.2f, %.2f, %.2f] m/s\n', X0(4:6));
fprintf('\n');

%% Test single evaluation of Sat_template

fprintf('Testing single evaluation of Sat_template...\n');
tic;
try
    Xd = Sat_template(0, X0);
    elapsed = toc;
    fprintf('  SUCCESS! Completed in %.4f seconds\n\n', elapsed);
    
    fprintf('State derivatives:\n');
    fprintf('  Velocity (Xd 1:3):     [%.6e, %.6e, %.6e] m/s\n', Xd(1:3));
    fprintf('  Acceleration (Xd 4:6): [%.6e, %.6e, %.6e] m/s^2\n', Xd(4:6));
    fprintf('  Quat rate (Xd 7:10):   [%.6e, %.6e, %.6e, %.6e]\n', Xd(7:10));
    fprintf('  Ang accel (Xd 11:13):  [%.6e, %.6e, %.6e] rad/s^2\n', Xd(11:13));
    fprintf('\n');
    
    % Check expected values
    a_gravity = -mu/r0^2;  % Expected gravitational acceleration
    fprintf('Expected gravity acceleration: %.6e m/s^2\n', a_gravity);
    fprintf('Computed radial acceleration:  %.6e m/s^2\n', Xd(4));
    fprintf('Difference (includes J2+aero): %.6e m/s^2\n', Xd(4) - a_gravity);
    
catch ME
    elapsed = toc;
    fprintf('  FAILED after %.4f seconds\n', elapsed);
    fprintf('  Error: %s\n', ME.message);
    fprintf('  Stack:\n');
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

%% Short integration test

fprintf('\n=== Running short integration (10 seconds) ===\n');
tspan = [0, 10];  % 10 seconds

try
    tic;
    [t, X] = ode45(@Sat_template, tspan, X0);
    elapsed = toc;
    
    fprintf('Integration completed in %.4f seconds\n', elapsed);
    fprintf('Number of steps: %d\n', length(t));
    fprintf('\nFinal state:\n');
    fprintf('  Position: [%.2f, %.2f, %.2f] km\n', X(end,1:3)/1e3);
    fprintf('  Velocity: [%.2f, %.2f, %.2f] m/s\n', X(end,4:6));
    fprintf('  Quaternion: [%.6f, %.6f, %.6f, %.6f]\n', X(end,7:10));
    fprintf('  Angular vel: [%.6e, %.6e, %.6e] rad/s\n', X(end,11:13));
    
    % Check orbit didn't explode
    r_final = norm(X(end,1:3));
    alt_final = r_final - R_E;
    fprintf('\nFinal altitude: %.2f km (change: %.2f m)\n', alt_final/1e3, alt_final - alt);
    
    fprintf('\n=== TEST PASSED ===\n');
    
catch ME
    fprintf('Integration FAILED\n');
    fprintf('Error: %s\n', ME.message);
end
