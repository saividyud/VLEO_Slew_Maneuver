clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

fprintf('=== Testing sat_dynamics_nonlinear with aerodynamics ===\n\n');

c = vleo.util.constants();
earthRadius = c.R_earth;
mu = c.mu_earth;
altitude = 200e3;
r0 = earthRadius + altitude;
v0 = sqrt(mu / r0);

X0 = zeros(13, 1);
X0(1) = r0;
X0(5) = v0;
X0(7) = 1;

fprintf('Initial conditions:\n');
fprintf('  Altitude: %.0f km\n', altitude / 1e3);
fprintf('  Orbital velocity: %.2f m/s\n', v0);
fprintf('  Position: [%.2f, %.2f, %.2f] km\n', X0(1:3) / 1e3);
fprintf('  Velocity: [%.2f, %.2f, %.2f] m/s\n\n', X0(4:6));

fprintf('Testing single evaluation of sat_dynamics_nonlinear...\n');
tic;
try
    Xd = vleo.dynamics.sat_dynamics_nonlinear(0, X0);
    elapsed = toc;
    fprintf('  SUCCESS! Completed in %.4f seconds\n\n', elapsed);
    fprintf('State derivatives:\n');
    fprintf('  Velocity (Xd 1:3):     [%.6e, %.6e, %.6e] m/s\n', Xd(1:3));
    fprintf('  Acceleration (Xd 4:6): [%.6e, %.6e, %.6e] m/s^2\n', Xd(4:6));
    fprintf('  Quat rate (Xd 7:10):   [%.6e, %.6e, %.6e, %.6e]\n', Xd(7:10));
    fprintf('  Ang accel (Xd 11:13):  [%.6e, %.6e, %.6e] rad/s^2\n\n', Xd(11:13));
    aGravity = -mu / r0^2;
    fprintf('Expected gravity acceleration: %.6e m/s^2\n', aGravity);
    fprintf('Computed radial acceleration:  %.6e m/s^2\n', Xd(4));
    fprintf('Difference (includes J2+aero): %.6e m/s^2\n', Xd(4) - aGravity);
catch ME
    elapsed = toc;
    fprintf('  FAILED after %.4f seconds\n', elapsed);
    fprintf('  Error: %s\n', ME.message);
    for idx = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(idx).name, ME.stack(idx).line);
    end
end

fprintf('\n=== Running short integration (10 seconds) ===\n');
try
    tic;
    [tOut, XOut] = ode45(@vleo.dynamics.sat_dynamics_nonlinear, [0, 10], X0);
    fprintf('  SUCCESS! Integrated %d time steps in %.4f seconds\n', length(tOut), toc);
    fprintf('  Final altitude: %.2f km\n', (norm(XOut(end, 1:3)) - earthRadius) / 1e3);
catch ME
    fprintf('  FAILED: %s\n', ME.message);
end
