% optimize_thruster_orientation.m
% Finds the optimal thruster canting angle to maximize the minimum angular 
% acceleration (equalizing agility across all 3 axes) for a 6U CubeSat.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc;

%% 6U CubeSat Parameters
dims_m  = [0.10, 0.20, 0.30];       
I_body  = diag([0.13, 0.10, 0.05]); 

a = dims_m(1) / 2; % 0.05
b = dims_m(2) / 2; % 0.10
c = dims_m(3) / 2; % 0.15

Ix = I_body(1,1);
Iy = I_body(2,2);
Iz = I_body(3,3);

% We search over the first octant of a unit sphere
% Exhaust direction = [dx, dy, dz]
% Reaction force = [-dx, -dy, -dz]
% Torque from T1 at [a, b, c] is cross([a,b,c], [-dx,-dy,-dz])
% tau_x = abs(c*dy - b*dz)
% tau_y = abs(a*dz - c*dx)
% tau_z = abs(b*dx - a*dy)

N = 500;
azimuth = linspace(0, pi/2, N);
elevation = linspace(0, pi/2, N);

[AZ, EL] = meshgrid(azimuth, elevation);

DX = cos(EL) .* cos(AZ);
DY = cos(EL) .* sin(AZ);
DZ = sin(EL);

% Base torque arrays (for 1 thruster, we multiply by 4 later for the firing group)
TAU_X = abs(c.*DY - b.*DZ);
TAU_Y = abs(a.*DZ - c.*DX);
TAU_Z = abs(b.*DX - a.*DY);

ALPHA_X = (4 * TAU_X) / Ix;
ALPHA_Y = (4 * TAU_Y) / Iy;
ALPHA_Z = (4 * TAU_Z) / Iz;

% We want to maximize the MINIMUM angular acceleration
MIN_ALPHA = min(min(ALPHA_X, ALPHA_Y), ALPHA_Z);

% Constraint: We must maintain full 6-DOF translation. 
% Standard [1,1,1] gives 2.309 N per N of thrust.
% Let's enforce that no translation axis drops below 1.5 N per N of thrust (meaning dx, dy, dz >= 1.5/4 = 0.375)
MIN_TRANS_FORCE = 1.5;
valid_mask = (4*DX >= MIN_TRANS_FORCE) & (4*DY >= MIN_TRANS_FORCE) & (4*DZ >= MIN_TRANS_FORCE);

% Set invalid configurations to a terrible score
MIN_ALPHA(~valid_mask) = -1;

[max_min_alpha, max_idx] = max(MIN_ALPHA(:));

opt_dx = DX(max_idx);
opt_dy = DY(max_idx);
opt_dz = DZ(max_idx);

opt_ax = ALPHA_X(max_idx);
opt_ay = ALPHA_Y(max_idx);
opt_az = ALPHA_Z(max_idx);

fprintf('=== OPTIMIZATION RESULTS ===\n');
fprintf('Objective: Maximize minimum angular acceleration (Minimax Agility)\n\n');

fprintf('Standard [1,1,1] diagonal configuration:\n');
dx_std = 1/sqrt(3); dy_std = 1/sqrt(3); dz_std = 1/sqrt(3);
ax_std = 4*abs(c*dy_std - b*dz_std)/Ix;
ay_std = 4*abs(a*dz_std - c*dx_std)/Iy;
az_std = 4*abs(b*dx_std - a*dy_std)/Iz;
fprintf('  Exhaust Dir: [%.4f, %.4f, %.4f]\n', dx_std, dy_std, dz_std);
fprintf('  Alpha [rad/s^2 per N]: X = %.4f, Y = %.4f, Z = %.4f\n', ax_std, ay_std, az_std);
fprintf('  Min Alpha: %.4f\n\n', min([ax_std, ay_std, az_std]));

fprintf('Optimal Asymmetric configuration:\n');
fprintf('  Exhaust Dir: [%.4f, %.4f, %.4f]\n', opt_dx, opt_dy, opt_dz);
% Normalize to make the max component 1.0 for easier comparison with [1,1,1]
norm_dir = [opt_dx, opt_dy, opt_dz] / max([opt_dx, opt_dy, opt_dz]);
fprintf('  Normalized Dir: [%.4f, %.4f, %.4f]\n', norm_dir(1), norm_dir(2), norm_dir(3));
fprintf('  Alpha [rad/s^2 per N]: X = %.4f, Y = %.4f, Z = %.4f\n', opt_ax, opt_ay, opt_az);
fprintf('  Min Alpha: %.4f (%.1f%% improvement in bottleneck axis!)\n\n', max_min_alpha, ((max_min_alpha/min([ax_std, ay_std, az_std])) - 1)*100);

fprintf('Translational efficiency check (Force per N per axis):\n');
fprintf('  Standard Tx, Ty, Tz: %.4f, %.4f, %.4f\n', 4*dx_std, 4*dy_std, 4*dz_std);
fprintf('  Optimal  Tx, Ty, Tz: %.4f, %.4f, %.4f\n', 4*opt_dx, 4*opt_dy, 4*opt_dz);

