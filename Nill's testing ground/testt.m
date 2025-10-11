% =========================================================================
% main.m
%
% Description:
% This is the main script to run a three-way orbit propagator comparison
% between SGP4, a custom numerical propagator, and the High-Precision
% Orbit Propagator (HPOP). It serves as the primary driver for the
% simulation, initializing parameters, calling propagation routines,
% and generating final plots and analysis.
%
% For your 6U CubeSat project, you can adapt this script to test different
% orbit scenarios and evaluate the required propagation accuracy for your
% mission's needs.
% =========================================================================

clear; clc; close all;

fprintf('========================================\n');
fprintf(' Three-Way Propagator Comparison\n');
fprintf(' SGP4 vs Custom vs HPOP\n');
fprintf('========================================\n\n');

%% Setup Paths and Constants
addpath('SGP4');
addpath('HPOP');

% MATLAB constants for our propagators
TWOPI = 2*pi;
MINUTES_PER_DAY = 1440;
MINUTES_PER_DAY_SQUARED = MINUTES_PER_DAY^2;

%% Load HPOP Data
% This function loads all necessary data files for HPOP
[hpop_available] = load_hpop_data();
if ~hpop_available
    fprintf('Warning: HPOP data not loaded. HPOP will be skipped.\n\n');
end

%% TLE Definition for ISS
tle_line1 = '1 25544U 98067A   25269.21152988  .00014701  00000-0  26152-3 0  9992';
tle_line2 = '2 25544  51.6332 166.3292 0002995  10.9396 349.1657 15.50384251530832';

%% Parse TLE
% This function parses the TLE strings into a satellite data structure
[satdata] = parse_tle(tle_line1, tle_line2);

fprintf('Satellite TLE Parameters:\n');
fprintf(' Epoch: 20%s, day %.8f\n', tle_line1(19:20), str2double(tle_line1(21:32)));
fprintf(' Inclination: %.4f deg\n', str2double(tle_line2(9:16)));
fprintf(' RAAN: %.4f deg\n', str2double(tle_line2(18:25)));
fprintf(' Eccentricity: %.7f\n', str2double(['0.' tle_line2(27:33)]));
fprintf(' Mean Motion: %.8f rev/day\n', satdata.no_kozai);
fprintf(' B*: %.6e\n\n', satdata.bstar);

%% Propagation Settings
propagation_hours = 24;
time_step_minutes = 2;
tsince_vec = 0:time_step_minutes:(propagation_hours*60);
n_steps = length(tsince_vec);

fprintf('Propagation: %.0f hours in %d steps\n\n', propagation_hours, n_steps);

%% 1. SGP4 Propagation
fprintf('1. Running SGP4 propagation...\n');
r_sgp4 = zeros(3, n_steps);
v_sgp4 = zeros(3, n_steps);

tic;
for idx = 1:n_steps
    [rteme, vteme] = sgp4(tsince_vec(idx), satdata);
    r_sgp4(:, idx) = rteme * 1000; % km to m
    v_sgp4(:, idx) = vteme * 1000; % km/s to m/s
end
t_sgp4 = toc;

fprintf(' ✓ Completed in %.4f s\n\n', t_sgp4);

%% 2. Custom Numerical Propagator
fprintf('2. Running custom numerical propagator (J2 + drag)...\n');

params.mu = 3.986004418e14;
params.R_e = 6378137;
params.J2 = 1.08263e-3;
params.bstar = satdata.bstar;

r0_custom = r_sgp4(:, 1);
v0_custom = v_sgp4(:, 1);

X0_custom = [r0_custom; v0_custom];
X_custom = zeros(6, n_steps);
X_custom(:, 1) = X0_custom;
Xi = X0_custom;

ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);

tic;
for idx = 2:n_steps
    dt = time_step_minutes * 60;
    [~, X_traj] = ode45(@(t, X) two_body_j2_drag(t, X, params), [0, dt], Xi, ode_options);
    Xi = X_traj(end, :)';
    X_custom(:, idx) = Xi;
end
t_custom = toc;

fprintf(' ✓ Completed in %.4f s\n\n', t_custom);


%% 3. HPOP (High Precision Orbit Propagator)
if hpop_available
    fprintf('3. Running HPOP (high-precision numerical propagator)...\n');
    fprintf('   Configuration: 20x20 gravity, Sun/Moon, SRP, drag\n');

    year = str2double(tle_line1(19:20));
    doy = str2double(tle_line1(21:32));
    if year < 57
        full_year = 2000 + year;
    else
        full_year = 1900 + year;
    end

    Mjd0_UTC = Mjday(full_year, 1, 0) + doy;
    Y0_hpop = [r_sgp4(:,1); v_sgp4(:,1)]; % Initial state from SGP4

    % Setup HPOP parameters
    global AuxParam
    AuxParam = struct('Mjd_UTC',Mjd0_UTC,'area_solar',1,'area_drag',1,'mass',12,'Cr',1.5,...
                     'Cd',2.2,'n',20,'m',20,'sun',1,'moon',1,'sRad',1,'drag',1,...
                     'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

    Step = time_step_minutes * 60;
    N_Step = n_steps - 1;

    tic;
    Eph_hpop = Ephemeris(Y0_hpop, N_Step, Step);
    t_hpop = toc;

    r_hpop = Eph_hpop(:, 2:4)';
    v_hpop = Eph_hpop(:, 5:7)';

    fprintf('   ✓ Completed in %.4f s\n\n', t_hpop);
else
    t_hpop = NaN;
    r_hpop = NaN(3, n_steps);
end

%% Final Analysis and Plotting
% This function calculates errors and generates all plots
plot_results(tsince_vec, r_sgp4, X_custom(1:3,:), r_hpop, t_sgp4, t_custom, t_hpop);

fprintf('\n========================================\n');
if hpop_available
    fprintf('Three-Way Benchmark Complete!\n');
else
    fprintf('Two-Way Benchmark Complete!\n');
end
fprintf('========================================\n');


%% Supporting Function for Custom Propagator
function dX = two_body_j2_drag(t, X, params)
    r_vec = X(1:3);
    v_vec = X(4:6);
    r_norm = norm(r_vec);

    % Two-body acceleration
    a_2body = -params.mu * r_vec / r_norm^3;

    % J2 perturbation
    z2 = r_vec(3)^2;
    tx = r_vec(1)/r_norm * (5*z2/r_norm^2 - 1);
    ty = r_vec(2)/r_norm * (5*z2/r_norm^2 - 1);
    tz = r_vec(3)/r_norm * (5*z2/r_norm^2 - 3);
    a_j2 = 1.5 * params.J2 * params.mu * params.R_e^2 / r_norm^4 * [tx; ty; tz];
    
    % Simple atmospheric drag
    rho0 = 1.225; % kg/m^3
    H = 8500; % m
    h = r_norm - params.R_e;
    rho = rho0 * exp(-h/H);
    a_drag = -0.5 * rho * norm(v_vec) * v_vec * (params.bstar); 

    % Total acceleration
    a_total = a_2body + a_j2 + a_drag;
    
    dX = [v_vec; a_total];
end
