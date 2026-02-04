% Test script for computeAeroForces function
% This script demonstrates how to use the consolidated aero force calculator
%
% Author: Generated from ADBSat toolkit
%------------- BEGIN CODE --------------

clear; clc;

% Add ADBSat paths (needed for atmosnrlmsise00)
ADBSat_path = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(ADBSat_path));

%% Define inputs

% 1. Path to .obj file
objFile = fullfile(ADBSat_path, 'inou', 'obj_files', 'cube.obj');

% 2. Location (3D position)
location.altitude = 200e3;      % 200 km in meters
location.latitude = 25.8;       % Geodetic latitude [deg]
location.longitude = 0;         % Longitude [deg]

% 3. Attitude (orientation)
attitude.aoa = 0;               % Angle of attack [deg]
attitude.aos = 0;               % Angle of sideslip [deg]

% Alternative: Use quaternion
% attitude.quaternion = [1, 0, 0, 0]; % Identity quaternion (no rotation)

% 4. Time
time.dayOfYear = 106;           % Day of year (1-365)
time.UTseconds = 0;             % Seconds of the day

% 5. Options (optional - these are the defaults)
options.f107Average = 65;       % 81-day average F10.7 flux
options.f107Daily = 65;         % Daily F10.7 flux
options.magneticIndex = ones(1,7)*3;  % Magnetic indices
options.anomalousOxygen = 0;    % No anomalous oxygen
options.gsi_model = 'cook';     % Gas-Surface Interaction model
options.alpha = 1;              % Accommodation coefficient
options.Tw = 300;               % Wall temperature [K]
options.enableShadow = 1;       % Enable shadow analysis
options.enableSolar = 1;        % Enable solar radiation
options.sol_cR = 0.15;          % Specular reflectivity
options.sol_cD = 0.25;          % Diffuse reflectivity

%% Run the computation
fprintf('Computing aerodynamic forces...\n');
tic;
results = computeAeroForces(objFile, location, attitude, time, options);
elapsed = toc;
fprintf('Computation completed in %.4f seconds\n\n', elapsed);

%% Display results

fprintf('==================== RESULTS ====================\n\n');

fprintf('--- Reference Quantities ---\n');
fprintf('Reference Area:  %.6f m^2\n', results.AreaRef);
fprintf('Projected Area:  %.6f m^2\n', results.AreaProj);
fprintf('Reference Length: %.6f m\n', results.LenRef);
fprintf('\n');

fprintf('--- Atmospheric Conditions ---\n');
fprintf('Density (rho):    %.6e kg/m^3\n', results.rho);
fprintf('Velocity (V_inf): %.2f m/s\n', results.vinf);
fprintf('\n');

fprintf('--- Force Coefficients ---\n');
fprintf('Wind Frame (Cf_w):   [%.6f, %.6f, %.6f]\n', results.Cf_w);
fprintf('Flight Frame (Cf_f): [%.6f, %.6f, %.6f]\n', results.Cf_f);
fprintf('Body Frame (Cf_b):   [%.6f, %.6f, %.6f]\n', results.Cf_b);
fprintf('\n');

fprintf('--- Moment Coefficients (Body Frame) ---\n');
fprintf('Cm_B: [%.6f, %.6f, %.6f]\n', results.Cm_B);
fprintf('\n');

fprintf('--- Aerodynamic Forces [N] (Wind Frame) ---\n');
fprintf('F_aero: [%.6e, %.6e, %.6e] N\n', results.F_aero);
fprintf('|F_aero|: %.6e N\n', norm(results.F_aero));
fprintf('\n');

fprintf('--- Aerodynamic Moments [N*m] (Body Frame) ---\n');
fprintf('M_aero: [%.6e, %.6e, %.6e] N*m\n', results.M_aero);
fprintf('|M_aero|: %.6e N*m\n', norm(results.M_aero));
fprintf('\n');

if isfield(results, 'Cf_s')
    fprintf('--- Solar Force Coefficients ---\n');
    fprintf('Cf_s: [%.6f, %.6f, %.6f]\n', results.Cf_s);
    fprintf('\n');
    
    fprintf('--- Solar Forces [N] (Wind Frame) ---\n');
    fprintf('F_solar: [%.6e, %.6e, %.6e] N\n', results.F_solar);
    fprintf('|F_solar|: %.6e N\n', norm(results.F_solar));
    fprintf('\n');
    
    fprintf('--- Solar Moments [N*m] (Body Frame) ---\n');
    fprintf('M_solar: [%.6e, %.6e, %.6e] N*m\n', results.M_solar);
    fprintf('|M_solar|: %.6e N*m\n', norm(results.M_solar));
end

fprintf('=================================================\n');

%% Example: Sweep angle of attack

fprintf('\n\nRunning AoA sweep...\n');
aoa_range = -30:5:30;  % degrees
Cd_values = zeros(size(aoa_range));
Cl_values = zeros(size(aoa_range));

for i = 1:length(aoa_range)
    attitude.aoa = aoa_range(i);
    attitude.aos = 0;
    res = computeAeroForces(objFile, location, attitude, time, options);
    Cd_values(i) = -res.Cf_w(1);  % Drag is negative x in wind frame
    Cl_values(i) = res.Cf_w(3);   % Lift is z in wind frame
end

fprintf('AoA Sweep Complete!\n');
fprintf('AoA [deg]\tCd\t\tCl\n');
fprintf('----------------------------------------\n');
for i = 1:length(aoa_range)
    fprintf('%+6.1f\t\t%.4f\t\t%.4f\n', aoa_range(i), Cd_values(i), Cl_values(i));
end

%------------- END OF CODE --------------
