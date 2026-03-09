% Example use of computeAeroForces for the 6U CubeSat mesh.

clear;
clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

objFile = fullfile(projectRoot, 'src', 'dynamics', '6U CubeSat.obj');

location.altitude = 200e3;
location.latitude = 25.8;
location.longitude = 0;

attitude.aoa = 0;
attitude.aos = 0;

time.year = 2002;
time.dayOfYear = 106;
time.UTseconds = 0;

commonOptions = struct( ...
    'f107Average', 65, ...
    'f107Daily', 65, ...
    'magneticIndex', ones(1, 7) * 3, ...
    'anomalousOxygen', 0, ...
    'alpha', 1, ...
    'Tw', 300, ...
    'enableShadow', 1, ...
    'enableSolar', 1, ...
    'sol_cR', 0.15, ...
    'sol_cD', 0.25 ...
);

cookOptions = commonOptions;
cookOptions.gsi_model = 'cook';

sentmanOptions = commonOptions;
sentmanOptions.gsi_model = 'sentman';

cookResults = computeAeroForces(objFile, location, attitude, time, cookOptions);
sentmanResults = computeAeroForces(objFile, location, attitude, time, sentmanOptions);

fprintf('6U CubeSat aerodynamic comparison at 200 km\n');
fprintf('------------------------------------------\n');
fprintf('Reference area: %.6f m^2\n', cookResults.AreaRef);
fprintf('Projected area: %.6f m^2\n', cookResults.AreaProj);
fprintf('Density: %.6e kg/m^3\n', cookResults.rho);
fprintf('Velocity: %.2f m/s\n\n', cookResults.vinf);

fprintf('Cook drag coefficient:    %.6f\n', -cookResults.Cf_w(1));
fprintf('Sentman drag coefficient: %.6f\n', -sentmanResults.Cf_w(1));
fprintf('Cook solar coefficient:   %.6f\n', -cookResults.Cf_s(1));
fprintf('Sentman solar coefficient %.6f\n', -sentmanResults.Cf_s(1));

fprintf('\nTip: location.positionECI can also be supplied directly when calling from the nonlinear dynamics models.\n');
