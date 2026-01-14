clc
close all

if ~exist('ts', 'var')
    if exist('orbit_data.mat', 'file')
        load('orbit_data.mat');
    else
        fprintf('Data not found. Running orbit_testing...\n');
        orbit_testing;
    end
end

%% Plotting rotation rates
fig = figure(1);

plot(ts/60, omegas, LineWidth=1)

xlabel('Time [min]')
ylabel('Rotation Rates [rad/s]')

legend('\omega_x', '\omega_y', '\omega_z')