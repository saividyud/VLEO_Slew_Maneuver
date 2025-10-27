clc
close all

%% Plotting rotation rates
fig = figure(1);

plot(ts/60, omegas, LineWidth=1)

xlabel('Time [min]')
ylabel('Rotation Rates [rad/s]')

legend('\omega_x', '\omega_y', '\omega_z')