clc; clear; close all;

%Polar circular orbit with simple initial orientation values
Xi = [0 0 6678e3 -7789 0 0 1 1 1 1 0 0 1]';

%propogating orbit with ode15s
[t,X] = ode15s(@Sat_template,[0, 50000],Xi);

%plotting
figure(1)
plot3(X(:,1),X(:,2),X(:,3),'o',MarkerFaceColor=[.5 0 .5]);
figure(2)
plot3(X(:,11),X(:,12),X(:,13),'o',MarkerFaceColor=[0 .5 .5]);
axis equal