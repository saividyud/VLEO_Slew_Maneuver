clc; clear; close all;

%Polar circular orbit with simple initial orientation values

Xi = [0 0 6678e3 -7789 0 0 1 1 1 1 0 0 1]';
Tf = 1000;

%animating orbit with ode15s

X = Xi;
for i = 1:Tf
    [t,X] = ode45(@Sat_template,[0, 1],X);
    X = X(end,:)';
    %plotting
    hold on
    figure(1)
    plot3(X(1),X(2),X(3),'o');
    axis padded
    view(45,45)
end