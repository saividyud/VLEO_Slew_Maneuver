clc; clear; close all;

%Polar circular orbit with random initial orientation and initial angular velocity values
Xi = [0, 0, 6678e3, -7789, 0, 0, 20*rand, 20*rand-10, 20*rand-10, 20*rand-10, 2*rand-1, 2*rand-1, 2*rand-1]';
Tf = 40;

%initial state
X = Xi;

%animating orbit with ode45
for i = 1:Tf
    %integrating to next time step
    [t,X] = ode45(@Sat_template,[0, 1],X);
    X = X(end,:)';

    %plotting beta1,beta2 and beta3
    hold on
    figure(1)
    plot(i,X(8),'o');
    plot(i,X(9),'*');
    plot(i,X(10),'.');
    legend('beta1','beta2','beta3');
    title('Attitude of Satellite in Quaternions')
    ylabel('Rad')
    xlabel('Time(s)')
    axis equal

    %pausing for animation realism
    pause(.1)
end