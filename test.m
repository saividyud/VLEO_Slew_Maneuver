clear
clc
close all

[t,y] = ode45(@Sat_template,[0.1000],[0,0,6678e3,-7750,0,0,0,0,0,0,0,0,0]');
plot3(y(:,1),y(:,2),y(:,3));
