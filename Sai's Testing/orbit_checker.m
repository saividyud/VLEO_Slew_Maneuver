% clear
clc
close all

num = length(ts);

%% Testing radial distance error (should be constant for a circular orbit)
r_norms = (vecnorm(rs, 2, 2) - a) / a;

fig = figure(1);

plot(ts, r_norms);

xlabel('Time [s]')
ylabel('Error from Semimajor Axis')
title('Normalized Radial Error');
grid on;

%% Testing rate of change of orbital elements
as = zeros(1, num);
es = zeros(1, num);
is = zeros(1, num);
raans = zeros(1, num);
aops = zeros(1, num);
tas = zeros(1, num);

for i = 1 : 1 : length(r_norms)

    orbit = OEfromRV(rs(i, :), vs(i, :));

    as(i) = orbit(1);
    es(i) = orbit(2);
    is(i) = orbit(3);
    raans(i) = orbit(4);
    aops(i) = orbit(5);
    tas(i) = orbit(6);

end

fig = figure(2);

plot(ts, raans)


