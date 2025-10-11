%--------------------------------------------------------------------------
%
% Density_Jacchia70: Computes the atmospheric density
%
% Inputs:
%   r_Sun       Geocentric equatorial position of the Sun (in [m]) referred
%               to the mean equator and equinox of J2000 (EME2000, ICRF)
%   Mjd_UTC     Modified Julian Date UTC
%   r_ECEF      Satellite position vector in ECEF system [m]
%   gast        Greenwich apparent sidereal time [rad]
%
% Output:
%   dens        Density [kg/m^3]
%
% Last modified:   2022/11/19   Meysam Mahooti
%
%--------------------------------------------------------------------------
function dens = Density_Jacchia70(r_Sun,Mjd_UTC,r_ECEF,gast)

global swdata

% Define Location and Altitude
[lon, lat, height] = Geodetic(r_ECEF);

j70iniz;

indata(1) = height/1000;
indata(2) = lat;
indata(3) = lon;

% Define Time of Interest
[year, month, day, hour, minute, sec] = invjday(Mjd_UTC);
indata(4) = year;
indata(5) = month;
indata(6) = day;
indata(7) = hour;
indata(8) = minute;
indata(9) = 2; % geomagnetic index type(1 = indata(12) is Kp, 2 = indata(12) is Ap)
i = find((year==swdata(1,:)) & (month==swdata(2,:)) & (day==swdata(3,:)),1,'first');
sw = swdata(:,i);
sw_1 = swdata(:,i-1);
sw_2 = swdata(:,i-81);

% Define Solar Flux Values
if (indata(9) == 1)
    if hour < 3
        Kp = sw(6);
    elseif hour < 6
        Kp = sw(7);
    elseif hour < 9
        Kp = sw(8);
    elseif hour < 12
        Kp = sw(9);
    elseif hour < 15
        Kp = sw(10);
    elseif hour < 18
        Kp = sw(11);
    elseif hour < 21
        Kp = sw(12);
    else
        Kp = sw(13);
    end
    indata(12) = Kp; % geomagnetic activity index    
else    
    if hour < 3
        Ap = sw(15);
    elseif hour < 6
        Ap = sw(16);
    elseif hour < 9
        Ap = sw(17);
    elseif hour < 12
        Ap = sw(18);
    elseif hour < 15
        Ap = sw(19);
    elseif hour < 18
        Ap = sw(20);
    elseif hour < 21
        Ap = sw(21);
    else
        Ap = sw(22);
    end
    indata(12) = Ap;              % geomagnetic activity index
end
indata(10) = sw_1(30);            % daily F10.7 flux for previous day (observed).
indata(11) = (sw(28)+sw_2(28))/2; % Centered 162-day arithmetic average of F10.7 (adjusted).
indata(13) = gast;
outdata = jatmos70(r_Sun, indata);
dens = outdata(10);
%  outdata(1)  = exospheric temperature (deg K)
%  outdata(2)  = temperature at altitude (deg K)
%  outdata(3)  = N2 number density (per meter-cubed)
%  outdata(4)  = O2 number density (per meter-cubed)
%  outdata(5)  = O number density (per meter-cubed)
%  outdata(6)  = A number density (per meter-cubed)
%  outdata(7)  = He number density (per meter-cubed)
%  outdata(8)  = H number density (per meter-cubed)
%  outdata(9)  = average molecular weight
%  outdata(10) = total density (kilogram/meter-cubed)
%  outdata(11) = log10(total density)
%  outdata(12) = total pressure (pascals)

