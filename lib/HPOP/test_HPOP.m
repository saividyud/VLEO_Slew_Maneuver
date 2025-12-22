%--------------------------------------------------------------------------
%
%               High Precision Orbit Propagator
%
% References:
% Montenbruck O., and Gill E.; Satellite Orbits: Models, Methods, and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 4th edition (2013).
%
% McCarthy D. D.; IERS Conventions (1996).
% 
% http://celestrak.org/SpaceData/EOP-All.txt
% http://celestrak.org/SpaceData/SW-All.txt
% https://sol.spacenvironment.net/JB2008/indices/SOLFSMY.TXT
% https://sol.spacenvironment.net/JB2008/indices/SOLRESAP.TXT
% https://sol.spacenvironment.net/JB2008/indices/DTCFILE.TXT
% https://ssd.jpl.nasa.gov/planets/eph_export.html
%
% Last modified:   2025/06/24   Meysam Mahooti
%
%--------------------------------------------------------------------------
clc
clear all
format longG
close all
profile on
tic

global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC

SAT_Const
constants
load DE440Coeff.mat
PC = DE440Coeff;

% read Earth gravity field coefficients
Cnm = zeros(361,361);
Snm = zeros(361,361);
fid = fopen('GGM03C.txt','r');
for n=0:360
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end
fclose(fid);

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

% read space weather data
fid = fopen('SW-All.txt','r');
%  ---------------------------------------------------------------------------------------------------------------------------------
% |                                                                                             Adj     Adj   Adj   Obs   Obs   Obs 
% | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
%  ---------------------------------------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k ~= 0) 
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %i %f %f %f %f %f',[33 1]);
        end
        % remove the row of the Q parameter
        swdata = [swdata(1:27,:);swdata(29:33,:)];
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(28:end));
        tline = fgetl(fid);
        %  -------------------------------------------------------------------------------------------------------------------------------
        % |                                                                                             Adj   Adj   Adj   Obs   Obs   Obs 
        % | yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN F10.7 Ctr81 Lst81 F10.7 Ctr81 Lst81
        %  -------------------------------------------------------------------------------------------------------------------------------
        for i=numrecsobs+1:numrecsobs+numrecspred
            swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %f %f %f %f %f',[32 1]);
        end
        break
    end
end
fclose(fid);

% read solar storm indices
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% read Ap data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d',[12 inf]);
fclose(fid);

% read geomagnetic storm indices
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

% epoch state (Envisat) in ITRF
fid = fopen('InitialState.txt','r');
tline = fgetl(fid);
year = str2num(tline(1:4));
month = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
minute = str2num(tline(15:16));
sec = str2num(tline(18:end));
Y0 = zeros(1,6);
for ii=1:6
    tline = fgetl(fid);
    Y0(ii) = str2num(tline); % [m, m/s]
end
tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:end)); % [m^2]
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:end)); % [m^2]
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:end)); % [kg]
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:end));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:end));
fclose(fid);

% epoch
Mjd0_UTC = Mjday(year, month, day, hour, minute, sec);
Y0 = ECEF2ECI(Mjd0_UTC, Y0);

AuxParam.Mjd_UTC = Mjd0_UTC;
AuxParam.n       = 70;
AuxParam.m       = 70;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;

Step   = 60;   % [s] integration step size
N_Step = 1588; % number of integration steps (26.47 hours)

% propagation
Eph_eci = Ephemeris(Y0, N_Step, Step);

fid = fopen('SatelliteStates.txt','w');
for i=1:N_Step+1
    [year,month,day,hour,minute,sec] = invjday(Mjd0_UTC+Eph_eci(i,1)/86400);
    fprintf(fid,'  %4d/%2.2d/%2.2d  %2d:%2d:%6.3f',year,month,day,hour,minute,sec);
    fprintf(fid,'  %14.3f%14.3f%14.3f%12.3f%12.3f%12.3f\n',...
            Eph_eci(i,2),Eph_eci(i,3),Eph_eci(i,4),Eph_eci(i,5),Eph_eci(i,6),Eph_eci(i,7));
end
fclose(fid);

[n, m] = size(Eph_eci);
Eph_ecef = zeros(n,m);
for i=1:n
    Eph_ecef(i,1) = Eph_eci(i,1);
    Eph_ecef(i,2:7) = ECI2ECEF(Mjd0_UTC+Eph_ecef(i,1)/86400, Eph_eci(i,2:7));    
end

True_EnvisatStates
dd = True_Eph-Eph_ecef(:,2:7);

% Plot orbit in ECI reference
figure(1)
plot3(Eph_eci(:,2),Eph_eci(:,3),Eph_eci(:,4),'o-r')
grid;
title('Orbit ECI (inertial) (m)')

% Plot orbit in ECEF reference
figure(2)
plot3(Eph_ecef(:,2),Eph_ecef(:,3),Eph_ecef(:,4),'-')
title('Orbit ECEF (m)')
xlabel('X');ylabel('Y');zlabel('Z');
grid

% Plot Discrepancy of Precise and Propagated orbits
figure(3)
subplot(3,1,1);
plot(dd(:,1));
title('Discrepancy of Precise and Propagated Envisat Positions for 26.47 hours');
axis tight
xlabel('Time')
ylabel('dX[m]')
hold on
subplot(3,1,2);
plot(dd(:,2));
axis tight
xlabel('Time')
ylabel('dY[m]')
subplot(3,1,3);
plot(dd(:,3));
axis tight
xlabel('Time')
ylabel('dZ[m]')
toc

tic
lamda = zeros(n,1);
phi = zeros(n,1);
height = zeros(n,1);
for i=1:n
    [lamda(i),phi(i),height(i)] = Geodetic(Eph_ecef(i,2:4));
end
figure(4)
geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
title('Satellite''s Ground Track')
hold on
plot(lamda*(180/pi),phi*(180/pi),'.r')
% animation
an = animatedline('Marker','*');
for k = 1:n
    addpoints(an,lamda(k)*(180/pi),phi(k)*(180/pi));
    drawnow
    pause(0.01);
    clearpoints(an);
end
toc

profile viewer
profile off

