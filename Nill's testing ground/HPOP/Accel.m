%--------------------------------------------------------------------------
%
% Accel: Computes the acceleration of an Earth orbiting satellite due to 
%    	 - Earth's harmonic gravity field (including Solid Earth Tides and
%      	   Ocean Tides), 
%    	 - gravitational perturbations of the Sun, Moon and planets
%    	 - solar radiation pressure
%    	 - atmospheric drag and
%	 	 - relativity
%
% Inputs:
%   Mjd_UTC     Modified Julian Date (UTC)
%   Y           Satellite state vector in the ICRF/EME2000 system
%   Area        Cross-section 
%   mass        Spacecraft mass
%   Cr          Radiation pressure coefficient
%   Cd          Drag coefficient
%
% Output:
%   dY          Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2025/02/19   Meysam Mahooti
% 
%--------------------------------------------------------------------------
function dY = Accel(t, Y)

global const AuxParam eopdata

MJD_UTC = AuxParam.Mjd_UTC+t/86400;
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400;

% Form bias-precession-nutation matrix
NPB = iauPnm06a(const.DJM0, MJD_TT);
% Form Earth rotation matrix
gast = iauGst06(const.DJM0, MJD_UT1, const.DJM0, MJD_TT, NPB);
Theta  = iauRz(gast, eye(3));
% Polar motion matrix (TIRS->ITRS, IERS 2003)
Pi = iauPom00(x_pole, y_pole, iauSp00(const.DJM0, MJD_TT));
% ICRS to ITRS transformation
E = Pi*Theta*NPB;

MJD_TDB = Mjday_TDB(MJD_TT);
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE440(MJD_TDB+2400000.5);

% Acceleration due to harmonic gravity field
if (AuxParam.SolidEarthTides || AuxParam.OceanTides)
    a = AccelHarmonic_AnelasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole*const.Arcs,y_pole*const.Arcs);
    % a = AccelHarmonic_ElasticEarth(MJD_UTC,r_Sun,r_Moon,Y(1:3),E,UT1_UTC,TT_UTC,x_pole*const.Arcs,y_pole*const.Arcs);
else
    a = AccelHarmonic(Y(1:3), E, AuxParam.n, AuxParam.m);
end

% Luni-solar perturbations
if (AuxParam.sun)
    a = a + AccelPointMass(Y(1:3),r_Sun,const.GM_Sun);
end

if (AuxParam.moon)
    a = a + AccelPointMass(Y(1:3),r_Moon,const.GM_Moon);
end

% Planetary perturbations
if (AuxParam.planets)
    a = a + AccelPointMass(Y(1:3),r_Mercury,const.GM_Mercury);
    a = a + AccelPointMass(Y(1:3),r_Venus,const.GM_Venus);
    a = a + AccelPointMass(Y(1:3),r_Mars,const.GM_Mars);
    a = a + AccelPointMass(Y(1:3),r_Jupiter,const.GM_Jupiter);
    a = a + AccelPointMass(Y(1:3),r_Saturn,const.GM_Saturn);
    a = a + AccelPointMass(Y(1:3),r_Uranus,const.GM_Uranus);    
    a = a + AccelPointMass(Y(1:3),r_Neptune,const.GM_Neptune);
    a = a + AccelPointMass(Y(1:3),r_Pluto,const.GM_Pluto);
end

% Solar radiation pressure
if (AuxParam.sRad)
    a = a + AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB, ...
        AuxParam.area_solar,AuxParam.mass,AuxParam.Cr,const.P_Sol,const.AU,'conical');
end

% Atmospheric drag
if (AuxParam.drag)
    % Atmospheric density
	Omega = const.omega_Earth-0.843994809*1e-9*LOD; % [rad/s]; IERS
    dens = nrlmsise00(MJD_UTC,E*Y(1:3),UT1_UTC,TT_UTC);
    % [~,dens] = JB2008(MJD_UTC,r_Sun,Y(1:3),E);
    % [~,dens] = JB2006(MJD_UTC,r_Sun,Y(1:3),E);
    % [d,~] = msis86(MJD_UTC,E*Y(1:3),gast);
    % dens = 1e3*d(6);
    % dens = Density_Jacchia70(r_Sun,MJD_UTC,E*Y(1:3),gast);
    % dens = Density_HP(r_Sun,NPB*Y(1:3));
    a = a + AccelDrag(dens,Y(1:3),Y(4:6),NPB,AuxParam.area_drag,AuxParam.mass,AuxParam.Cd,Omega);
end

% Relativistic Effects
if (AuxParam.Relativity)
    a = a + Relativity(Y(1:3),Y(4:6));
end

dY = [Y(4:6);a];

