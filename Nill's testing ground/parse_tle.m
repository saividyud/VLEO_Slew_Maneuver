% =========================================================================
% parse_tle.m
%
% Description:
% This function parses two-line element (TLE) set strings to extract the
% orbital elements required by the SGP4 propagator. It converts the raw
% TLE values into the units and format expected by the SGP4 algorithm.
%
% Inputs:
%   line1 (string): The first line of the TLE.
%   line2 (string): The second line of the TLE.
%
% Returns:
%   satdata (struct): A structure containing the parsed and converted
%                     orbital elements for use with SGP4.
% =========================================================================

function [satdata] = parse_tle(line1, line2)

    % Define constants
    TWOPI = 2*pi;
    MINUTES_PER_DAY = 1440;
    MINUTES_PER_DAY_SQUARED = MINUTES_PER_DAY^2;

    % Parse TLE Line 1
    satdata.epochyr = str2double(line1(19:20));
    satdata.epochdays = str2double(line1(21:32));
    satdata.ndot = str2double(line1(34:43));
    satdata.nddot = str2double(line1(45:50)) * 10^str2double(line1(51:52));
    bstar_val = str2double(line1(54:59));
    bstar_exp = str2double(line1(60:61));
    satdata.bstar = bstar_val * 1e-5 * 10^bstar_exp;
    
    % Parse TLE Line 2
    satdata.inclo = str2double(line2(9:16));   % Inclination [deg]
    satdata.nodeo = str2double(line2(18:25));   % RAAN [deg]
    satdata.ecco = str2double(['0.' line2(27:33)]); % Eccentricity
    satdata.argpo = str2double(line2(35:42));   % Argument of Perigee [deg]
    satdata.mo = str2double(line2(44:51));      % Mean Anomaly [deg]
    satdata.no_kozai = str2double(line2(53:63)); % Mean Motion [rev/day]

    % Convert to SGP4 units (radians and minutes)
    satdata.xmo = satdata.mo * (pi/180);
    satdata.xnodeo = satdata.nodeo * (pi/180);
    satdata.omegao = satdata.argpo * (pi/180);
    satdata.xincl = satdata.inclo * (pi/180);
    satdata.eo = satdata.ecco;
    satdata.xno = satdata.no_kozai * TWOPI / MINUTES_PER_DAY;
    satdata.xndt2o = satdata.ndot * TWOPI / MINUTES_PER_DAY_SQUARED;
    satdata.xndd6o = satdata.nddot * TWOPI / (MINUTES_PER_DAY^3);
    
    % The SGP4 library has some different field names
    satdata.epoch = satdata.epochyr*1000 + satdata.epochdays;

end
