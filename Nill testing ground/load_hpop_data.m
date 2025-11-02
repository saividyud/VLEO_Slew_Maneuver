% =========================================================================
% load_hpop_data.m
%
% Description:
% This function loads all external data files required by the HPOP
% propagator. It initializes global variables that store constants,
% Earth orientation parameters, space weather data, gravity models,
% and planetary ephemerides.
%
% Returns:
%   success (logical): True if all data files are loaded successfully,
%                      false otherwise.
% =========================================================================

function [] = load_hpop_data()
    global const Cnm Snm AuxParam eopdata swdata SOLdata DTCdata APdata PC
    
    % Initialize constants
    SAT_Const
    constants % HPOP constants
    
    % MATLAB constants for our propagators
    TWOPI = 2*pi;
    MINUTES_PER_DAY = 1440;
    MINUTES_PER_DAY_SQUARED = (MINUTES_PER_DAY * MINUTES_PER_DAY);
    
    % Load HPOP data files
    fprintf('Loading HPOP data files...\n');
    
    % Load JPL ephemeris
    load('HPOP/DE440Coeff.mat');
    PC = DE440Coeff;
    
    % Read Earth gravity field coefficients
    Cnm = zeros(361,361);
    Snm = zeros(361,361);
    fid = fopen('HPOP/GGM03C.txt','r');
    for n=0:360
        for m=0:n
            temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);
            Cnm(n+1,m+1) = temp(3);
            Snm(n+1,m+1) = temp(4);
        end
    end
    fclose(fid);
    
    % Read Earth orientation parameters
    fid = fopen('HPOP/EOP-All.txt','r');
    eopdata = [];
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
    
    % Read space weather data
    fid = fopen('HPOP/SW-All.txt','r');
    swdata = [];
    while ~feof(fid)
        tline = fgetl(fid);
        k = strfind(tline,'NUM_OBSERVED_POINTS');
        if (k ~= 0)
            numrecsobs = str2num(tline(21:end));
            tline = fgetl(fid);
            for i=1:numrecsobs
                swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %i %f %f %f %f %f',[33 1]);
            end
            swdata = [swdata(1:27,:);swdata(29:33,:)];
            for i=1:4
                tline = fgetl(fid);
            end
            numrecspred = str2num(tline(28:end));
            tline = fgetl(fid);
            for i=numrecsobs+1:numrecsobs+numrecspred
                swdata(:,i) = fscanf(fid,'%i %d %d %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %f %i %i %f %f %f %f %f %f',[32 1]);
            end
            break
        end
    end
    fclose(fid);
    
    % Read solar storm indices
    fid = fopen('HPOP/SOLFSMY.txt','r');
    SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
    fclose(fid);
    
    % Read Ap data
    fid = fopen('HPOP/SOLRESAP.txt','r');
    APdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d',[12 inf]);
    fclose(fid);
    
    % Read geomagnetic storm indices
    fid = fopen('HPOP/DTCFILE.txt','r');
    DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
    fclose(fid);
end
