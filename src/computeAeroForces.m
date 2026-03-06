function [results] = computeAeroForces(objFilePath, location, attitude, time, options)
% COMPUTEAEROFORCES Computes aerodynamic and solar radiation forces/moments
%
% This function consolidates the ADBSat workflow to compute total forces and
% moments without saving intermediate files. It directly returns the results
% as a structure.
%
% Inputs:
%   objFilePath : Path to .obj mesh file (string)
%   location    : Structure with fields:000000
%                   .altitude   - Altitude [m]
%                   .latitude   - Geodetic latitude [deg]
%                   .longitude  - Longitude [deg]
%   attitude    : Structure with fields:
%                   .aoa        - Angle of attack [deg]
%                   .aos        - Angle of sideslip [deg]
%                 OR (alternative representation):
%                   .quaternion - [q0, q1, q2, q3] attitude quaternion (body to wind)
%   time        : Structure with fields:
%                   .dayOfYear  - Day of year (1-365)
%                   .UTseconds  - Seconds of the day (0-86400)
%   options     : (Optional) Structure with fields:
%                   .f107Average    - 81-day average F10.7 flux (default: 150)
%                   .f107Daily      - Daily F10.7 flux (default: 150)
%                   .magneticIndex  - 7-element Ap array (default: ones(1,7)*3)
%                   .anomalousOxygen- Flag for anomalous oxygen (default: 0)
%                   .gsi_model      - GSI model name (default: 'cook')
%                   .alpha          - Accommodation coefficient (default: 1)
%                   .Tw             - Wall temperature [K] (default: 300)
%                   .enableShadow   - Enable shadow analysis (default: 1)
%                   .enableSolar    - Enable solar radiation (default: 1)
%                   .sol_cR         - Specular reflectivity (default: 0.15)
%                   .sol_cD         - Diffuse reflectivity (default: 0.25)
%
% Outputs:.F_aero     - Aerodynamic force vector [N] in wind frame
%                   .M_aero     - Aerodynamic moment vector [N*m] in body frame
%   results     : Structure with fields:
%                   .Cf_w       - Force coefficient in wind frame [3x1]
%                   .Cf_f       - Force coefficient in flight frame [3x1]
%                   .Cf_b       - Force coefficient in body frame [3x1]
%                   .Cm_B       - Moment coefficient in body frame [3x1]
%                   .Cf_s       - Solar force coefficient (if solar enabled) [3x1]
%                   .Cm_S       - Solar moment coefficient (if solar enabled) [3x1]
%                   .AreaRef    - Reference area [m^2]
%                   .AreaProj   - Projected area [m^2]
%                   .LenRef     - Reference length [m]
%                   .rho        - Atmospheric density [kg/m^3]
%                   .vinf       - Free-stream velocity [m/s]
%                   
%                   .F_solar    - Solar force vector [N] in wind frame (if enabled)
%                   .M_solar    - Solar moment vector [N*m] in body frame (if enabled)
%                   .param_eq   - Full parameter structure used in calculations
%
% Example:
%   objFile = 'path/to/cube.obj';
%   location.altitude = 200e3;  % 200 km
%   location.latitude = 25.8;
%   location.longitude = 0;
%   attitude.aoa = 0;
%   attitude.aos = 0;
%   time.dayOfYear = 106;
%   time.UTseconds = 0;
%   results = computeAeroForces(objFile, location, attitude, time);
%
% Author: Generated from ADBSat toolkit
% Based on work by David Mostaza Prieto, Nicholas H. Crisp, Luciana Sinpetru, Sabrina Livadiotti
% The University of Manchester
%
%------------- BEGIN CODE --------------

% Set default options
if nargin < 5
    options = struct();
end

% Default environment parameters
if ~isfield(options, 'f107Average'),     options.f107Average = 150; end
if ~isfield(options, 'f107Daily'),       options.f107Daily = 150; end
if ~isfield(options, 'magneticIndex'),   options.magneticIndex = ones(1,7)*3; end
if ~isfield(options, 'anomalousOxygen'), options.anomalousOxygen = 0; end

% Default GSI parameters
if ~isfield(options, 'gsi_model'),       options.gsi_model = 'cook'; end
if ~isfield(options, 'alpha'),           options.alpha = 1; end
if ~isfield(options, 'Tw'),              options.Tw = 300; end

% Default simulation flags
if ~isfield(options, 'enableShadow'),    options.enableShadow = 1; end
if ~isfield(options, 'enableSolar'),     options.enableSolar = 1; end

% Default solar parameters
if ~isfield(options, 'sol_cR'),          options.sol_cR = 0.15; end
if ~isfield(options, 'sol_cD'),          options.sol_cD = 0.25; end

%% Step 1: Import mesh directly (without saving to file)
meshdata = importMeshDirect(objFilePath);

%% Step 2: Setup parameter structure
param_eq.gsi_model = options.gsi_model;
param_eq.alpha = options.alpha;
param_eq.Tw = options.Tw;
param_eq.sol_cR = options.sol_cR;
param_eq.sol_cD = options.sol_cD;

%% Step 3: Compute environment parameters
param_eq = computeEnvironment(param_eq, location.altitude, location.latitude, ...
    location.longitude, time.dayOfYear, time.UTseconds, ...
    options.f107Average, options.f107Daily, options.magneticIndex, ...
    options.anomalousOxygen);

%% Step 4: Convert attitude to angles of attack/sideslip
if isfield(attitude, 'quaternion')
    % Convert quaternion to aoa/aos
    [aoa, aos] = quaternionToAeroAngles(attitude.quaternion);
else
    aoa = deg2rad(attitude.aoa);
    aos = deg2rad(attitude.aos);
end

%% Step 5: Calculate coefficients
results = calculateCoefficients(meshdata, aoa, aos, param_eq, ...
    options.enableShadow, options.enableSolar);

%% Step 6: Compute actual forces and moments
% Dynamic pressure: q = 0.5 * rho * V^2
rho_total = param_eq.rho(6); % Total mass density [kg/m^3]
q = 0.5 * rho_total * param_eq.vinf^2;

% Aerodynamic forces: F = q * Aref * Cf
results.F_aero = q * results.AreaRef * results.Cf_w;

% Aerodynamic moments: M = q * Aref * Lref * Cm
results.M_aero = q * results.AreaRef * results.LenRef * results.Cm_B;

% Solar forces (if enabled)
if options.enableSolar
    % Solar pressure at 1 AU ~ 4.56e-6 N/m^2
    P_solar = 4.56e-6;
    results.F_solar = P_solar * results.AreaRef * results.Cf_s;
    results.M_solar = P_solar * results.AreaRef * results.LenRef * results.Cm_S;
end

% Store additional info
results.rho = rho_total;
results.vinf = param_eq.vinf;
results.param_eq = param_eq;

end

%% ======================== HELPER FUNCTIONS ========================

function meshdata = importMeshDirect(objFilePath)
% Import mesh from .obj file without saving to disk
% Based on importobjtri.m and obj_fileTri2patch.m

    [V, F, X, Y, Z, M] = readObjFile(objFilePath);
    [surfN, areas, bariC] = computeSurfaceNormals(X, Y, Z);
    
    Lref = max(max(X)) - min(min(X));
    
    meshdata.XData = X;
    meshdata.YData = Y;
    meshdata.ZData = Z;
    meshdata.MatID = M;
    meshdata.Areas = areas;
    meshdata.SurfN = surfN;
    meshdata.BariC = bariC;
    meshdata.Lref = Lref;
end

function [V, F, X, Y, Z, M] = readObjFile(fileIn)
% Read .obj file and extract vertices/faces
% Based on obj_fileTri2patch.m

    V = zeros(0,3);
    F = zeros(0,3);
    M = zeros(0,1);
    vertex_index = 1;
    face_index = 1;
    mat_id = 1;
    
    fid = fopen(fileIn, 'r');
    if fid == -1
        error('Cannot open file: %s', fileIn);
    end
    
    line = fgets(fid);
    
    while ischar(line)
        vertex = sscanf(line, 'v %f %f %f');
        face = sscanf(line, 'f %d %d %d');
        face_long = sscanf(line, 'f %d/%d/%d %d/%d/%d %d/%d/%d');
        face_long2 = sscanf(line, 'f %d/%d %d/%d %d/%d');
        face_long3 = sscanf(line, 'f %d//%d %d//%d %d//%d');
        material = sscanf(line, 'usemtl %s;');
        
        if size(vertex) > 0
            V(vertex_index,:) = vertex;
            vertex_index = vertex_index + 1;
        elseif size(material) > 0
            mat_id = mat_id + 1;
        elseif size(face) > 0
            facelength = [length(face), length(face_long), length(face_long2), length(face_long3)];
            [~, indFfor] = max(facelength);
            
            switch indFfor
                case 1
                    F(face_index,:) = face;
                case 2
                    face_long = face_long(1:3:end);
                    F(face_index,:) = face_long;
                case 3
                    face_long = face_long2(1:2:end);
                    F(face_index,:) = face_long;
                case 4
                    face_long = face_long3(1:2:end);
                    F(face_index,:) = face_long;
            end
            M(face_index,:) = mat_id;
            face_index = face_index + 1;
        end
        
        line = fgets(fid);
    end
    fclose(fid);
    
    % Convert to patch structure
    FI = F';
    X = [V(FI(1,:),1)'; V(FI(2,:),1)'; V(FI(3,:),1)'];
    Y = [V(FI(1,:),2)'; V(FI(2,:),2)'; V(FI(3,:),2)'];
    Z = [V(FI(1,:),3)'; V(FI(2,:),3)'; V(FI(3,:),3)'];
    M = M';
end

function [surfN, areas, bariC] = computeSurfaceNormals(x, y, z)
% Calculate surface normals, areas, and barycenters
% Based on surfaceNormals.m

    xV1 = x(2,:) - x(1,:);
    xV2 = x(3,:) - x(1,:);
    yV1 = y(2,:) - y(1,:);
    yV2 = y(3,:) - y(1,:);
    zV1 = z(2,:) - z(1,:);
    zV2 = z(3,:) - z(1,:);
    
    V1 = [xV1; yV1; zV1];
    V2 = [xV2; yV2; zV2];
    
    % Barycenters
    bariC = 1/3 .* [x(1,:)' + x(2,:)' + x(3,:)', ...
                    y(1,:)' + y(2,:)' + y(3,:)', ...
                    z(1,:)' + z(2,:)' + z(3,:)'];
    bariC = bariC';
    
    % Normals
    surfN = cross(V1, V2);
    ModN = sqrt(dot(surfN, surfN));
    surfN = surfN ./ [ModN; ModN; ModN];
    
    % Areas
    areas = 0.5 .* ModN;
end

function param_eq = computeEnvironment(param_eq, h, lat, lon, dayOfYear, UTseconds, ...
    f107Average, f107Daily, magneticIndex, AnO)
% Compute atmospheric parameters using NRLMSISE-00
% Based on environment.m

    % Constants
    data = getConstants();
    
    if AnO
        Oflag = 'Oxygen';
    else
        Oflag = 'NoOxygen';
    end
    
    % Atmospheric properties using MATLAB's atmosnrlmsise00
    [T, param_eq.rho] = atmosnrlmsise00(h, lat, lon, 2002, dayOfYear, UTseconds, ...
        (UTseconds/3600 + lon/15), f107Average, f107Daily, magneticIndex, Oflag);
    
    % Temperature data
    param_eq.Texo = T(1);
    param_eq.Tinf = T(2);
    
    % Calculate mean molecular mass [g mol^-1]
    if AnO
        ND_T = (param_eq.rho(1) + param_eq.rho(2) + param_eq.rho(3) + param_eq.rho(4) + ...
                param_eq.rho(5) + param_eq.rho(7) + param_eq.rho(8) + param_eq.rho(9));
        param_eq.massConc(1,8) = param_eq.rho(8)/param_eq.rho(6) * (data.mAnO/data.NA/1000);
        
        param_eq.mmean = (data.mHe * param_eq.rho(1) + data.mO * param_eq.rho(2) + ...
            data.mN2 * param_eq.rho(3) + data.mO2 * param_eq.rho(4) + ...
            data.mAr * param_eq.rho(5) + data.mH * param_eq.rho(7) + ...
            data.mN * param_eq.rho(8) + data.mAnO * param_eq.rho(9)) / ND_T;
    else
        ND_T = (param_eq.rho(1) + param_eq.rho(2) + param_eq.rho(3) + param_eq.rho(4) + ...
                param_eq.rho(5) + param_eq.rho(7) + param_eq.rho(8));
        
        param_eq.mmean = (data.mHe * param_eq.rho(1) + data.mO * param_eq.rho(2) + ...
            data.mN2 * param_eq.rho(3) + data.mO2 * param_eq.rho(4) + ...
            data.mAr * param_eq.rho(5) + data.mH * param_eq.rho(7) + ...
            data.mN * param_eq.rho(8)) / ND_T;
    end
    
    param_eq.massConc(1,1) = param_eq.rho(1)/param_eq.rho(6) * (data.mHe/data.NA/1000);
    param_eq.massConc(1,2) = param_eq.rho(2)/param_eq.rho(6) * (data.mO/data.NA/1000);
    param_eq.massConc(1,3) = param_eq.rho(3)/param_eq.rho(6) * (data.mN2/data.NA/1000);
    param_eq.massConc(1,4) = param_eq.rho(4)/param_eq.rho(6) * (data.mO2/data.NA/1000);
    param_eq.massConc(1,5) = param_eq.rho(5)/param_eq.rho(6) * (data.mAr/data.NA/1000);
    param_eq.massConc(1,6) = param_eq.rho(7)/param_eq.rho(6) * (data.mH/data.NA/1000);
    param_eq.massConc(1,7) = param_eq.rho(8)/param_eq.rho(6) * (data.mN/data.NA/1000);
    
    % Specific gas constant [J kg^-1 K^-1]
    param_eq.Rmean = (data.R / param_eq.mmean) * 1000;
    
    % Orbital velocity [m/s]
    param_eq.vinf = sqrt(data.mu_E / (data.R_E + h));
    
    % Thermal velocity
    param_eq.vth = sqrt(2 * data.kb * param_eq.Tinf / (param_eq.mmean / data.NA / 1000));
    
    % Speed ratio
    param_eq.s = param_eq.vinf / param_eq.vth;
end

function data = getConstants()
% Return physical constants
% Based on ADBSatConstants.m

    data.mu_E = 3.986004418e14;  % GM Earth [m^3/s^2]
    data.R_E = 6.37813649e6;     % Equatorial Earth radius [m]
    data.R = 8.31446261815324;   % Molar Gas constant [J K^-1 mol^-1]
    data.kb = 1.3806503e-23;     % Boltzmann constant [m^2 kg^-2 K^-1]
    data.NA = 6.02214076e23;     % Avogadro constant [n mol^-1]
    
    data.mHe = 4.002602;         % Molecular mass of He [g mol^-1]
    data.mO = 15.9994;           % Molecular mass of O [g mol^-1]
    data.mN2 = 28.0134;          % Molecular mass of N2 [g mol^-1]
    data.mO2 = 31.9988;          % Molecular mass of O2 [g mol^-1]
    data.mAr = 39.948;           % Molecular mass of Ar [g mol^-1]
    data.mH = 1.0079;            % Molecular mass of H [g mol^-1]
    data.mN = 14.0067;           % Molecular mass of N [g mol^-1]
    data.mAnO = 15.99;           % Molecular mass of anomalous oxygen [g mol^-1]
end

function [aoa, aos] = quaternionToAeroAngles(q)
% Convert attitude quaternion to angle of attack and sideslip
% q = [q0, q1, q2, q3] where q0 is scalar part
% Assumes quaternion represents rotation from body to wind frame

    R_body2wind = quat2dcm(q);
    
    % Extract aoa and aos from DCM
    % Using standard aerospace convention
    aoa = atan2(R_body2wind(3, 1), R_body2wind(1, 1));
    aos = asin(-R_body2wind(2, 1));
end

function results = calculateCoefficients(meshdata, aoa, aos, param_eq, flag_shad, flag_sol)
% Calculate aerodynamic and solar coefficients
% Based on calc_coeff.m

    x = meshdata.XData;
    y = meshdata.YData;
    z = meshdata.ZData;
    areas = meshdata.Areas;
    surfN = meshdata.SurfN;
    barC = meshdata.BariC;
    LenRef = meshdata.Lref;
    matID = meshdata.MatID;
    
    % Rotation matrices
    L_wb = [cos(aos)*cos(aoa), sin(aos), sin(aoa)*cos(aos); ...
            -sin(aos)*cos(aoa), cos(aos), -sin(aoa)*sin(aos); ...
            -sin(aoa), 0, cos(aoa)];
    
    L_gb = [1 0 0; 0 -1 0; 0 0 -1]; % Body to Geometric
    L_gw = L_gb * L_wb';             % Wind to Geometric
    L_fb = [-1 0 0; 0 1 0; 0 0 -1]; % Body to Flight
    
    % Flow direction
    vdir = L_gw * [-1; 0; 0];
    vdir = vdir / norm(vdir);
    
    % Surface normals and angles between faces and flow
    vMatrix = [vdir(1) * ones(1, length(surfN(1,:))); ...
               vdir(2) * ones(1, length(surfN(1,:))); ...
               vdir(3) * ones(1, length(surfN(1,:)))];
    
    % Angles between flow and normals
    delta = real(acos(dot(-vMatrix, surfN)));
    
    uD = vMatrix; % Unit drag vector
    uL = -cross(cross(uD, surfN), uD) ./ vecnorm(cross(cross(uD, surfN), uD));
    
    % Undefined uL panels (plates normal to flow)
    col = find(all(isnan(uL), 1));
    uL(:, col) = -surfN(:, col);
    
    % Dot products for GSI model
    param_eq.gamma = dot(-uD, surfN);
    param_eq.ell = dot(-uL, surfN);
    
    % Local flat plate coefficients
    [cp, ctau, cd, cl, param_eq] = computeMainCoeff(param_eq, delta, matID);
    
    % Solar coefficients
    if flag_sol
        [cn, cs] = computeSolarCoeff(delta, param_eq);
    end
    
    % Backwards facing panels
    areaB = areas;
    areaB(delta * 180/pi > 90) = 0;
    
    shadow = zeros(size(areas));
    shadow(areaB == 0) = 1;
    
    % Shadow analysis
    if flag_shad
        shadPan = computeShadow(x, y, z, barC, delta, L_gw, surfN);
        
        cp(shadPan) = 0;
        ctau(shadPan) = 0;
        cd(shadPan) = 0;
        cl(shadPan) = 0;
        areaB(shadPan) = 0;
        
        if flag_sol
            cn(shadPan) = 0;
            cs(shadPan) = 0;
        end
        shadow(shadPan) = 0.5;
    end
    
    % Areas
    AreaProj = areaB * cos(delta)';
    AreaT = sum(areas);
    AreaRef = AreaT / 2;
    
    % Shear direction
    tauDir = cross(surfN, cross(vMatrix, surfN));
    ModT = sqrt(dot(tauDir, tauDir));
    tauDir = tauDir ./ [ModT; ModT; ModT];
    
    TF = isnan(tauDir(1,:));
    tauDir(:, TF) = 0;
    
    % Force coefficients
    Cf_g = 1/(AreaRef) .* (tauDir * (ctau' .* areas') - surfN * (cp' .* areas'));
    Cf_w = L_gw' * Cf_g;
    Cf_f = L_fb * L_gb' * Cf_g;
    Cf_b = L_gb' * Cf_g;
    
    % Moment coefficients
    Cmom_G = 1/(AreaRef * LenRef) .* (cross(barC, tauDir) * (ctau' .* areas') + ...
             cross(barC, -surfN) * (cp' .* areas'));
    Cm_B = L_gb' * Cmom_G;
    
    % Store results
    results.aoa = aoa;
    results.aos = aos;
    results.Cf_w = Cf_w;
    results.Cf_f = Cf_f;
    results.Cf_b = Cf_b;
    results.Cm_B = Cm_B;
    results.AreaRef = AreaRef;
    results.AreaProj = AreaProj;
    results.LenRef = LenRef;
    
    % Solar coefficients
    if flag_sol
        cstau = cs .* sin(delta);
        cs_n = cs .* cos(delta);
        Cs_G = 1/(AreaRef) .* (tauDir * (cstau' .* areas') - surfN * ((cn + cs_n)' .* areas'));
        Cf_s = L_gw' * Cs_G;
        Cms_G = 1/(AreaRef * LenRef) .* (cross(barC, tauDir) * (cstau' .* areas') + ...
                cross(barC, -surfN) * ((cn + cs_n)' .* areas'));
        Cm_S = L_gb' * Cms_G;
        
        results.Cf_s = Cf_s;
        results.Cm_S = Cm_S;
    end
end

function [cp, ctau, cd, cl, param_eq] = computeMainCoeff(param_eq, delta, matID)
% Calculate main GSI coefficients
% Based on mainCoeff.m

    % Check material references
    nMat = max(matID);
    if nMat == 0
        matID(matID == 0) = 1;
        nMat = 1;
    end
    
    % Apply material properties per mesh face
    if isfield(param_eq, 'alpha')
        nParam = numel(param_eq.alpha);
        if nParam > 1 && isequal(nMat, nParam)
            param_eq.alpha = param_eq.alpha(matID);
        elseif nParam == 1
            param_eq.alpha = param_eq.alpha * ones(size(matID));
        end
    end
    
    % Calculate GSI based on model
    switch param_eq.gsi_model
        case 'cook'
            [cp, ctau, cd, cl] = coeffCook(param_eq, delta);
        otherwise
            % Default to Cook model
            [cp, ctau, cd, cl] = coeffCook(param_eq, delta);
    end
end

function [cp, ctau, cd, cl] = coeffCook(param_eq, delta)
% Cook's hyperthermal FMF coefficients
% Based on coeff_cook.m

    alpha = param_eq.alpha;
    Tw = param_eq.Tw;
    Rmean = param_eq.Rmean;
    Vinf = param_eq.vinf;
    Tinf = Vinf^2 / (3 * Rmean);
    
    cd = 2 .* cos(delta) .* (1 + (2/3) * sqrt(1 + alpha .* (Tw/Tinf - 1)) .* cos(delta));
    cl = 4/3 .* sqrt(1 + alpha .* (Tw/Tinf - 1)) .* sin(delta) .* cos(delta);
    
    ind = find(delta >= pi/2);
    cd(ind) = 0;
    cl(ind) = 0;
    
    cp = cd .* cos(delta) + cl .* sin(delta);
    ctau = cd .* sin(delta) - cl .* cos(delta);
end

function [cn, cs] = computeSolarCoeff(delta, param_eq)
% Solar radiation pressure coefficients
% Based on coeff_solar.m

    cn = 2 * ((param_eq.sol_cD/3) * cos(delta) + param_eq.sol_cR * cos(delta).^2);
    cs = (1 - param_eq.sol_cR) * cos(delta);
    cs(cs < 0) = 0;
end

function shadPan = computeShadow(x, y, z, barC, delta, L_gw, surfN)
% Shadow analysis for mesh elements
% Based on shadowAnaly.m

    pAw = L_gw' * [x(1,:); y(1,:); z(1,:)];
    pBw = L_gw' * [x(2,:); y(2,:); z(2,:)];
    pCw = L_gw' * [x(3,:); y(3,:); z(3,:)];
    
    barCw = L_gw' * barC;
    
    xW = [pAw(1,:); pBw(1,:); pCw(1,:)];
    
    xWmax = max(xW);
    xWmin = min(xW);
    
    indB = find(delta * 180/pi > 90.0001);
    indF = find(delta * 180/pi <= 90.0001);
    
    if isempty(indB) || isempty(indF)
        shadPan = [];
        return;
    end
    
    minXwF = min(xWmin(indF));
    maxXwB = max(xWmin(indB));
    
    indFPot = find(xWmin(indF) - maxXwB < 10e-5);
    indBPot = find(xWmax(indB) - minXwF > 10e-5);
    
    if isempty(indFPot) || isempty(indBPot)
        shadPan = [];
        return;
    end
    
    yWC = [pAw(2, indB(indBPot)); pBw(2, indB(indBPot)); pCw(2, indB(indBPot))];
    zWC = [pAw(3, indB(indBPot)); pBw(3, indB(indBPot)); pCw(3, indB(indBPot))];
    
    shadPanVec = zeros(1, length(barCw));
    
    tolB = 1e-5;
    
    for i = 1:length(indFPot)
        transYcord = yWC - barCw(2, indF(indFPot(i)));
        yCoord = abs(sum(sign(transYcord)));
        yC_chang = find(yCoord < 3);
        
        if isempty(yC_chang)
            continue;
        else
            transZcord = zWC(:, yC_chang) - barCw(3, indF(indFPot(i)));
            zCoord = abs(sum(sign(transZcord)));
            zC_chang = find(zCoord < 3);
            
            if isempty(zC_chang)
                continue;
            else
                p1 = [transYcord(1, yC_chang(zC_chang)); transZcord(1, zC_chang)];
                p2 = [transYcord(2, yC_chang(zC_chang)); transZcord(2, zC_chang)];
                p3 = [transYcord(3, yC_chang(zC_chang)); transZcord(3, zC_chang)];
                
                Ftest = checkInsideTri(p1, p2, p3, zeros(2, length(zC_chang)));
                
                if isempty(Ftest)
                    continue;
                else
                    if any(dot(L_gw' * surfN(:, indB(indBPot(yC_chang(zC_chang(Ftest))))), ...
                            barCw(:, indF(indFPot(i))) - pAw(:, indB(indBPot(yC_chang(zC_chang(Ftest)))))) > tolB)
                        shadPanVec(indF(indFPot(i))) = 1;
                    end
                end
            end
        end
    end
    
    shadPan = find(shadPanVec > 0);
end

function Ftest = checkInsideTri(p1, p2, p3, p0)
% Check if points are inside triangles
% Based on insidetri.m

    Ftest = [];
    
    for i = 1:size(p1, 2)
        v0 = p3(:,i) - p1(:,i);
        v1 = p2(:,i) - p1(:,i);
        v2 = p0(:,i) - p1(:,i);
        
        dot00 = dot(v0, v0);
        dot01 = dot(v0, v1);
        dot02 = dot(v0, v2);
        dot11 = dot(v1, v1);
        dot12 = dot(v1, v2);
        
        invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;
        
        if (u >= 0) && (v >= 0) && (u + v <= 1)
            Ftest = [Ftest, i];
        end
    end
end

%------------- END OF CODE --------------
