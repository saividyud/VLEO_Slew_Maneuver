function Xd = sat_template_gui(t, X, subFigHandle)
% Sat_template calculates the time rate of change of the state X at a time
% t.
%
% Parameters
% ----------
% t : float
%   Time [s]
% X : 13x1 vector
%   State of the system with the following attributes:
%       - Position (1:3) [m]
%       - Velocity (4:6) [m/s]
%       - Quaternion (7:10)
%       - Angular velocity (11:13)
% subFigHandle : figure handle
%   Handle to the figure containing the simulation GUI, used to access
%
% Returns
% -------
% Xd : 13x1 vector
%   Dotted state vector [m/s, m/s^2, -, rad/s^2]

% Xd and X are 13 dimensional state vectors
% X has to be vertical for function to work
% t is time(used for numerical integration)

    %% Satellite configuration (persistent to avoid reloading)
    persistent objFilePath aeroOptions satMass ICB mu simParams day_start seconds_start
    if isempty(objFilePath)
        % Reading in settings from GUI
        simParams = guidata(subFigHandle);

        % Path to satellite mesh file (relative to this file's location)
        thisFilePath = fileparts(mfilename('fullpath'));
        objFilePath = fullfile(thisFilePath, 'dynamics', '6U CubeSat.obj');
        
        % Aerodynamic options
        aeroOptions.f107Average = simParams.initParams.Environment.F107Average;      % 81-day average F10.7 flux
        aeroOptions.f107Daily = simParams.initParams.Environment.F107Daily;        % Daily F10.7 flux
        aeroOptions.magneticIndex = simParams.initParams.Environment.magneticIndices;
        aeroOptions.anomalousOxygen = simParams.initParams.Environment.enableAnomalousOxygen;
        aeroOptions.gsi_model = simParams.initParams.Environment.gasSurfaceInteractionModel; % Gas-surface interaction model
        aeroOptions.alpha = simParams.initParams.Environment.accommodationCoefficient;              % Accommodation coefficient
        aeroOptions.Tw = simParams.initParams.Environment.wallTemperature;               % Wall temperature [K]
        aeroOptions.enableShadow = simParams.initParams.Environment.enableShadowAnalysis;
        aeroOptions.enableSolar = simParams.initParams.Environment.enableSolarRadiationPressure;
        aeroOptions.sol_cR = simParams.initParams.Environment.specularReflectivity;          % Specular reflectivity
        aeroOptions.sol_cD = simParams.initParams.Environment.diffuseReflectivity;          % Diffuse reflectivity
        
        satMass = 83;  % Satellite mass [kg]
        
        % Moment of inertia tensor
        % Can start with approximating a sphere (think Sputnik)
        ICB = 2/5*83*(.58/2)^2*[1 0 0 ; 0 1 0 ;0 0 1]; % [kg m^2]

        % Earth gravitational parameter [m^3/s^2]
        mu = simParams.initParams.Environment.mu;

        % Days of year
        day_start = simParams.initParams.Environment.dayOfYear;

        % Seconds of day
        seconds_start = simParams.initParams.Environment.secondsOfDay;
    end

    %% Compute aerodynamic forces and moments
    % Constants
    R_E = earthRadius;     % Earth equatorial radius [m]
    
    r = norm(X(1:3));       % Distance from Earth center [m]
    
    % Convert ECI position to geodetic coordinates
    altitude = r - R_E;                           % Altitude [m] (spherical approx)
    latitude = asind(X(3)/r);                     % Latitude [deg]
    longitude = atan2d(X(2), X(1));               % Longitude [deg]
    
    % Set up location structure for computeAeroForces
    location.altitude = altitude;
    location.latitude = latitude;
    location.longitude = longitude;
    
    % Compute angle of attack and sideslip from velocity and attitude
    % First, get velocity in body frame using quaternion
    q = X(7:10);  % [q0, q1, q2, q3] - scalar first
    v_eci = X(4:6);
    
    % Rotation matrix from ECI to body frame (using quaternion)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    R_eci2body = [1-2*(q2^2+q3^2),   2*(q1*q2+q0*q3),   2*(q1*q3-q0*q2);
                  2*(q1*q2-q0*q3),   1-2*(q1^2+q3^2),   2*(q2*q3+q0*q1);
                  2*(q1*q3+q0*q2),   2*(q2*q3-q0*q1),   1-2*(q1^2+q2^2)];
    
    % Velocity in body frame
    v_body = R_eci2body * v_eci;
    v_mag = norm(v_body);
    
    % Compute aerodynamic angles (angle of attack and sideslip)
    % Assuming body x-axis is forward
    if v_mag > 1  % Avoid division by zero
        aoa = atan2d(v_body(3), v_body(1));      % Angle of attack [deg]
        aos = atan2d(v_body(2), v_body(1));      % Angle of sideslip [deg]
    else
        aoa = 0;
        aos = 0;
    end
    
    % Set up attitude structure
    attitude.aoa = aoa;
    attitude.aos = aos;
    
    % Set up time structure (using simulation time)
    % Assuming simulation starts at day 106 (mid-April)
    time.dayOfYear = day_start + floor((seconds_start + t)/86400);     % Day of year
    time.UTseconds = seconds_start + mod(t, 86400);              % Seconds of the day
    
    % Compute aerodynamic forces and moments
    try
        aeroResults = computeAeroForces(objFilePath, location, attitude, time, aeroOptions);
        
        % Aerodynamic force in wind frame [N]
        F_aero_wind = aeroResults.F_aero;
        
        % Aerodynamic moment in body frame [N*m]
        M_aero_body = aeroResults.M_aero;
        
        % Transform force from wind frame to body frame
        % Wind frame: x = -drag direction (velocity), z = lift
        aoa_rad = deg2rad(aoa);
        aos_rad = deg2rad(aos);
        
        % Rotation matrix from wind to body frame
        L_bw = [cos(aos_rad)*cos(aoa_rad), -sin(aos_rad)*cos(aoa_rad), -sin(aoa_rad);
                sin(aos_rad),               cos(aos_rad),               0;
                cos(aos_rad)*sin(aoa_rad), -sin(aos_rad)*sin(aoa_rad),  cos(aoa_rad)];
        
        F_aero_body = L_bw * F_aero_wind;
        
        % Transform force from body frame to ECI frame
        R_body2eci = R_eci2body';
        F_aero_eci = R_body2eci * F_aero_body;
        
        % Aerodynamic acceleration [m/s^2]
        a_aero = F_aero_eci / satMass;
        
        % Add solar radiation pressure if available
        if isfield(aeroResults, 'F_solar')
            F_solar_body = L_bw * aeroResults.F_solar;
            F_solar_eci = R_body2eci * F_solar_body;
            a_aero = a_aero + F_solar_eci / satMass;
            M_aero_body = M_aero_body + aeroResults.M_solar;
        end
    catch
        % If aero computation fails, set to zero
        a_aero = [0; 0; 0];
        M_aero_body = [0; 0; 0];
    end
    
    %% Main computation
    % LC is control torques

    %initialization
    Xd = zeros(13,1);
    
    % Perturbation functions
    a_J2 = J2(X);
    
    % 2BP(states 1:6)
    Xd(1:3) = X(4:6);
    
    % Translational dynamics with perturbations
    % Includes: Two-body gravity + J2 + Aerodynamic drag + Solar radiation pressure
    Xd(4:6) = -mu*X(1:3)/r^3 + a_J2 + a_aero;
    
    % quaternion kinematics (states 7:10) 
    B = [X(7) -X(8) -X(9) -X(10); X(8) X(7) -X(10) X(9); X(9) X(10) X(7) -X(8); X(10) -X(9) X(8) X(7)];
    Xd(7:10) = .5*B*[0;X(11);X(12);X(13)];

    % Kinetics(states 11:13)
    % Includes: Control torques + Aerodynamic moments
    LC = simParams.finalParams.Controller(t, X);
    
    % Total external torque (control + aerodynamic)
    L_total = LC + M_aero_body;

    WX = [0 -X(13) X(12); X(13) 0 -X(11); -X(12) X(11) 0];

    Xd(11:13) = ICB \ (L_total - WX * ICB * X(11:13));

    % Add a controlling term, which u = K*X(11:13), where K the gain (need
    % to tune), can also add an integration term to reduce steady state
    % error
    
    % Xdot = AX + Bu, B changes with the type of actuator we use

end
