function Xd = Sat_template(t,X)
% Sat_template calculates the time rate of change of the state X at a time
% t.
%
% Parameters
% ----------
% t : float
%   Time [s]
% X : 13x1 vector
%   State of the system with the following attributes:
%       - Position (1:3)
%       - Velocity (4:6)
%       - Quaternion (7:10)
%       - Angular velocity (11:13)
%
% Returns
% -------
% Xd : 13x1 vector
%   Dotted state vector

% Xd and X are 13 dimensional state vectors
% X has to be vertical for function to work
% t is time(used for numerical integration)
% LC is initial torques
    %initialization
    Xd = zeros(13,1);
    
    % Perturbation functions
    a_J2 = J2(X);
    
    % Constants
    mu = 3.986e14;
    r = norm(X(1:3));

    % Moment of inertia tensor
    % Can start with approximating a sphere (think Sputnik)
    ICB = 2/5*83*(.58/2)^2*[1 0 0 ; 0 1 0 ;0 0 1]; % [kg m^2]
    
    % 2BP(states 1:6)
    Xd(1:3) = X(4:6);
    
    % extra forces and perturbations can be added here
    Xd(4:6) = -mu*X(1:3)/r^3 + a_J2;
    
    % quaternion kinematics (states 7:10) 
    B = [X(7) -X(8) -X(9) -X(10); X(8) X(7) -X(10) X(9); X(9) X(10) X(7) -X(8); X(10) -X(9) X(8) X(7)];
    Xd(7:10) = .5*B*[0;X(11);X(12);X(13)];

    %calculating u
    w_r = [0;0;0];
    wdot_r = [0;0;0];
    P = [10 0 0; 0 10 0; 0 0 10];
    Kp = [.5 0 0 ; 0 .5 0 ; 0 0 .5];
    delw = X(11:13) - w_r;
    u = -Kp * ([5,5,5]' -X(8:10)) - P * delw + ICB * wdot_r - cross(X(11:13),w_r) + X(11:13)' * ICB * X(11:13);

    
    % Kinetics(states 11:13)
    % can add extra perturbations/Control inputs here 
    %LC = control_torques(t, X);
    LC = u;

    WX = [0 -X(13) X(12); X(13) 0 -X(11); -X(12) X(11) 0];

    Xd(11:13) = inv(ICB)*(LC - WX* ICB * X(11:13));

    % Add a controlling term, which u = K*X(11:13), where K the gain (need
    % to tune), can also add an integration term to reduce steady state
    % error
    
    % Xdot = AX + Bu, B changes with the type of actuator we use

end
