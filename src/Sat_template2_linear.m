function Xd = Sat_template2_linear(t,X)
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


    %% Main computation
    % LC is initial torques

    %initialization
    Xd = zeros(13,1);
    x = X(1);
    y = X(2);
    z = X(3);
    Vx = X(4);
    Vy = X(5);
    Vz = X(6);
    q0 = X(7);
    q1 = X(8);
    q2 = X(9);
    q3 = X(10);
    wx = X(11);
    wy = X(12);
    wz = X(13);

    % Constants
    mu = 3.986e14;
    r = norm(X(1:3));

    % Moment of inertia tensor
    % Can start with approximating a sphere (think Sputnik)
    I = 2/5*83*(.58/2)^2*[1 .001 .001 ; .001 1 .001 ;.001 .001 1]; % [kg m^2]
    
    % 2BP(states 1:6)
    Xd(1:3) = [1 0 0; 0 1 0; 0 0 1]*X(4:6);
    
    %Creating Jacobian for Velocities
    J11 = (3*mu*x^2)/r^5 - mu/r^3;
    J22 = (3*mu*y^2)/r^5 - mu/r^3;
    J33 = (3*mu*z^2)/r^5 - mu/r^3;
    J4to6 = [J11 0 0; 0 J22 0; 0 0 J33];
    Xd(4:6) = J4to6*X(1:3);
    
    % quaternion kinematics (states 7:10) 
    B = [q0 -q1 -q2 -q3; q1 q0 -q3 q2; q2 q3 q0 -q1; q3 -q2 q1 q0];
    Xd(7:10) = .5*B*[0;wx;wy;wz];

    %calculating u
    w_r = [0;0;0];
    wdot_r = [0;0;0];
    P = [.15 0 0; 0 .15 0; 0 0 .15];
    Kp = [.0025 0 0 ; 0 .0025 0 ; 0 0 .0025];
    delw = X(11:13) - w_r;
    u = -Kp * X(8:10);
    u = u- P * delw;
    u = u + I * wdot_r;
    u = u - cross(X(11:13),w_r);
    % u = u + X(11:13)' * ICB * X(11:13)
    % if abs(norm(u)) >= 1
    %     disp("nonsense")
    %     disp(norm(u))
    % end

    %Creating Jacobian for Angular Velocities
    J11 = (-I(3,1)*wy + I(2,1)*wz)/I(1,1) + (-2*I(2,1)*wx + I(1,1)*wy - I(2,2)*wy -I(2,3)*wz)/I(1,3) + (2*I(3,1)*wx + I(3,2)*wy - I(1,1)*wz + I(3,3)*wz)/I(1,2);
    J12 = (I(3,2)*wx - I(1,2)*wz)/I(1,2) + (I(1,1)*wx - I(2,2)*wx + 2*I(1,2)*wy + I(1,3)*wz)/I(1,3) +(-I(3,1)*wx - 2*I(3,2)*wy +I(2,2)*wz - I(3,3)*wz)/I(1,1);
    J13 = (-I(2,3)*wx + I(1,3)*wy)/I(1,3) + (-I(1,1)*wx + I(3,3)*wx - I(1,2)*wy - 2*I(1,3)*wz)/I(1,2) + (I(2,1)*wx + I(2,2)*wy - I(3,3)*wy + 2*I(2,3)*wz)/I(1,1);
    J21 = (-I(3,1)*wy + I(2,1)*wz)/I(2,1) + (-2*I(2,1)*wx + I(1,1)*wy - I(2,2)*wy -I(2,3)*wz)/I(2,3) + (2*I(3,1)*wx + I(3,2)*wy - I(1,1)*wz + I(3,3)*wz)/I(2,2);
    J22 = (I(3,2)*wx - I(1,2)*wz)/I(2,2) + (I(1,1)*wx - I(2,2)*wx + 2*I(1,2)*wy + I(1,3)* wz)/I(2,3) + (-I(3,1)*wx - 2*I(3,2)*wy + I(2,2)*wz - I(3,3)*wz)/I(2,1);
    J23 = (-I(2,3)*wx + I(1,3)*wy)/I(2,3) + (-I(1,1)*wx + I(3,3)*wx - I(1,2)*wy -2*I(1,3)*wz)/I(2,2) + (I(2,1)*wx + I(2,2)*wy - I(3,3)*wy + 2*I(2,3)*wz)/I(2,1);
    J31 = (-I(3,1)*wy + I(2,1)*wz)/I(3,1) + (-2*I(2,1)*wx + I(1,1)*wy - I(2,2)*wy - I(2,3)*wz)/I(3,3) + (2*I(3,1)*wx + I(3,2)*wy - I(1,1)*wz + I(3,3)*wz)/I(3,2);
    J32 = (I(3,2)*wx - I(1,2)*wz)/I(3,2) + (I(1,1)*wx - I(2,2)*wx + 2*I(1,2)*wy + I(1,2)*wz)/I(3,3) + (-I(3,1)*wx - 2*I(3,2)*wy +I(2,2)*wz - I(3,3)*wz)/I(3,1);
    J33 = (-I(2,3)*wx + I(1,3)*wy)/I(3,3)+ (-I(1,1)*wx + I(3,3)*wx - I(1,2)*wy - 2*I(1,3)*wz)/I(3,2) + (I(2,1)*wx + I(2,2)*wy - I(3,3)*wy + 2*I(2,3)*wz)/I(3,1);

    J11to13X = [ J11 J12 J13;
        J21 J22 J23;
        J31 J32 J22];
    J11to13U = [1/I(1,1) 0 0; 0 1/I(2,2) 0; 0 0 1/I(3,3)];


    Xd(11:13) = J11to13X *X(11:13) + J11to13U*u;


    % Add a controlling term, which u = K*X(11:13), where K the gain (need
    % to tune), can also add an integration term to reduce steady state
    % error
    
    % Xdot = AX + Bu, B changes with the type of actuator we use

end
