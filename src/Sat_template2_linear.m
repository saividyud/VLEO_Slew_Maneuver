function Xd = Sat_template2_linear(t,X,Xr)
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

    %initialize derivative vector
    Xd = zeros(13,1);

    %reference conditions
    q0 = Xr(1);
    q1 = Xr(2);
    q2 = Xr(3);
    q3 = Xr(4);
    wx = Xr(5);
    wy = Xr(6);
    wz = Xr(7);

    % Constants
    mu = 3.986e14;
    r = norm(X(1:3));

    % Moment of inertia tensor
    % Can start with approximating a sphere (think Sputnik)
    I = 2/5*83*(.58/2)^2*[1 .001 .001 ; .001 1 .001 ;.001 .001 1]; % [kg m^2]
    
    % 2BP(states 1:6)
    Xd(1:3) = [1 0 0; 0 1 0; 0 0 1]*X(4:6);

    % extra forces and perturbations can be added here
    Xd(4:6) = -mu*X(1:3)/r^3;

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

    %Creating full jacobian for rotational dynamics
    A = [ 0 -wx/2 -wy/2 wz/2 -q1/2 -q2/2 q3/2;
        wx/2 0 wz/2 -wy/2 q0/2 -q3/2 q2/2;
        wy/2 -wz/2 0 wx/2 q3/2 q0/2 -q1/2;
        wz/2 wy/2 -wx/2 0 -q2/2 q1/2 q0/2;
        0 0 0 0 J11 J12 J13;
        0 0 0 0 J21 J22 J23;
        0 0 0 0 J31 J32 J33];

    %Proportional Term
    Kp = [.0025 0 0 ; 0 .0025 0 ; 0 0 .0025];
    KpOm = [.1 0 0 ; 0 .1 0; 0 0 .1];
    u = -Kp * X(8:10);
    u = u - KpOm*X(11:13);
    
    %integral term
    % Ki = .001;
    % u = u - Ki * [sum(xHist(8,:));sum(xHist(9,:));sum(xHist(10,:))];
    
    %Derivative Term
    Kd = .1;
    Xd2 = A*X(7:13);
    u = u - Kd*Xd2(5:7);
    
    B = [ zeros(4,3); inv(I)];

    Xd(7:13) = A*X(7:13) + B*u;

    % Add a controlling term, which u = K*X(11:13), where K the gain (need
    % to tune), can also add an integration term to reduce steady state
    % error
    
    % Xdot = AX + Bu, B changes with the type of actuator we use

end
