clear
clc
close all

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

%% Initializing video writer
save = false;

if save
    v = VideoWriter(fullfile(fileparts(mfilename('fullpath')), '/../assets/Final Report Rough Draft.mp4'), 'MPEG-4');
    v.FrameRate = 30; % Set frame rate to 60 frames per second
    v.Quality = 75;   % Set video quality (0-100)
    open(v);
end

c = vleo.util.constants();

%% Initial state (in Keplerian orbital elements)
a = 250e3 + c.R_earth; % 250 km above Earth semimajor axis
e = 0; % Eccentricity
i = 0; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]

% Computing intial body frame such that b_1 is in the direction of motion,
% b_3 is pointing towards surface (camera pointing), and b_2 is normal to
% both b_3 and b_1
b_1_i = v_i' / norm(v_i);
b_3_i = -r_i' / norm(r_i);
b_2_i = cross(b_3_i, b_1_i); % Calculate b_2 as the cross product of b_3 and b_1

% Computing initial quaternion from these initial body axes
R_BI_i = [
    b_1_i';
    b_2_i';
    b_3_i'
];

beta_i = QfromDCM(R_BI_i)'; % Initial quaternion

% Zero initial angular spin
omega_i = [0, 0, 0]; % Initial angular rate

X_i = [r_i, v_i, beta_i, omega_i, zeros(1,4)]';
Xr = [1; 0; 0; 0; 0; 0; 0]; % Reference: identity quaternion, zero angular rate
%% Simulating
% Simualation bounds
t0 = 0;
t_span = 6*60; % 0.5 hour
dt = 1;

ts = t0 : dt : t_span;

opts = odeset('RelTol', 1e-6,'AbsTol', 1e-6);
[t, X] = ode45(@(t,X) vleo.dynamics.sat_dynamics_linearized_augmented(t,X,Xr) , ts, X_i, opts);

% Extract position and velocity from the state vector
rs = X(:, 1:3);
vs = X(:, 4:6);
betas = X(:, 7:10);
omegas = X(:, 11:13);

% Extract torques
torques = zeros(length(rs), 3);
torques_norm = zeros(length(rs),1);
for i = 1 : 1 : length(rs)
    torques(i, :) = vleo.control.attitude_pd_controller(ts(i), X(i, :)')';
    torques_norm(i) = norm(torques(i,:));
end

%% Reading plot
fig = figure(Name="CubeSat");
set(fig, 'CloseRequestFcn', @myCloseRequestFunction);

% Splitting figure into two plots: left is position with the Earth for
% reference and right is orientation of satellite
t = tiledlayout(fig, 6, 2, TileSpacing="compact", Padding="compact");

% fig.WindowState = "maximized";
set(fig, 'Position', [150, 150, 1280, 720]);

%% Plotting position
ax1 = nexttile(t, 1, [3, 1]);

% Defining Earth using WGS84 ellipsoid
E = wgs84Ellipsoid;
[x,y,z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x, y, z, FaceAlpha="texturemap", FaceColor="texturemap", EdgeAlpha="texturemap", Parent=ax1); % Plot the Earth
hold on;

% Plotting full trajectory of satellite
plot3(rs(:, 1), rs(:, 2), rs(:, 3), LineWidth=2, Color='k', LineStyle='--', Parent=ax1)

% Plotting initial position of satellite
pos_marker = scatter3(rs(1, 1), rs(1, 2), rs(1, 3), 'r', 'filled', Parent=ax1);
hold off;

% Formatting plot
axis equal;

view(ax1, 30, 30)

xlabel(ax1, 'X Position (m)');
ylabel(ax1, 'Y Position (m)');
zlabel(ax1, 'Z Position (m)');
title(ax1, 'Satellite Trajectory');

bounds = 1.1 * c.R_earth;
xlim(ax1, [-bounds, bounds])
ylim(ax1, [-bounds, bounds])
zlim(ax1, [-bounds, bounds])

grid on;

%% Plotting orientation
ax2 = nexttile(t, 7, [3, 1]);

% Defining origin
origin = scatter3([0, 0, 0], [0, 0, 0], [0, 0, 0], 40, 'black', '+', Parent=ax2);
hold on;

% Defining inertial reference frame axes
i_1 = [1, 0, 0]'; % Vernal equinox
i_2 = [0, 1, 0]'; % Normal to i_1 and i_3
i_3 = [0, 0, 1]'; % North Pole

inertial_arrows = draw_frame([i_1, i_2, i_3], 'black', ax2);
i_1_text = text(i_1(1), i_1(2), i_1(3), 'i_1', 'color', 'k', Parent=ax2);
i_2_text = text(i_2(1), i_2(2), i_2(3), 'i_2', 'color', 'k', Parent=ax2);
i_3_text = text(i_3(1), i_3(2), i_3(3), 'i_3', 'color', 'k', Parent=ax2);

% Defining orbital reference frame axes
step = 1;

pos = rs(step, :);
vel = vs(step, :);

% o_1 is the radially outward direction
o_1 = pos' / norm(pos);

% Guessing o_2 is the direction of velocity
guess_o_2 = vel' / norm(vel);

% o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
o_3 = cross(o_1, guess_o_2);

% o_2 is the direction perpendicular to o_1 and o_2 (cross product)
o_2 = cross(o_3, o_1);

orbital_arrows = draw_frame([o_1, o_2, o_3], 'red', ax2);

o_1_text = text(o_1(1), o_1(2), o_1(3), 'o_1', 'color', 'r', Parent=ax2);
o_2_text = text(o_2(1), o_2(2), o_2(3), 'o_2', 'color', 'r', Parent=ax2);
o_3_text = text(o_3(1), o_3(2), o_3(3), 'o_3', 'color', 'r', Parent=ax2);

% Reading in CubeSat model
data = stlread(fullfile(fileparts(mfilename('fullpath')), '../data/geometry/6U CubeSat.STL'));
triangle_mat = data.ConnectivityList;
points = 2 * data.Points;

% Ensuring the origin is at the center of the CubeSat
CoM = mean(points, 1);
points = points - CoM;

% The STL file has weird axes, so adjusting axes such that it lines up with
% the inertial frame
stl_1 = [0, 0, 1]';
stl_2 = [0, 1, 0]';
stl_3 = [-1, 0, 0]';

R_stlI = [
    stl_1';
    stl_2';
    stl_3'
];

points = R_stlI * points';

% Translating CubeSat into body frame
points_i = R_BI_i' * points;

% Defining CubeSat mesh
cubesat = trimesh(triangle_mat, ...
    points_i(1, :), ...
    points_i(2, :), ...
    points_i(3, :), ...
    FaceColor="blue", FaceAlpha=0.1, EdgeColor='blue', EdgeAlpha=0.5, ...
    Parent=ax2);

% Drawing initial body reference frame axes in inertial frame
body_arrows = draw_frame([b_1_i, b_2_i, b_3_i], 'blue', ax2);

b_1_text = text(b_1_i(1), b_1_i(2), b_1_i(3), 'b_1', 'color', 'b', Parent=ax2);
b_2_text = text(b_2_i(1), b_2_i(2), b_2_i(3), 'b_2', 'color', 'b', Parent=ax2);
b_3_text = text(b_3_i(1), b_3_i(2), b_3_i(3), 'b_3', 'color', 'b', Parent=ax2);

% Plotting full trajectory of satellite
traj = plot3(rs(:, 1) - rs(1, 1), rs(:, 2) - rs(1, 2), rs(:, 3) - rs(1, 3), LineWidth=2, Color='k', Parent=ax2, LineStyle='--');

% Formatting figure
axis equal;
grid on;

view(ax2, 30, 30)

bounds = 1;
xlim(ax2, [-bounds, bounds])
ylim(ax2, [-bounds, bounds])
zlim(ax2, [-bounds, bounds])

xlabel(ax2, 'X')
ylabel(ax2, 'Y')
zlabel(ax2, 'Z')

title(ax2, 'Inertial and Orbital Reference Frames');
lgd = legend(ax2, [inertial_arrows{1}, orbital_arrows{1}, body_arrows{1}], {'Inertial Frame', 'Orbital Frame', 'Body Frame'});
lgd.Position = [0.01, 0.35, 0.14, 0.1];

time_str = ['Time = ', num2str(round(ts(1) / 60, 2)), ' min'];
time_ann = annotation(fig, 'textbox', [0.08 0.92 0.1 0.06], 'String', {time_str}, FitBoxToText="off");
time_ann.HorizontalAlignment = "center";
time_ann.VerticalAlignment = "middle";
time_ann.BackgroundColor = "white";

%% Plotting orientation
ax3 = nexttile(t, 2, [2, 1]);

[yaws, pitches, rolls] = quat2angle(betas, "ZYX");
rolls = rad2deg(rolls);
pitches = rad2deg(pitches);
yaws = rad2deg(yaws);

% Plotting the angles over time
plot(ax3, ts/60, rolls, LineWidth=1, Color="r")
hold on
plot(ax3, ts/60, pitches, LineWidth=1, Color="b");
plot(ax3, ts/60, yaws, LineWidth=1, Color="k");

roll_marker = scatter(ax3, ts(1)/60, rolls(1), "r", 'filled');
pitch_marker = scatter(ax3, ts(1)/60, pitches(1), "b", 'filled');
yaw_marker = scatter(ax3, ts(1)/60, yaws(1), "k", 'filled');

ax3line = xline(ax3, ts(1), "k", LineStyle="--", LineWidth=1);
hold off

xlabel(ax3, "Time [min]")
ylabel(ax3, "Angle [deg]")

grid on

range = abs(max([rolls; pitches; yaws]) - min([rolls; pitches; yaws]));
if range == 0
    range = 0.01;
end
ylim(ax3, [min([rolls; pitches; yaws]) - 0.1*range, max([rolls; pitches; yaws]) + 0.1*range])

legend(ax3, [roll_marker, pitch_marker, yaw_marker], {"Roll \phi", "Pitch \theta", "Yaw \psi"}, 'NumColumns', 3, Location="best")

%% Plotting angular velocity
ax4 = nexttile(t, 6, [2, 1]);

omegas = rad2deg(omegas);

% Plotting the angular rate over time
plot(ax4, ts/60, omegas(:, 1), LineWidth=1, Color="r")
hold on
plot(ax4, ts/60, omegas(:, 2), LineWidth=1, Color="b");
plot(ax4, ts/60, omegas(:, 3), LineWidth=1, Color="k");

roll_rate_marker = scatter(ax4, ts(1)/60, omegas(1, 1), "r", 'filled');
pitch_rate_marker = scatter(ax4, ts(1)/60, omegas(1, 2), "b", 'filled');
yaw_rate_marker = scatter(ax4, ts(1)/60, omegas(1, 3), "k", 'filled');

ax4line = xline(ax4, ts(1), "k", LineStyle="--", LineWidth=1);
hold off

xlabel(ax4, "Time [min]")
ylabel(ax4, "Rotation Rate [deg/s]")

grid on

range = abs(max(omegas(:)) - min(omegas(:)));
if range == 0
    range = 0.01;
end
ylim(ax4, [min(omegas(:)) - 0.1*range, max(omegas(:)) + 0.1*range])

legend(ax4, [roll_rate_marker, pitch_rate_marker, yaw_rate_marker], {"\omega_1", "\omega_2", "\omega_3"}, 'NumColumns', 3, Location="best")

%% Plotting torques
ax5 = nexttile(t, 10, [2, 1]);

% Plotting the angular rate over time
plot(ax5, ts/60, torques(:, 1), LineWidth=1, Color="r")
hold on
plot(ax5, ts/60, torques(:, 2), LineWidth=1, Color="b");
plot(ax5, ts/60, torques(:, 3), LineWidth=1, Color="k");

tau_1_marker = scatter(ax5, ts(1)/60, torques(1, 1), "r", 'filled');
tau_2_marker = scatter(ax5, ts(1)/60, torques(1, 2), "b", 'filled');
tau_3_marker = scatter(ax5, ts(1)/60, torques(1, 3), "k", 'filled');

ax5line = xline(ax5, ts(1), "k", LineStyle="--", LineWidth=1);
hold off

xlabel(ax5, "Time [min]")
ylabel(ax5, "Torque [Nm]")

grid on

range = abs(max(torques(:)) - min(torques(:)));
if range == 0
    range = 0.01;
end
ylim(ax5, [min(torques(:)) - 0.1*range, max(torques(:)) + 0.1*range])

legend(ax5, [tau_1_marker, tau_2_marker, tau_3_marker], {"\tau_1", "\tau_2", "\tau_3"}, 'NumColumns', 3, Location="best")

set(findall(fig, '-property', 'FontSize'), 'FontSize', 14);

%% Animating
start = 1;
step_size = 5;

for i = start : step_size : length(rs)
    if isvalid(fig)
        % Updating time
        time_str = ['Time = ', num2str(round(ts(i)/60, 2)), ' min'];
        time_ann.String = {time_str};
        
        % Extract current position and velocity vectors
        pos = rs(i, :);
        vel = vs(i, :);

        % Update position of CubeSat
        pos_marker.XData = pos(1);
        pos_marker.YData = pos(2);
        pos_marker.ZData = pos(3);
        
        % o_1 is the radially outward direction
        o_1 = pos' / norm(pos);
        
        % o_2 is the direction of velocity
        o_2 = vel' / norm(vel);
        
        % o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
        o_3 = cross(o_1, o_2);
        
        % Changing orbital frame vectors
        set(orbital_arrows{1}, 'UData', o_1(1), 'VData', o_1(2), 'WData', o_1(3));
        set(orbital_arrows{2}, 'UData', o_2(1), 'VData', o_2(2), 'WData', o_2(3));
        set(orbital_arrows{3}, 'UData', o_3(1), 'VData', o_3(2), 'WData', o_3(3));

        o_1_text.Position = o_1;
        o_2_text.Position = o_2;
        o_3_text.Position = o_3;
        
        % Rotating body frame
        beta_now = betas(i, :)';

        R_BI = DCMfromQ(beta_now);
        b_1 = R_BI(1, :)';
        b_2 = R_BI(2, :)';
        b_3 = R_BI(3, :)';

        set(body_arrows{1}, 'UData', b_1(1), 'VData', b_1(2), 'WData', b_1(3));
        set(body_arrows{2}, 'UData', b_2(1), 'VData', b_2(2), 'WData', b_2(3));
        set(body_arrows{3}, 'UData', b_3(1), 'VData', b_3(2), 'WData', b_3(3));

        b_1_text.Position = b_1;
        b_2_text.Position = b_2;
        b_3_text.Position = b_3;
        
        % Rotating CubeSat model
        cubesat.Vertices = (R_BI' * points)'; % Update the vertices of the mesh

        % Following satellite in first axes
        bounds = 3000e3;
        xlim(ax1, [pos(1) - bounds, pos(1) + bounds]);
        ylim(ax1, [pos(2) - bounds, pos(2) + bounds]);
        zlim(ax1, [pos(3) - bounds, pos(3) + bounds]);

        % Updating trajectory in orientation figure
        traj.XData = rs(:, 1) - pos(1);
        traj.YData = rs(:, 2) - pos(2);
        traj.ZData = rs(:, 3) - pos(3);

        % Updating third axis
        roll_marker.XData = ts(i)/60;
        roll_marker.YData = rolls(i);

        pitch_marker.XData = ts(i)/60;
        pitch_marker.YData = pitches(i);

        yaw_marker.XData = ts(i)/60;
        yaw_marker.YData = yaws(i);

        ax3line.Value = ts(i)/60;

        % Updating fourth axis
        roll_rate_marker.XData = ts(i)/60;
        roll_rate_marker.YData = omegas(i, 1);

        pitch_rate_marker.XData = ts(i)/60;
        pitch_rate_marker.YData = omegas(i, 2);

        yaw_rate_marker.XData = ts(i)/60;
        yaw_rate_marker.YData = omegas(i, 3);

        ax4line.Value = ts(i)/60;

        % Updating fifth axis
        tau_1_marker.XData = ts(i)/60;
        tau_1_marker.YData = torques(i, 1);

        tau_2_marker.XData = ts(i)/60;
        tau_2_marker.YData = torques(i, 2);

        tau_3_marker.XData = ts(i)/60;
        tau_3_marker.YData = torques(i, 3);

        ax5line.Value = ts(i)/60;
    
        drawnow;
        
        if save
            % 2. Use 'fig' specific handle, not 'gcf'
            frame = getframe(fig); 
            
            if isempty(frame.cdata)
                return; % Skip if frame is empty (window minimized/closed)
            end

            % 3. Initialize target size on the very first frame
            if i == start 
                [h, w, ~] = size(frame.cdata);
                
                % Force dimensions to be even for MPEG-4
                h = h - mod(h, 2);
                w = w - mod(w, 2);
                
                % Store this "perfect" size for all future frames
                targetSize = [h, w];
                
                % Crop the first frame
                frame.cdata = frame.cdata(1:h, 1:w, :);
            else
                % 4. FORCE subsequent frames to match the first frame's size exactly
                % This prevents "Frame must be X by Y" errors if the window shifts by 1px
                frame.cdata = imresize(frame.cdata, targetSize);
            end

            writeVideo(v, frame);  % Write the captured frame to the video object
        end

    end
end

if save
    close(v);
end

% exportgraphics(fig, './assets/Final Report Rough Draft 2.png', 'Resolution', 300)

%% Functions
function myCloseRequestFunction(src, ~)
    disp('Figure is closing!');
    delete(src); 
end

function arrows = draw_frame(vectors, color, parent)
    
    arrows = cell(1, 3);
    for i = 1:3
        arrows{i} = quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), color, 'LineWidth', 2, Parent=parent);
    end

end

function RV = RVfromOE(orbit, mu)
    arguments
        orbit
        mu double = 3.986004e14
    end
    
    %this function calculates the classical two-body orbital parameters 
    %using R and V vectors at one instant as input
    % global mu_sun; global mu_earth; global mu;
    % mu_sun = 132712440000.00002;
    % mu_earth = 3.986004e14;
    % erad = 6378.14e3;
    % dtr = pi/180;
    % mu = mu_earth;
    
    a = orbit(1);
    e = orbit(2);
    i = orbit(3);
    Om = orbit(4);
    om = orbit(5);
    phi = orbit(6);
    
    cphi = cos(phi);
    sphi = sin(phi);
    
    R313 = FRE(3, om)*FRE(1, i)*FRE(3, Om);  %this is ROI
    R313T = R313';      %this is RIO
    
    E = -mu/2/a;
    p = a*(1 - e^2);
    h = sqrt(mu*p);
    r = p/(1 + e*cphi);
    v = sqrt(2*(E + mu/r));
    phid = h/r^2;
    rd = p*e*sphi*sqrt(mu/p^3);
    
    
    phat = R313T(:,1);
    qhat = R313T(:,2);
    hhat = R313T(:,3);
    
    rbar = r*(phat*cphi + qhat*sphi);
    vbar = (rd*cphi - r*phid*sphi)*phat + (rd*sphi + r*phid*cphi)*qhat;
    
    RV = [rbar vbar];
end

function beta = QfromDCM(R_BI)
    % QfromDCM calculates the quaternion given a DCM from inertial to body 
    % frame.
    %
    % Parameters
    % ----------
    % R_BI : 3x3 matrix
    %   Rotation DCM for converting from inertial frame to body frame
    %
    % Returns
    % -------
    % beta : 4x1 quaternion
    %   Output quaternion
    
    tr = trace(R_BI);
    
    
    %%%first compute squares of each quaternion
    bb(1) = 0.25*(1 + tr);
    bb(2) = 0.25*(1 + 2*R_BI(1,1) - tr);
    bb(3) = 0.25*(1 + 2*R_BI(2,2) - tr);
    bb(4) = 0.25*(1 + 2*R_BI(3,3) - tr);
    
    %%%find the max!
    [~, mxpos] = max(bb);
    
    %%%identify the quaternion with max sq-value
    beta(mxpos) = sqrt(bb(mxpos));
    
    stan(1) = (R_BI(2,3) - R_BI(3,2))/4;
    stan(2) = (R_BI(3,1) - R_BI(1,3))/4;
    stan(3) = (R_BI(1,2) - R_BI(2,1))/4;
    stan(4) = (R_BI(1,2) + R_BI(2,1))/4;
    stan(5) = (R_BI(1,3) + R_BI(3,1))/4;
    stan(6) = (R_BI(2,3) + R_BI(3,2))/4;
    
    if mxpos == 1
        sindex = [1;2;3];
        betaindex = [2;3;4];
    elseif mxpos == 2
        sindex = [1;4;5];
        betaindex = [1;3;4];
    elseif mxpos == 3
        sindex = [2;4;6];
        betaindex = [1;2;4];
    else
        sindex = [3;5;6];
        betaindex = [1;2;3];
    end
    
    beta(betaindex) = stan(sindex)/beta(mxpos);
    if beta(1) < 0
        beta = -beta;
    end
    beta = beta';
end

function R_BI = DCMfromQ(beta)
    % DCMfromQ calculates the DCM from inertial to body frame given a 
    % quaternion.
    %
    % Parameters
    % ----------
    % beta : 4x1 quaternion
    %   Input quaternion
    %
    % Returns
    % -------
    % R_BI : 3x3 matrix
    %   Rotation DCM for converting from inertial frame to body frame
        
        % Extracting quaternions
        b0 = beta(1);
        b1 = beta(2);
        b2 = beta(3);
        b3 = beta(4);
    
        % Calculating DCM
        R_BI = [
            b0^2 + b1^2 - b2^2 - b3^2, 2*(b1*b2 + b0*b3), 2*(b1*b3 - b0*b2);
            2*(b1*b2 - b0*b3), b0^2 - b1^2 + b2^2 - b3^2, 2*(b2*b3 + b0*b1);
            2*(b1*b3 + b0*b2), 2*(b2*b3 - b0*b1), b0^2 - b1^2 - b2^2 + b3^2
        ];
    
end
