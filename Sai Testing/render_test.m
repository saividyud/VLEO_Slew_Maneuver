% clear
clc
close all

fig = figure(1);
t = tiledlayout(fig, 2, 2);

fig.WindowState = "maximized";
% fig.Position = [100, 100, 1060, 1060*0.8];

%% First axis
ax1 = nexttile;

E = wgs84Ellipsoid;
[x,y,z] = ellipsoid(0, 0, 0, E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis);
surf(x, y, z, FaceAlpha="texturemap", FaceColor="texturemap", EdgeAlpha="texturemap", Parent=ax1); % Plot the Earth
hold on

pos = scatter3(rs(1, 1), rs(1, 2), rs(1, 3), Parent=ax1);

hold off
axis equal;

view(30, 30)

xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Satellite Trajectory');

bounds = 2 * earthRadius;
xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])

%% Second axis
ax2 = nexttile;

a = 250e3 + earthRadius; % 250 km above Earth semimajor axis
e = 0; % Eccentricity
i = 0; % Inclination
raan = 0; % Right ascension of ascending node
aop = 0; % Argument of periapse
ta = 0; % True anomaly

orbit = [a, e, deg2rad(i), deg2rad(raan), deg2rad(aop), deg2rad(ta)];

RV = RVfromOE(orbit);

r_i = RV(:, 1)'; % [m]
v_i = RV(:, 2)'; % [m/s]

b_1 = v_i' / norm(v_i);
b_3 = -r_i' / norm(r_i);
b_2 = cross(b_3, b_1); % Calculate b_2 as the cross product of b_3 and b_1

R_BI = [
    b_1';
    b_2';
    b_3'
];

beta_i = QfromDCM(R_BI);

data = stlread('./Sai Testing/ADBSat/inou/obj_files/6U CubeSat.stl');
triangle_mat = data.ConnectivityList;
points = 2 * data.Points;

CoM = mean(points, 1);
points = points - CoM;

stl_1 = [0, 0, 1]';
stl_2 = [0, 1, 0]';
stl_3 = [-1, 0, 0]';

R_stlI = [
    stl_1';
    stl_2';
    stl_3'
];

points = R_BI' * R_stlI * points';
% points = R_BI * points';

stl_1 = R_BI' * R_stlI * stl_1;
stl_2 = R_BI' * R_stlI * stl_2;
stl_3 = R_BI' * R_stlI * stl_3;

h = trimesh(triangle_mat, points(1, :), points(2, :), points(3, :), FaceColor="red", FaceAlpha=0.1, EdgeColor='red', EdgeAlpha=0.5, Parent=ax2);
hold on
scatter3(0, 0, 0, "black", Marker="+", Parent=ax2)

bounds = 1;

% o_1 is the radially outward direction
o_1 = r_i' / norm(r_i);

% o_2 is the direction of velocity
o_2 = v_i' / norm(v_i);

% o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
o_3 = cross(o_1, o_2);

orbital_arrows = draw_frame([o_1, o_2, o_3], 'red', ax2);

o_1_text = text(o_1(1), o_1(2), o_1(3), 'o_1', 'color', 'r', Parent=ax2);
o_2_text = text(o_2(1), o_2(2), o_2(3), 'o_2', 'color', 'r', Parent=ax2);
o_3_text = text(o_3(1), o_3(2), o_3(3), 'o_3', 'color', 'r', Parent=ax2);

% b_1 = [0, 0, 1]';
% b_2 = [0, 1, 0]';
% b_3 = [-1, 0, 0]';

body_arrows = draw_frame([b_1, b_2, b_3], 'blue', ax2);

b_1_text = text(b_1(1), b_1(2), b_1(3), 'b_1', 'color', 'b', Parent=ax2);
b_2_text = text(b_2(1), b_2(2), b_2(3), 'b_2', 'color', 'b', Parent=ax2);
b_3_text = text(b_3(1), b_3(2), b_3(3), 'b_3', 'color', 'b', Parent=ax2);

stl_arrows = draw_frame([stl_1, stl_2, stl_3], 'magenta', ax2);

stl_1_text = text(stl_1(1), stl_1(2), stl_1(3), 'stl_1', 'color', 'm', Parent=ax2);
stl_2_text = text(stl_2(1), stl_2(2), stl_2(3), 'stl_2', 'color', 'm', Parent=ax2);
stl_3_text = text(stl_3(1), stl_3(2), stl_3(3), 'stl_3', 'color', 'm', Parent=ax2);

% for i = 1 : 1 : length(data.Points)
%     scatter3(h.Vertices(i, 1), h.Vertices(i, 2), h.Vertices(i, 3))
% end
p = scatter3(h.Vertices(:, 1), h.Vertices(:, 2), h.Vertices(:, 3), Parent=ax2);

axis equal;

xlabel("X")
ylabel("Y");
zlabel("Z");

xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])

title("3D Model of 6U CubeSat");

view(30, 30)
% view(3);
grid on;

step = 1;

% for angle = linspace(0, 45, 500)
%     R = [
%         cosd(angle), sind(angle), 0;
%         -sind(angle), cosd(angle), 0;
%         0, 0, 1
%         ];
%     rotatedPoints = (R * points)';
% 
%     h.Vertices = rotatedPoints; % Update the vertices of the mesh
% 
%     p.XData = h.Vertices(:, 1);
%     p.YData = h.Vertices(:, 2);
%     p.ZData = h.Vertices(:, 3);
% 
%     pos.XData = rs(step, 1);
%     pos.YData = rs(step, 2);
%     pos.ZData = rs(step, 3);
% 
%     drawnow; % Update the figure window
% 
%     step = step + 1;
% end
hold off

%% Third axis
ax3 = nexttile;