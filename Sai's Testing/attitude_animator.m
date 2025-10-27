clc
close all

%% Initializing video writer
save = true;

if save
    v = VideoWriter('myAnimation.mp4', 'MPEG-4');
    v.FrameRate = 30; % Set frame rate to 60 frames per second
    v.Quality = 75;   % Set video quality (0-100)
    open(v);
end

%% Readying plot
vec_fig = figure(1);
set(vec_fig, 'CloseRequestFcn', @myCloseRequestFunction);
ax = axes;

% Defining origin
origin = scatter3([0, 0, 0], [0, 0, 0], [0, 0, 0], 40, 'black', '+');
hold on;

%% Defining inertial reference frame axes
i_1 = [1, 0, 0]'; % Vernal equinox
i_2 = [0, 1, 0]'; % Normal to i_1 and i_3
i_3 = [0, 0, 1]'; % North Pole

inertial_arrows = draw_frame([i_1, i_2, i_3], 'black', ax);
i_1_text = text(i_1(1), i_1(2), i_1(3), 'i_1', 'color', 'k');
i_2_text = text(i_2(1), i_2(2), i_2(3), 'i_2', 'color', 'k');
i_3_text = text(i_3(1), i_3(2), i_3(3), 'i_3', 'color', 'k');

%% Defining orbital reference frame axes
step = 1;

pos = rs(step, :);
vel = vs(step, :);

% o_1 is the radially outward direction
o_1 = pos' / norm(pos);

% o_2 is the direction of velocity
o_2 = vel' / norm(vel);

% o_3 is the direction perpendicular to both o_1 and o_2 (cross product)
o_3 = cross(o_1, o_2);

% orbital_transform = hgtransform(Parent=ax);
orbital_arrows = draw_frame([o_1, o_2, o_3], 'red', ax);

o_1_text = text(o_1(1), o_1(2), o_1(3), 'o_1', 'color', 'r');
o_2_text = text(o_2(1), o_2(2), o_2(3), 'o_2', 'color', 'r');
o_3_text = text(o_3(1), o_3(2), o_3(3), 'o_3', 'color', 'r');

%% Defining body reference frame axes in inertial frame
b_1 = [0, 1, 0]';
b_2 = [1, 0, 0]';
b_3 = [0, 0, -1]';

% body_transform = hgtransform(Parent=ax);
body_arrows = draw_frame([b_1, b_2, b_3], 'blue', ax);

b_1_text = text(b_1(1), b_1(2), b_1(3), 'b_1', 'color', 'b');
b_2_text = text(b_2(1), b_2(2), b_2(3), 'b_2', 'color', 'b');
b_3_text = text(b_3(1), b_3(2), b_3(3), 'b_3', 'color', 'b');

%% Formatting figure
axis equal;
grid on;

view(30, 30)

bounds = 1;
xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])

xlabel('X')
ylabel('Y')
zlabel('Z')

% Update the plot title and legend
title('Inertial and Orbital Reference Frames');
legend([inertial_arrows{1}, orbital_arrows{1}, body_arrows{1}], {'Inertial Frame', 'Orbital Frame', 'Body Frame'});

t_start = 2.5 * 60; % Start 5 minutes into simulation
duration = 5 * 60; % Lasts for 5 minutes
omega_str = "\omega_x = " + compose("%4.3e", rad2deg(omegas(1, 1))) + " deg/s";
test_ann = annotation(vec_fig, 'textbox', [0.15 0.85 0.2 0.05], 'String', {"Torque about i_1:", "0 Nm", omega_str}, FitBoxToText="on");
test_ann.HorizontalAlignment = "center";
test_ann.BackgroundColor = "white";

start = 1;
step_size = 1;

for i = start : step_size : length(rs)
    if isvalid(vec_fig)
        % Create title string with the current frame number
        title_str = ['Inertial and Orbital Reference Frames', newline, 'Time = ', num2str(round(ts(i) / 60, 2)), ' min'];
    
        % Set the title
        title(title_str);

        omega_str = "\omega_x = " + compose("%4.3e", rad2deg(omegas(i, 1))) + " deg/s";

        if (t_start < ts(i)) && (ts(i) <= t_start + duration/2)
            test_ann.String = {"Torque about i_1:", "0.0001 Nm", omega_str, "Accelerating"};
            test_ann.BackgroundColor = "green";
        elseif (t_start + duration/2 < ts(i)) && (ts(i) <= t_start + duration)
            test_ann.String = {"Torque about i_1:", "-0.0001 Nm", omega_str, "Decelerating"};
            test_ann.BackgroundColor = "red";
        else
            test_ann.String = {"Torque about i_1:", "0 Nm", omega_str};
            test_ann.BackgroundColor = "white";
        end
        
        % Extract current position and velocity vectors
        pos = rs(i, :);
        vel = vs(i, :);
        
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
        
        beta_now = betas(i, :)';
        % [e_hat, phi] = PATfromQ(beta_now);

        % body_matrix = makehgtform("axisrotate", e_hat', deg2rad(phi));
        % body_transform.Matrix = body_matrix;

        R_BI = DCMfromQ(beta_now);
        b_1 = R_BI' * [0; 1; 0];
        b_2 = R_BI' * [1; 0; 0];
        b_3 = R_BI' * [0; 0; -1];

        set(body_arrows{1}, 'UData', b_1(1), 'VData', b_1(2), 'WData', b_1(3));
        set(body_arrows{2}, 'UData', b_2(1), 'VData', b_2(2), 'WData', b_2(3));
        set(body_arrows{3}, 'UData', b_3(1), 'VData', b_3(2), 'WData', b_3(3));

        b_1_text.Position = b_1;
        b_2_text.Position = b_2;
        b_3_text.Position = b_3;
    
        drawnow;
        
        if save
            frame = getframe(gcf); % Capture the current figure as a frame
            writeVideo(v, frame);  % Write the captured frame to the video object
        end

    end
end

if save
    close(v);
end

%% Functions
function arrows = draw_frame(vectors, color, parent)
    
    arrows = cell(1, 3);
    for i = 1:3
        arrows{i} = quiver3(0, 0, 0, vectors(1, i), vectors(2, i), vectors(3, i), color, 'LineWidth', 2, Parent=parent);
    end

end

function myCloseRequestFunction(src, ~)
    disp('Figure is about to close!');
    delete(src); 
end