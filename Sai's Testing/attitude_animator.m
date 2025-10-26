clc
close all

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

orbital_transform = hgtransform(Parent=ax);
orbital_arrows = draw_frame([o_1, o_2, o_3], 'red', orbital_transform);

%% Defining body reference frame axes
b_1 = [1, 0, 0]';
b_2 = [0, 1, 0]';
b_3 = [0, 0, 1]';

body_transform = hgtransform(Parent=ax);
body_arrows = draw_frame([b_1, b_2, b_3], 'blue', body_transform);

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

step_size = 100;

for i = 1 : step_size : length(rs)
    if isvalid(vec_fig)
        % Create title string with the current frame number
        title_str = ['Inertial and Orbital Reference Frames', newline, 'Time = ', num2str(round(ts(i) / 3600, 2)), ' hrs'];
    
        % Set the title
        title(title_str);
        
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
        
        beta_now = betas(i, :)';
        [e_hat, phi] = PATfromQ(beta_now);

        body_matrix = makehgtform("axisrotate", e_hat', phi);
        body_transform.Matrix = body_matrix;
    
        drawnow;
    end
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