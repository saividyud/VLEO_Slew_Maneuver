clear
clc
close all

test_cell = cell(1, 3);

test_cell{1} = 1;
test_cell{2} = 'Hello';
test_cell{3} = rand(5); % Assign a random 5x5 matrix

collection_of_cells = cell(2, 3);

collection_of_cells(1, :) = test_cell;

figure;

arrow = quiver3(0, 0, 0, 1, 0, 0, LineWidth=1);
arrow.UData
arrow.VData
arrow.WData

axis equal;

view(30, 30)

bounds = 1;
xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])