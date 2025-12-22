clear
clc
close all

fig = figure;
ax = axes;

grp = hgtransform(Parent=ax);

arrow = quiver3(0, 0, 0, 1, 0, 0, LineWidth=1, Parent=grp);

some_text = text(0.05, 0.05, 0.05, 'Hello');

axis equal;

view(30, 30)

bounds = 1;
xlim([-bounds, bounds])
ylim([-bounds, bounds])
zlim([-bounds, bounds])

for ang = linspace(0, 0.5*pi, 1000)
   tm = makehgtform("axisrotate", [0, 1, 0], ang);
   grp.Matrix = tm;

   some_text.Position = [cos(ang) sin(ang) 0.05];


   drawnow
end