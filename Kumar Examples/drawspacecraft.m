function drawspacecraft(RBI)

dtr = pi/180;
%%%
% draw spacecraft
%%%
orgn = [0 0 0]';
scl = 0.25;
plen = 1;   %panel length
pan = 30;   %panel angle
scol = 'r';
can = cos(pan*dtr);
san = sin(pan*dtr);
pcol = 'c';
sc = [-scl/2, -scl/2, -scl/2;
    -scl/2, scl/2, -scl/2;
    scl/2, scl/2, -scl/2;
    scl/2, -scl/2, -scl/2;
    -scl/2, scl/2, scl/2;
    -scl/2, -scl/2, scl/2;
    scl/2, -scl/2, scl/2;
    scl/2, scl/2, scl/2];
panel1 = [-scl/2, -scl/2, 0;
    -scl/2, scl/2, 0;
    -(scl/2 + plen*can), scl/2, plen*san;
    -(scl/2 + plen*can), -scl/2, plen*san];
panel2 = [scl/2, -scl/2, 0;
    scl/2, scl/2, 0;
    (scl/2 + plen*can), scl/2, plen*san;
    (scl/2 + plen*can), -scl/2, plen*san];
for i = 1:8
    scrot(i,:) = (RBI*(sc(i,:)'))';
end
for i = 1:4
    p1rot(i,:) = (RBI*panel1(i,:)')';
    p2rot(i,:) = (RBI*panel2(i,:)')';
end

fill3(scrot(1:4,1), scrot(1:4,2), scrot(1:4,3), scol); hold on;
fill3(scrot(5:8,1), scrot(5:8,2), scrot(5:8,3), scol);
fill3([scrot(1:2,1); scrot(5:6,1)], [scrot(1:2,2); scrot(5:6,2)], [scrot(1:2,3); scrot(5:6,3)], scol);
fill3([scrot(3:4,1); scrot(7:8,1)], [scrot(3:4,2); scrot(7:8,2)], [scrot(3:4,3); scrot(7:8,3)], scol);
fill3([scrot(3:4,1); scrot(7:8,1)], [scrot(3:4,2); scrot(7:8,2)], [scrot(3:4,3); scrot(7:8,3)], scol);
fill3([scrot(2:3,1); scrot(8:-3:5,1)], [scrot(2:3,2); scrot(8:-3:5,2)], [scrot(2:3,3); scrot(8:-3:5,3)], scol);
fill3([scrot(1:3:4,1); scrot(7:-1:6,1)], [scrot(1:3:4,2); scrot(7:-1:6,2)], [scrot(1:3:4,3); scrot(7:-1:6,3)], scol);
fill3(p1rot(:,1), p1rot(:,2), p1rot(:,3), pcol);
fill3(p2rot(:,1), p2rot(:,2), p2rot(:,3), pcol);
