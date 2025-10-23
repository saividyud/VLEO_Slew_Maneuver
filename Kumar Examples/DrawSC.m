clear all; clc; close all;

global dtr;

dtr = pi/180;

%%%
% draw spacecraft
%%%
orgn = [0 0 0]';
scl = 0.25;
plen = 1;   %panel length
pan = 30;   %panel angle
sc = [-scl/2, -scl/2, -scl/2;
    -scl/2, scl/2, -scl/2;
    scl/2, scl/2, -scl/2;
    scl/2, -scl/2, -scl/2;
    -scl/2, scl/2, scl/2;
    -scl/2, -scl/2, scl/2;
    scl/2, -scl/2, scl/2;
    scl/2, scl/2, scl/2];
scol = 'r';
fill3(sc(1:4,1), sc(1:4,2), sc(1:4,3), scol); hold on;
fill3(sc(5:8,1), sc(5:8,2), sc(5:8,3), scol);
fill3([sc(1:2,1); sc(5:6,1)], [sc(1:2,2); sc(5:6,2)], [sc(1:2,3); sc(5:6,3)], scol);
fill3([sc(3:4,1); sc(7:8,1)], [sc(3:4,2); sc(7:8,2)], [sc(3:4,3); sc(7:8,3)], scol);
fill3([sc(3:4,1); sc(7:8,1)], [sc(3:4,2); sc(7:8,2)], [sc(3:4,3); sc(7:8,3)], scol);
fill3([sc(2:3,1); sc(8:-3:5,1)], [sc(2:3,2); sc(8:-3:5,2)], [sc(2:3,3); sc(8:-3:5,3)], scol);
fill3([sc(1:3:4,1); sc(7:-1:6,1)], [sc(1:3:4,2); sc(7:-1:6,2)], [sc(1:3:4,3); sc(7:-1:6,3)], scol);
can = cos(pan*dtr);
san = sin(pan*dtr);
pcol = 'k';
panel1 = [-scl/2, -scl/2, 0;
    -scl/2, scl/2, 0;
    -(scl/2 + plen*can), scl/2, plen*san;
    -(scl/2 + plen*can), -scl/2, plen*san];
panel2 = [scl/2, -scl/2, 0;
    scl/2, scl/2, 0;
    (scl/2 + plen*can), scl/2, plen*san;
    (scl/2 + plen*can), -scl/2, plen*san];
fill3(panel1(:,1), panel1(:,2), panel1(:,3), pcol);
fill3(panel2(:,1), panel2(:,2), panel2(:,3), pcol);
