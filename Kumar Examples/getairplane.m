function [af, wi, tl, vt] = getairplane(BCM, aclen, acwid, achgt, wspan, wchord, tchord, tspan)

[xa, ya, za] = ellipsoid(BCM(1), BCM(2), BCM(3), aclen/2, acwid/2, achgt/2, 50);    %fuselage
[xw, yw, zw] = ellipsoid(BCM(1) + aclen/20, BCM(2), BCM(3), wchord/2, wspan/2, wchord/10, 50);   %wing
[xt, yt, zt] = ellipsoid(BCM(1) - aclen/2.5, BCM(2), BCM(3), tchord/2, tspan/2, tchord/5, 50);  %tail
[xv, yv, zv] = ellipsoid(BCM(1) - aclen/2.5, BCM(2), BCM(3), tchord/2, tchord/5, tspan/2);  %vertical fin
fz = find(zv < BCM(3));
xv(fz) = BCM(1); yv(fz) =  BCM(2); zv(fz) = BCM(3);

af.x = xa; af.y = ya; af.z = za;
wi.x = xw; wi.y = yw; wi.z = zw;
tl.x = xt; tl.y = yt; tl.z = zt;
vt.x = xv; vt.y = yv; vt.z = zv;
