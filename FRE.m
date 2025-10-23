function R = FRE(ax, phi)

cph = cos(phi);
sph = sin(phi);

if ax == 1    
    R = [1 0 0; 0 cph sph; 0 -sph cph];
elseif ax == 2
    R = [cph 0 -sph; 0 1 0; sph 0 cph];
elseif ax == 3
    R = [cph sph 0; -sph cph 0; 0 0 1];
else
    R = eye(3);
end