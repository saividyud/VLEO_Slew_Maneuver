function TBorbit = OEfromRV(rbar, vbar)

%this function calculates the classical two-body orbital parameters 
%using R and V vectors at one instant as input
global mu_sun; global mu_earth;
mu_sun = 132712440000.00002;
mu_earth = 3.986004e14;
mu = mu_earth;

% Ensure column vectors
rbar = rbar(:);
vbar = vbar(:);

R = norm(rbar);
V = norm(vbar);
a = 1/(2/R - V^2/mu);  %SEMI-MAJOR AXIS
hbar = cross(rbar, vbar);   %ANGULAR MOMENTTUM VECTOR
h = norm(hbar);
p = h^2/mu; %PARAMETER
e = real(sqrt(1 - p/a));  %ECCENTRICITY
En = -mu/2/a;   %Energy
sig = dot(rbar, vbar)/sqrt(mu); %QUANTITY 'SIGMA'
cbar = cross(vbar, hbar) - mu*rbar/R;
ihhat = hbar/h;

% Handle small eccentricity
if e < 1e-9
    e = 0;
    iehat = [1; 0; 0]; % Arbitrary reference for circular orbit (Column vector)
else
    iehat = cbar/(mu*e);
end

imhat = cross(ihhat, iehat);
imhat = imhat/norm(imhat);
i = acos(ihhat(3)); %INCLINATION;
si = sin(i);

% Handle small inclination
if abs(si) < 1e-9
    OMEGA = 0;
else
    OMEGA = atan2(ihhat(1)/si, -ihhat(2)/si);
end

if OMEGA < 0
    OMEGA = OMEGA + 2*pi;
end

% Handle undefined Argument of Periapsis for circular/equatorial
if e < 1e-9 || abs(si) < 1e-9
    omega = 0;
else
    omega = atan2(iehat(3)/si, imhat(3)/si);   %ARGUMENT OF PERIAPSIS
end

%%%ROI
ROI = [iehat'; imhat'; ihhat'];

%%%compute true anomaly
brO = ROI*rbar;
cosf = brO(1)/R;
sinf = brO(2)/R;
if abs(brO(3)) > 1e-9
    fprintf('Please check for out of plane error: %1.2e \n', brO(3));
end

% E2 = atan2((sig/sqrt(a)), (1 - R/a));  %ECCENTRIC ANOMALY
% if E2 < 0
%     E2 = E2 + 2*pi;
% end
%Alternate way to find eccentric anomaly
% fpa = asin(dot(rbar, vbar)/(R*V));
% sfpa = sin(fpa);
% rdot = V*sfpa;
% sinf = rdot*p/(h*e);
% cosf = -(1 - p/R)/e;

f = atan2(sinf, cosf);
if f < 0
    f = f + 2*pi;
end
fac = sqrt((1-e)/(1+e));
E = 2*atan(fac*tan(f/2));
if E < 0
    E = E + 2*pi;
end
%M = E - e*sin(E);
TBorbit = [a; e; i; OMEGA; omega; f];