function RV = RVfromOE(orbit)

%this function calculates the classical two-body orbital parameters 
%using R and V vectors at one instant as input
global mu_sun; global mu_earth; global mu;
global mum;
mu_sun = 132712440000.00002;
mu_earth = 3.986004e14;
erad = 6378.14e3;
dtr = pi/180;
mu = mu_earth;

a = orbit(1);
e = orbit(2);
i = orbit(3);
Om = orbit(4);
om = orbit(5);
phi = orbit(6);

cphi = cos(phi);
sphi = sin(phi);

R313 = FRE(3, om)*FRE(1, i)*FRE(3, Om);  %this is ROI
R313T = R313';      %this is RIO

E = -mu/2/a;
p = a*(1 - e^2);
h = sqrt(mu*p);
r = p/(1 + e*cphi);
v = sqrt(2*(E + mu/r));
phid = h/r^2;
rd = p*e*sphi*sqrt(mu/p^3);


phat = R313T(:,1);
qhat = R313T(:,2);
hhat = R313T(:,3);

rbar = r*(phat*cphi + qhat*sphi);
vbar = (rd*cphi - r*phid*sphi)*phat + (rd*sphi + r*phid*cphi)*qhat;

RV = [rbar vbar];