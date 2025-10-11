%**************************************************************************
%     function [d,t] = msis86(day,sec,alt,glat,glong,stl,f107a,f107,ap)
%     implicit none
%       msis-86/cira 1986 neutral thermosphere model
%        a.e.hedin 3/15/85;2/26/87 (variable names shortened)
%        10/14/87 increase altitude limit of o mixing calculation
%            altl(2) from 300.0 to 400.0 km .
%    inputs:
%	    day - day number of the year.
%       sec - ut(sec)
%       alt - altitude(km) (greater than 85 km)
%       glat - geodetic latitude(deg)
%       glong - geodetic longitude(deg)
%       stl - local apparent solar time(hrs)
%       f107a - 3 month average of f10.7 flux
%       f107 - daily f10.7 flux for previous day
%       ap - magnetic index(daily) or when sw(9)=-1. :
%          - array containing:
%            (1) daily ap
%            (2) 3 hr ap index for current time
%            (3) 3 hr ap index for 3 hrs before current time
%            (4) 3 hr ap index for 6 hrs before current time
%            (5) 3 hr ap index for 9 hrs before current time
%            (6) average of eight 3 hr ap indicies from 12 to 33 hrs prior
%                   to current time
%            (7) average of eight 3 hr ap indicies from 36 to 57 hrs prior
%                   to current time
%    outputs: 
%       d(1) - He number density(cm-3)
%       d(2) - O number density(cm-3)
%       d(3) - N2 number density(cm-3)
%       d(4) - O2 number density(cm-3)
%       d(5) - Ar number density(cm-3)
%       d(6) - total mass density(gm/cm3)
%       d(7) - H number density(cm-3)
%       d(8) - N number density(cm-3)
%       t(1) - exospheric temperature
%       t(2) - temperature at alt
%
%     to get output in m-3 and kg/m3:   call meters(.true.) 
%
%         additional comments
%          (1) lower bound quantities in common/gts3c/
%          (2) to turn on and off particular variations call tselec(sw)
%              where sw is a 25 element array containing 0. for off, 1. 
%              for on, or 2. for main effects off but cross terms on
%              for the following variations
%              1 - f10.7 effect on mean  2 - time independent
%              3 - symmetrical annual    4 - symmetrical semiannual
%              5 - asymmetrical annual   6 - asymmetrical semiannual
%              7 - diurnal               8 - semidiurnal
%              9 - daily ap             10 - all ut/long effects
%             11 - longitudinal         12 - ut and mixed ut/long
%             13 - mixed ap/ut/long     14 - terdiurnal
%             15 - departures from diffusive equilibrium
%             16 - all tinf var         17 - all tlb var
%             18 - all t0 var           19 - all s var
%             20 - all z0 var           21 - all nlb var
%             22 - all tr12 var         23 - turbo scale height var
%
% Last modified:   2022/11/19   Meysam Mahooti
%
%**************************************************************************
function [d, t] = msis86(Mjd_UTC, r_ECEF, gast)

global swdata
global gsurf re sav sw swc pt pd ps pdl_0 pdl_1 ptm pdm
msis86_coeff

[year, month, day, hour, minute, sec] = invjday(Mjd_UTC);
doy = finddays(year, month, day, hour, minute, sec);

[glong, glat, alt] = Geodetic(r_ECEF);
alt = alt/1000;
stl = glong + gast;
stl = mod(stl,2*pi);
stl = (stl*24)/(2*pi); % hours

i = find((year==swdata(1,:)) & (month==swdata(2,:)) & (day==swdata(3,:)),1,'first');
swd_1 = swdata(:,i);
ap(1) = swd_1(23); % Arithmetic average of the 8 AP indices for the day
ap(2) = swd_1(15); % 3 hr AP index for current time
swd_2 = swdata(:,i-1);
ap(3) = swd_2(22); % 3 hr AP index for 3 hrs before current time
ap(4) = swd_2(21); % 3 hr AP index for 6 hrs before current time
ap(5) = swd_2(20); % 3 hr AP index for 9 hrs before current time
sum  = swd_2(19)+swd_2(18)+swd_2(17)+swd_2(16)+swd_2(15);
swd_3 = swdata(:,i-2);
sum = sum+swd_3(22)+swd_3(21)+swd_3(20);
ap(6) = sum/8; % Average of eight 3 hr AP indicies from 12 to 33 hrs prior to current time
swd_4 = swdata(:,i-3);
sum = swd_3(19)+swd_3(18)+swd_3(17)+swd_3(16)+swd_3(15)+swd_4(22)+swd_4(21)+swd_4(20);
ap(7) = sum/8; % Average of eight 3 hr AP indicies from 36 to 57 hrs prior to current time 

f107 = swd_2(30);  % Daily F10.7 flux for previous day (observed).
f107a = swd_1(28); % Centered 81-day arithmetic average of F10.7 (adjusted).
 
d=zeros(8,1);
t=zeros(2,1);
altl = [200.0, 400.0, 150.0, 200.0, 240.0, 450.0, 320.0, 450.0];

% eq. a7
% c!!old!! tinf=ptm(1)*(1.+sw(16)*globe5(yrd,sec,glat,glong,stl,f107a,f107,
% c!!old!!$ ap,pt))*pt(1)
gggg = globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pt);

tinf	= ptm(1)*(1.0 + sw(16)*gggg)*pt(1);
za		= ptm(5)*pdl_1(16);
% eq. a9
t0		= ptm(3)*pd(3,76)*(1.0 + sw(18)*glob5l(pd(3,76:100)));
% eq. a8
tlb		= ptm(2)*(1.0 + sw(17)*glob5l(pd(3,26:50)))*pd(3,26);
% eq. a10
z0		= ptm(7)*(1.0 + sw(20)*glob5l(pd(3,51:75)))*pd(3,51);
% eq. a6
g0		= ptm(4)*ps(1)*(1.0 + sw(19)*...
		globe5(doy, sec, glat, glong, stl, f107a, f107, ap, ps));
% eq. a5
s		= g0/(tinf - tlb);
% eq. a11
tr12	= pd(3,101)*(1.0 + sw(22)*glob5l(pd(3,101:125)));
t(1)	= tinf;
% eq. a18  n2
g28		= sw(21)*glob5l(pd(3,1:25));
t(1)	= tinf;
xmm		= pdm(5,3);

% **** N2 density ****

% eq. a18
db28	= pdm(1,3)*exp(g28)*pd(3,1);
% eq. a13 - a17
[tz, d(3)] = denss(alt, db28, tinf, tlb, 28.0, 0.0, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
% eq. a19
zh28	= pdm(3,3);
zhm28	= pdm(4,3)*pdl_1(6);
xmd		= 28.0 - xmm;
[~, b28] = denss(zh28, db28, tinf, tlb, xmd, -1.0, ptm(6), s, t0, za, z0, tr12);
if ((alt < altl(3)) && (sw(15) ~= 0.0))
	[~, dm28] = denss(alt, b28, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);	
	% eq. a12
	d(3) = dnet(d(3), dm28, zhm28, xmm, 28.0);
end

% **** He density ****

% eq. a18
g4		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(1,:));
db04	= pdm(1,1)*exp(g4)*pd(1,1);
% eq. a13 - a17
[tz, d(1)] = denss(alt, db04, tinf, tlb, 4.0, -.40, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
if ((alt < altl(1)) && (sw(15) ~= 0.0))
	% eq. a19
	zh04	= pdm(3,1);
	[tz, b04] = denss(zh04, db04, tinf, tlb, 4.0 - xmm, -1.40, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm04] = denss(alt, b04, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm04	= zhm28;
	d(1)	= dnet(d(1), dm04, zhm04, xmm, 4.0);
	% eq. a20b
	rl		= log(b28*pdm(2,1)/b04);
	% eq. a20a
	zc04	= pdm(5,1)*pdl_1(1);
	hc04	= pdm(6,2)*pdl_1(2);
	d(1)	= d(1)*ccor(alt, rl, hc04, zc04);
end

% **** O density ****

% eq. a18
g16		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(2,:));
db16	= pdm(1,2)*exp(g16)*pd(2,1);
% eq. a13 - a17
[tz, d(2)] = denss(alt, db16, tinf, tlb, 16.0, 0.0, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
if ((alt <= altl(2)) && (sw(15) ~= 0.0))
	% corrected from pdm(3,1) to pdm(3,2)  12/2/85
	% eq. a19
	zh16	= pdm(3,2);
	[tz, b16] = denss(zh16, db16, tinf, tlb, 16.0 - xmm, -1.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm16] = denss(alt, b16, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm16	= zhm28;
	d(2)	= dnet(d(2), dm16, zhm16, xmm, 16.0);
	% eq. a20b
	rl		= log(b28*pdm(2,2)*abs(pdl_1(17))/b16);
	% eq. a20a
	hc16	= pdm(6,2)*pdl_1(4);
	zc16	= pdm(5,2)*pdl_1(3);
	d(2)	= d(2)*ccor(alt, rl, hc16, zc16);
	% eq. a21
	hcc16	= pdm(8,2)*pdl_1(14);
	zcc16	= pdm(7,2)*pdl_1(13);
	rc16	= pdm(4,2)*pdl_1(15);
	d(2)	= d(2)*ccor(alt, rc16, hcc16, zcc16);
end

% **** O2 density ****

% eq. a18
g32		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(4,:));
db32	= pdm(1,4)*exp(g32)*pd(4,1);
% eq. a13 - a17
[tz, d(4)] = denss(alt, db32, tinf, tlb, 32.0, 0.0, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;

if ((alt <= altl(4)) && (sw(15) ~= 0.0))
	% eq. a19
	zh32	= pdm(3,4);
	[tz, b32] = denss(zh32, db32, tinf, tlb, 32.0 - xmm, -1.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm32] = denss(alt, b32, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm32	= zhm28;
	d(4)	= dnet(d(4), dm32, zhm32, xmm, 32.0);
	% eq. a20b
	rl		= log(b28*pdm(2,4)/b32);
	% eq. a20a
	hc32	= pdm(6,4)*pdl_1(8);
	zc32	= pdm(5,4)*pdl_1(7);
	d(4)	= d(4)*ccor(alt, rl, hc32, zc32);
end

% **** Ar density ****

% eq. a18
g40		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(5,:));
db40	= pdm(1,5)*exp(g40)*pd(5,1);
% eq. a13 - a17
[tz, d(5)] = denss(alt, db40, tinf, tlb, 40.0, 0.0, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
if ((alt <= altl(5)) && (sw(15) ~= 0.0))
	% eq. a19
	zh40	= pdm(3,5);
	[tz, b40] = denss(zh40, db40, tinf, tlb, 40.0 - xmm, -1.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm40] = denss(alt, b40, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm40	= zhm28;
	d(5)	= dnet(d(5), dm40, zhm40, xmm, 40.0);
	% eq. a20b
	rl		= log(b28*pdm(2,5)/b40);
	% eq. a20a
	hc40	= pdm(6,5)*pdl_1(10);
	zc40	= pdm(5,5)*pdl_1(9);
	d(5)	= d(5)*ccor(alt, rl, hc40, zc40);
end

% **** hydrogen density ****

% eq. a18
g1		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(6,:));
db01	= pdm(1,6)*exp(g1)*pd(6,1);
% eq. a13 - a17
[tz, d(7)] = denss(alt, db01, tinf, tlb, 1.0, -.40, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
if ((alt <= altl(7)) && (sw(15) ~= 0.0))
	% eq. a19
	zh01	= pdm(3,6);
	[tz, b01] = denss(zh01, db01, tinf, tlb, 1.0 - xmm, -1.40, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm01] = denss(alt, b01, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm01	= zhm28;
	d(7)	= dnet(d(7), dm01, zhm01, xmm, 1.0);
	% eq. a20b
	rl		= log(b28*pdm(2,6)*abs(pdl_1(18))/b01);
	% eq. a20a
	hc01	= pdm(6,6)*pdl_1(12);
	zc01	= pdm(5,6)*pdl_1(11);
	d(7)	= d(7)*ccor(alt, rl, hc01, zc01);
	% eq. a21
	hcc01	= pdm(8,6)*pdl_1(20);
	zcc01	= pdm(7,6)*pdl_1(19);
	rc01	= pdm(4,6)*pdl_1(21);
	d(7)	= d(7)*ccor(alt, rc01, hcc01, zcc01);
end

% **** atomic nitrogen density ****

% eq. a18
g14		= sw(21)*globe5(doy, sec, glat, glong, stl, f107a, f107, ap, pd(7,:));
db14	= pdm(1,7)*exp(g14)*pd(7,1);
% eq. a13 - a17
[tz, d(8)] = denss(alt, db14, tinf, tlb, 14.0, 0.0, ptm(6), s, t0, za, z0, tr12);
t(2)	= tz;
if ((alt <= altl(8)) && (sw(15) ~= 0.0))
	% eq. a19
	zh14	= pdm(3,7);
	[tz, b14] = denss(zh14, db14, tinf, tlb, 14.0 - xmm, -1.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	[tz, dm14] = denss(alt, b14, tinf, tlb, xmm, 0.0, ptm(6), s, t0, za, z0, tr12);
	t(2)	= tz;
	% eq. a12
	zhm14	= zhm28;
	d(8)	= dnet(d(8), dm14, zhm14, xmm, 14.0);
	% eq. a20b
	rl		= log(b28*pdm(2,7)*abs(pdl_0(3))/b14);
	% eq. a20a
	hc14	= pdm(6,7)*pdl_0(2);
	zc14	= pdm(5,7)*pdl_0(1);
	d(8)	= d(8)*ccor(alt, rl, hc14, zc14);
	% eq. a21
	hcc14	= pdm(8,7)*pdl_0(5);
	zcc14	= pdm(7,7)*pdl_0(4);
	rc14	= pdm(4,7)*pdl_0(6);
	d(8)	= d(8)*ccor(alt, rc14, hcc14, zcc14);
end

%  total mass density

d(6)	= 1.66e-24*(4.0*d(1) + 16.0*d(2) + 28.0*d(3) + 32.0*d(4) +...
	40.0*d(5) + d(7) + 14.*d(8));

end


function [tz, denss] = denss (alt, dlb, tinf, tlb, xm, alpha, zlb, s2, t0,...
						za, z0, tr12)
global gsurf re

% calculate temperature and density profiles for msis models
rgas = 831.4;

denss	= 1.;
z		= max(alt, za);
% eq. a4a
zg2		= zeta(z, zlb);
% eq. a1a
tt		= tinf - (tinf - tlb)*exp(-s2*zg2);
ta		= tt;
tz		= tt;
denss	= tz;
if (alt < za)
	% eq. a4b
	zg0	= zeta(z0, za);
	% eq. a2b
	dta	= (tinf - ta)*s2*((re + zlb)/(re + za))^2;
	% eq. a3e
	t12	= t0 + tr12*(ta - t0);
	% eq. a4b
	zg1	= zeta(alt, za);
	% calculate temperature below za
	% eq. a3a
	dd	= 0.666666*zg0*dta/ta/ta - 3.11111*(1./ta - 1./t0) + 7.11111*(1./t12 - 1./t0);
	% eq. a3b
	cc	= zg0*dta/(2.*ta*ta) - (1.0/ta - 1.0/t0) - 2.0*dd;
	% eq. a3c
	bb	= (1.0/ta - 1.0/t0) - cc - dd;
	% eq. a3d
	x	= (zg0-zg1)/zg0;
	% eq. a1b
	x2	= x*x;
	tz	= 1.0/(1.0/t0 + bb*x2 + cc*x2*x2 + dd*x2*x2*x2);
	denss	= tz;
	taf	= (t12 - t0)/(ta - t0);
end
if (xm ~= 0.0)
	if ((ta <= 0.0) || (tz <= 0.0))
		% write(6,*)alt,xm,tinf,tlb,t0,ta,ii,jg,n,dv(j),ifun,s2,zg0,tz
		tt	= tlb;
		ta	= tlb;
		tz	= tlb;
    end
	% calculate density above za
	% eq. a17a
	glb	= gsurf/(1.0 + zlb/re)^2;
	% eq. a16a
	gamma	= xm*glb/(s2*rgas*tinf);
	% eq. a13, a14a, & a15
	densa	= dlb*(tlb/tt)^(1.0 + alpha + gamma)*exp(-s2*gamma*zg2);
	denss	= densa;
	if (alt < za)
		% calculate density below za
		% eq. a17b
		glb	= gsurf/((1.0 + za/re)^2);
		% eq. a16b
		gamm	= xm*glb*zg0/rgas;
		% eq. a13, a14b, & a15
		denss	= densa*(ta/tz)^(1.0 + alpha)*...
			exp(gamm*((x - 1.0)/t0 + bb*(x*x2 - 1.0)/3.0 + cc*(x2*x2*x - 1.0)/5.0 +...
			dd*(x2*x2*x2*x - 1.0)/7.0));
    end
end
% write(6,100)cxm,alt,za,tinf,tlb,s2,t0,s1,ta,tz,dlb,densa,denss
% 100 format(' d',1p13e10.2)

end


function out = zeta (zz, zl)

global re

out = (zz - zl)*(re + zl)/(re + zz);

end


%function globe5(yrd,sec,lat,long,tloc,f107a,f107,ap,p)
% calculate g(l) function for msis-86/cira 1986
% upper thermosphere parameters
function tinf = globe5 (yrd, sec, lat, along, tloc, f107a, f107, ap, p)

global sw swc plg dfa ctloc stloc c2tloc s2tloc c3tloc s3tloc apdf day_ apt

nsw = 14;
t=zeros(15,1);
dgtr = 1.74533e-2;
dr = 1.72142e-2;
xl = 1000.;
tll = 1000.;
dayl = -1.0;
p14 = -1000.0;
p18 = -1000.0;
p32 = -1000.0;
hr = .2618;
sr = 7.2722e-5;
p39 = -1000.0;

day_	= yrd - fix(yrd/1000.0)*1000.0;

t(10)	= 0.0;
t(11)	= 0.0;
t(12)	= 0.0;
t(13)	= 0.0;

% eq. a22 (remainder of code)
if (xl ~= lat)
    % calculate legendre polynomials
    c	= sin(lat*dgtr);
	s	= cos(lat*dgtr);
	c2	= c*c;
	c4	= c2*c2;
	s2	= s*s;

	plg(2,1) = c;
	plg(3,1) = 0.5*(3.0*c2 - 1.0);
	plg(4,1) = 0.50*(5.0*c*c2 - 3.0*c);
	plg(5,1) = (35.0*c4 - 30.0*c2 + 3.0)/8.0;
	plg(6,1) = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0;
	plg(7,1) = (11.0*c*plg(6,1) - 5.0*plg(5,1))/6.0;
	plg(2,2) = s;
	plg(3,2) = 3.0*c*s;
	plg(4,2) = 1.50*(5.0*c2 - 1.0)*s;
	plg(5,2) = 2.50*(7.0*c2*c - 3.0*c)*s;
	plg(6,2) = 1.8750*(21.0*c4 - 14.0*c2 +1.0)*s;
	plg(7,2) = (11.0*c*plg(6,2) - 6.0*plg(5,2))/5.0;
	plg(3,3) = 3.0*s2;
	plg(4,3) = 15.0*s2*c;
	plg(5,3) = 7.50*(7.0*c2 - 1.0)*s2;
	plg(6,3) = 3.0*c*plg(5,3) - 2.0*plg(4,3);
	plg(7,3) = (11.0*c*plg(6,3) - 7.0*plg(5,3))/4.0;
	plg(8,3) = (13.0*c*plg(7,3) - 8.0*plg(6,3))/5.0;
	plg(4,4) = 15.0*s2*s;
	plg(5,4) = 105.0*s2*s*c;
	plg(6,4) = (9.0*c*plg(5,4) - 7.0*plg(4,4))/2.0;
	plg(7,4) = (11.0*c*plg(6,4) - 8.0*plg(5,4))/3.0;
	xl		  = lat;
end

if (tll ~= tloc)
    stloc  = sin(hr*tloc);
	ctloc  = cos(hr*tloc);
	s2tloc = sin(2.0*hr*tloc);
	c2tloc = cos(2.0*hr*tloc);
	s3tloc = sin(3.0*hr*tloc);
	c3tloc = cos(3.0*hr*tloc);
	tll    = tloc;
end

if ((day_ ~= dayl) || (p(14) ~= p14))
    cd14  = cos(dr*(day_ - p(14)));
	c2d14 = cos(dr*2.0*(day_ - p(14)));
	p14 = p(14);
end
if ((day_ ~= dayl) || (p(18) ~= p18))
	cd18  = cos(2.0*dr*(day_ - p(18)));
	p18 = p(18);
end
if ((day_ ~= dayl) || (p(32) ~= p32))
	cd32  = cos(dr*(day_ - p(32)));
	p32 = p(32);
end
if ((day_ ~= dayl) || (p(39) ~= p39))
	cd39  = cos(2.0*dr*(day_ - p(39)));
	p39 = p(39);
end
dayl = day_;

% f10.7 effect
df   = f107 - f107a;
dfa  = f107a - 150.0;
t(1) = p(20)*df + p(21)*df*df + p(22)*dfa + p(30)*dfa*dfa;
f1   = 1.0 + (p(48)*dfa + p(20)*df + p(21)*df*df)*swc(1);
f2   = 1.0 + (p(50)*dfa + p(20)*df + p(21)*df*df)*swc(1);

% time independent
t(2) = (p(2)*plg(3,1) + p(3)*plg(5,1) + p(23)*plg(7,1)) + (p(15)*plg(3,1))*dfa*swc(1)...
	+ p(27)*plg(2,1);

% symmetrical annual
t(3) = p(19)*cd32;

% symmetrical semiannual
t(4) = (p(16) + p(17)*plg(3,1))*cd18;

% asymmetrical annual
t(5) =  f1*(p(10)*plg(2,1) + p(11)*plg(4,1))*cd14;

% asymmetrical semiannual
t(6) =  p(38)*plg(2,1)*cd39;

% diurnal
t71  = (p(12)*plg(3,2) + p(36)*plg(2,2))*cd14*swc(5);
t72  = (p(13)*plg(3,2) + p(37)*plg(2,2))*cd14*swc(5);
t(7) = f2*((p(4)*plg(2,2) + p(5)*plg(4,2) + p(28)*plg(6,2) + t71)*ctloc...
	 + (p(7)*plg(2,2) + p(8)*plg(4,2) +p(29)*plg(6,2) + t72)*stloc);

% semidiurnal
t81  = (p(24)*plg(4,3))*cd14*swc(5);
t82  = (p(34)*plg(4,3))*cd14*swc(5);
t(8) = f2*((p(6)*plg(3,3) + p(42)*plg(5,3) + t81)*c2tloc + (p(9)*plg(3,3)...
	+ p(43)*plg(5,3) + t82)*s2tloc);

% terdiurnal
t(14) = f2*((p(40)*plg(3,3) + (p(94)*plg(5,4) + p(47)*plg(7,4))*cd14*swc(5))*s3tloc...
	  + (p(41)*plg(4,4) + (p(95)*plg(5,4) + p(49)*plg(7,4))*cd14*swc(5))*c3tloc);

% magnetic activity based on daily ap
if ((sw(9) == -1.0) && (p(52) ~= 0.0))
    exp1 = exp(-10800.0*abs(p(52))/(1.0 + p(139)*(45.0 - abs(lat))));
    if (exp1 > .99999)
        exp1 = .99999;
    end
    exp2 = exp(-10800.0*abs(p(54)));
    if (exp2 > .99999)
        exp2 = .99999;
    end
    if (p(25) < 1.e-4)
        p(25) = 1.e-4;
    end
    apt(1) = sg0(exp1, ap, p);
    apt(3) = sg0(exp2, ap, p);
	t(9)   = apt(1)*(p(51) + p(97)*plg(3,1) + p(55)*plg(5,1) +...
		(p(126)*plg(2,1) + p(127)*plg(4,1) + p(128)*plg(6,1))*cd14*swc(5)+...
		(p(129)*plg(2,2) + p(130)*plg(4,2) +...
		p(131)*plg(6,2))*swc(7)*cos(hr*(tloc - p(132))));
else
    apd  = (ap(1) - 4.0);
    p44  = p(44);
    p45  = p(45);
    if (p44 < 0.0)
        p44 = 1.e-5;
    end
    apdf = (apd + (p45 - 1.0)*(apd + (exp(-p44*apd) - 1.0)/p44));
	t(9) = apdf*(p(33) + p(46)*plg(3,1) + p(35)*plg(5,1) +...
		(p(101)*plg(2,1) + p(102)*plg(4,1) + p(103)*plg(6,1))*cd14*swc(5) +...
		(p(122)*plg(2,2) + p(123)*plg(4,2) +...
		p(124)*plg(6,2))*swc(7)*cos(hr*(tloc - p(125))));
end

if ((sw(10) ~= 0.0) && (along > -1000.0))
	% longitudinal
	t(11) = (1.0 + p(90)*plg(2,1))*(1.0 + p(81)*dfa*swc(1))*...
		((p(65)*plg(3,2) + p(66)*plg(5,2) + p(67)*plg(7,2) +...
		p(104)*plg(2,2) + p(105)*plg(4,2) + p(106)*plg(6,2) +...
		swc(5)*(p(110)*plg(2,2) + p(111)*plg(4,2) +...
		p(112)*plg(6,2))*cd14)* cos(dgtr*along) +...
		(p(91)*plg(3,2) + p(92)*plg(5,2) + p(93)*plg(7,2) +...
		p(107)*plg(2,2) + p(108)*plg(4,2) + p(109)*plg(6,2) +...
		swc(5)*(p(113)*plg(2,2) + p(114)*plg(4,2) +...
		p(115)*plg(6,2))*cd14)*sin(dgtr*along));
	% ut and mixed ut,longitude
	t(12) = (1.0 + p(96)*plg(2,1))*(1.0 + p(82)*dfa*swc(1))*...
		(1.0 + p(120)*plg(2,1)*swc(5)*cd14)*...
		((p(69)*plg(2,1) + p(70)*plg(4,1) +...
		p(71)*plg(6,1))*cos(sr*(sec-p(72)))) +...
		swc(11)*(p(77)*plg(4,3) + p(78)*plg(6,3) +...
		p(79)*plg(8,3))*cos(sr*(sec - p(80)) +...
		2.0*dgtr*along)*(1.0 + p(138)*dfa*swc(1));
	% ut,longitude magnetic activity
	if ((sw(9) == -1.0) && (p(52) ~= 0.0))
		t(13) = apt(1)*swc(11)*(1.0 + p(133)*plg(2,1))*...
			((p(53)*plg(3,2) + p(99)*plg(5,2) + p(68)*plg(7,2))*...
			cos(dgtr*(along - p(98)))) + apt(1)*swc(11)*swc(5)*...
			(p(134)*plg(2,2) + p(135)*plg(4,2) + p(136)*plg(6,2))*...
			cd14*cos(dgtr*(along - p(137))) + apt(1)*swc(12)*...
			(p(56)*plg(2,1) + p(57)*plg(4,1) + p(58)*plg(6,1))*...
			cos(sr*(sec - p(59)));
    else
		t(13) = apdf*swc(11)*(1.0 + p(121)*plg(2,1))*...
			((p(61)*plg(3,2) + p(62)*plg(5,2) + p(63)*plg(7,2))*...
			cos(dgtr*(along - p(64)))) + apdf*swc(11)*swc(5)*...
			(p(116)*plg(2,2) + p(117)*plg(4,2) + p(118)*plg(6,2))*...
			cd14*cos(dgtr*(along - p(119))) + apdf*swc(12)*...
			(p(84)*plg(2,1) + p(85)*plg(4,1) + p(86)*plg(6,1))*...
			cos(sr*(sec - p(76)));
    end
	% parms not used: 60,83,100,140-150
end
tinf = 0.0;
if (sw(9) == -1.0)
    tinf = p(31);
end

for i = 1:nsw
    tinf = tinf + abs(sw(i))*t(i);
end
  
end


function out = g0 (a, p25, p26)
	% eq. a24d
	% = (a - 4.0+(p(26)-1.0)*(a - 4.0+(exp(-abs(p(25))*(a - 4.0)) - 1.0)/abs(p(25))));
	out = a - 4.0 + (p26 - 1.0)*(a - 4.0 + (exp(-abs(p25)*(a - 4.0)) - 1.0)/abs(p25));
end


function out = sumex (ex)
	% eq. a24c
	% sumex(ex)=1.d0+(1.d0-ex**19)/(1.d0-ex)*ex**(.5d0);
	out = 1.0 + (1.0 - ex^19)/(1.0 - ex)*sqrt(ex);
end


function out = sg0 (ex, ap, p)
	% eq. a24a
	%	sg0(ex)=(g0(ap(2))+(g0(ap(3))*ex+g0(ap(4))*ex*ex+g0(ap(5))*ex**3
	%	+(g0(ap(6))*ex**4+g0(ap(7))*ex**12)*(1.d0-ex**8)/(1.d0-ex)))/sumex(ex);
	out = (g0(ap(1), p(24), p(25)) + ...
		(g0(ap(2), p(24), p(25))*ex + ...
		g0(ap(3), p(24), p(25))*ex*ex + g0(ap(4), p(24), p(25))*ex*ex*ex + ...
		(g0(ap(5), p(24), p(25))*(ex^4) + ...
		g0(ap(6), p(24), p(25))*(ex^12))*(1.0 - (ex^8))/(1.0 - ex)))/sumex(ex);
end


% limited parameter version of globe 9/2/82
% calculate g(l) function for msis-86/cira 1986
% lower thermosphere parameters
function tt = glob5l(p)

global sw swc plg dfa ctloc stloc c2tloc s2tloc c3tloc s3tloc apdf day_ apt

dr = 1.72142e-2;
t = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
dayl = -1.0;
p7 = -1000.0;
p9 = -1000.0;
p11 = -1000.0;

if ((day_ ~= dayl) || (p7 ~= p(7)))
    cd7 = cos(dr*(day_ - p(7)));
end
if ((day_ ~= dayl) || (p9 ~= p(9)))
    cd9 = cos(2.0*dr*(day_ - p(9)));
end
if ((day_ ~= dayl) || (p11 ~= p(11)))
    cd11 = cos(dr*(day_ - p(11)));
end

dayl	= day_;
p7		= p(7);
p9		= p(9);
p11		= p(11);

t(1)	= p(2)*dfa;
t(2)	= p(4)*plg(3,1);
t(3)	= p(6)*cd7;
t(4)	= p(8)*cd9;
t(5)	= (p(10)*plg(2,1) + p(22)*plg(4,1))*cd11;
t(6)	= 0.0;
t(7)	= p(14)*plg(2,2)*ctloc + p(15)*plg(2,2)*stloc;
t(8)	= (p(16)*plg(3,3) + p(18)*plg(5,3) + (p(20)*plg(6,3))*cd11*swc(5))*c2tloc +...
	(p(17)*plg(3,3) + p(19)*plg(5,3) + (p(21)*plg(6,3))*cd11*swc(5))*s2tloc;
t(14)	= p(12)*plg(4,4)*c3tloc + p(25)*plg(4,4)*s3tloc;

if (sw(9) == 1.0)
    t(9)  = apdf*(p(23) + p(24)*plg(3,1)*swc(2));
end
if (sw(9) == -1.0)
    t(9) = p(3)*apt(3) + p(5)*plg(3,1)*apt(3)*swc(2);
end

%parms not used: 13
tt		= 0.0;

for i = 1:14
    tt = tt + abs(sw(i))*t(i);
end

end


function dnet = dnet(dd, dm, zhm, xmm, xm)

% turbopause correction for msis models
% eq. a12b

a = zhm/(xmm - xm);

% eq. a12a
ylog = a*log(dm/dd);

if (ylog < -10.0)
    dnet = dd;
else
    if (ylog > 10.0)
        dnet = dm;
    else
        dnet = dd*(1.0 + exp(ylog))^(1/a);
    end
end

end


function out = ccor(alt, r, h1, zh)

% chemistry/dissociation correction for msis models
% eq. a20a or eq. a21
e = (alt - zh)/h1;
	
if (e > 70.0)
    ccor = 0.0;
else
    if (e < -70.0)
        ccor = r;
    else
        ex   = exp(e);
        ccor = r/(1.0 + ex);
    end
end

out = exp(ccor);

end

