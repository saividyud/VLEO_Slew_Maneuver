%---------------------------------------------------------------------------
%
% Density_HP: Computes the atmospheric density for the modified
%             Harris-Priester model.
%
% Inputs:
%   r_Sun       geocentric equatorial position referred to the International
% 				Celestial Reference Frame (ICRF) [m]
%   r_TOD       Satellite position vector in TOD system [m]
%
% Output:
%   density     Density [kg/m^3]
%
% Last modified:   2018/01/27   Meysam Mahooti
%
%---------------------------------------------------------------------------
function density = Density_HP(r_Sun,r_TOD)

HPcoeff

upper_limit = 2000;     % Upper height limit [km]
lower_limit = 110;      % Lower height limit [km]
ra_lag      = 0.523599; % Right ascension lag [rad]
n_prm       = 4;        % Harris-Priester parameter 
                        % 2(6) low(high) inclination
F107        = 175;      % solar radio noise flux

switch F107
    case 65
        h     = hpcoef(1,1:3:end);
        c_min = hpcoef(1,2:3:end);
        c_max = hpcoef(1,3:3:end);
    case 75
        h     = hpcoef(2,1:3:end);
        c_min = hpcoef(2,2:3:end);
        c_max = hpcoef(2,3:3:end);
    case 100
        h     = hpcoef(3,1:3:end);
        c_min = hpcoef(3,2:3:end);
        c_max = hpcoef(3,3:3:end);
    case 125
        h     = hpcoef(4,1:3:end);
        c_min = hpcoef(4,2:3:end);
        c_max = hpcoef(4,3:3:end);
    case 150
        h     = hpcoef(5,1:3:end);
        c_min = hpcoef(5,2:3:end);
        c_max = hpcoef(5,3:3:end);
    case 175
        h     = hpcoef(6,1:3:end);
        c_min = hpcoef(6,2:3:end);
        c_max = hpcoef(6,3:3:end);
    case 200
        h     = hpcoef(7,1:3:end);
        c_min = hpcoef(7,2:3:end);
        c_max = hpcoef(7,3:3:end);
    case 225
        h     = hpcoef(8,1:3:end);
        c_min = hpcoef(8,2:3:end);
        c_max = hpcoef(8,3:3:end);
    case 250
        h     = hpcoef(9,1:3:end);
        c_min = hpcoef(9,2:3:end);
        c_max = hpcoef(9,3:3:end);
    case 275
        h     = hpcoef(10,1:3:end);
        c_min = hpcoef(10,2:3:end);
        c_max = hpcoef(10,3:3:end);
end

% Harris-Priester atmospheric density model parameters
% Height [km], minimum density, maximum density [gm/km^3]
N_Coef = 59;

% Satellite height
[~, ~, height] = Geodetic(r_TOD);
height = height/1000;

% outside height model limits
if ( height >= upper_limit || height <= lower_limit )
    warning('satellite''s height is out of limits');
    density = 0;
    return
end

% Sun right ascension, declination
[S_phi, S_theta, ~] = CalcPolarAngles(r_Sun);
ra_Sun  = S_phi;
dec_Sun = S_theta;

% Unit vector u towards the apex of the diurnal bulge
% in inertial geocentric coordinates
c_dec = cos(dec_Sun);
u(1) = c_dec * cos(ra_Sun + ra_lag);
u(2) = c_dec * sin(ra_Sun + ra_lag);
u(3) = sin(dec_Sun);

% Cosine of half angle between satellite position vector and
% apex of diurnal bulge
c_psi2 = 0.5 + 0.5 * dot(r_TOD,u)/norm(r_TOD);

% Height index search and exponential density interpolation
ih = 1;                           % section index reset
for  i=1:N_Coef-1                 % loop over N_Coef height regimes
  if ( height >= h(i) && height < h(i+1) )
    ih = i;                       % ih identifies height section
    break
  end
end

h_min = ( h(ih) - h(ih+1) )/log( c_min(ih+1)/c_min(ih) );
h_max = ( h(ih) - h(ih+1) )/log( c_max(ih+1)/c_max(ih) );

d_min = c_min(ih) * exp( (h(ih)-height)/h_min );
d_max = c_max(ih) * exp( (h(ih)-height)/h_max );

% Density computation
density = d_min + (d_max-d_min)*(c_psi2^n_prm);
density = density * 1e-9; % [kg/m^3]

