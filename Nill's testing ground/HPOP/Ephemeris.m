%--------------------------------------------------------------------------
%
% Ephemeris computation using variable-order Radau IIA integrator with
% step-size control
%
% Last modified:   2023/08/24   Meysam Mahooti
%
%--------------------------------------------------------------------------
function Eph = Ephemeris(Y0, N_Step, Step)

% h = 0.01; % Step-size of integration [s]
% span = 0:Step:N_Step*Step;
% num = length(span);
% y = zeros(num,length(Y0));
% % Integration from t=t_0 to t=t_end
% y(1,:) = Y0;
% for ii = 1:num-1
%     [y_f, ~, h_next] = Runge_Kutta_Fehlberg_7_8(@Accel,y(ii,:)',span(ii),h,span(ii+1),1.0e-10);
%     h = h_next;
%     y(ii+1,:) = y_f;
% end
% Eph(:,1) = span;
% Eph(:,2:7) = y;

options = odeset(RelTol=1e-10,AbsTol=1e-12);
[t,yout] = ode113(@Accel,(0:Step:N_Step*Step),Y0,options);
% options = rdpset('RelTol',1e-13,'AbsTol',1e-16);
% [t,yout] = radau(@Accel,(0:Step:N_Step*Step),Y0,options);
Eph(:,1) = t;
Eph(:,2:7) = yout;

