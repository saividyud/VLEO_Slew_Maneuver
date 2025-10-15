clear
clc
close all

%constants
delt = 
Xref = 
Maxitr =
eps = 

% initial state
X0 = 
Xi = X0

hold on
for i = 1:Maxitr

  %Controller block here
  ui = f(err);

  %update step using integrator
  [t,Xi] = ode45(@Sat_template,[0,delt],[Xi,ui];
  plot3(Xi(:,1),Xi(:,2),Xi(:,3),'o');

  %Sensor block here
  Xhat = f(Xi);

  %Calculate error
  err = Xref - Xi;

  % Might need to calculate Jacobian and linearize the dynamics, but dont 
  % worry about it right now
end
