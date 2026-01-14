clc

% Defining model filepaths
modIn = "./src/dynamics/6U CubeSat.obj";
modOut = "./src/dynamics/6U CubeSat.mat";

inparam.gsi_model = 'cook';

[~,matName,~] = fileparts(modOut);

% Load mesh parameters
load(fiName,'meshdata');
x = meshdata.XData;
y = meshdata.YData;
z = meshdata.ZData;
areas = meshdata.Areas;
surfN = meshdata.SurfN;
barC = meshdata.BariC;
LenRef = meshdata.Lref;
matID = meshdata.MatID;