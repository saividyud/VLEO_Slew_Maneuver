projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
vleo.gui.open_simulation_gui();
