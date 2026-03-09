vleoProjectRootAutoSetup = fileparts(mfilename('fullpath'));

while ~isfile(fullfile(vleoProjectRootAutoSetup, 'setup_project.m'))
    vleoParentDirAutoSetup = fileparts(vleoProjectRootAutoSetup);
    if strcmp(vleoParentDirAutoSetup, vleoProjectRootAutoSetup)
        error('VLEO_Slew_Maneuver:SetupNotFound', ...
            'Could not locate setup_project.m for automatic project setup.');
    end
    vleoProjectRootAutoSetup = vleoParentDirAutoSetup;
end

vleoCurrentPathEntriesAutoSetup = strsplit(path, pathsep);
if ~any(strcmp(vleoCurrentPathEntriesAutoSetup, vleoProjectRootAutoSetup))
    addpath(vleoProjectRootAutoSetup);
end

setup_project();

clear vleoCurrentPathEntriesAutoSetup vleoParentDirAutoSetup vleoProjectRootAutoSetup
