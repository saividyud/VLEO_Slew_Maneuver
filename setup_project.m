function setup_project()
% SETUP_PROJECT Sets up the MATLAB path for the VLEO Slew Maneuver project.
% Run this function once when you start MATLAB to ensure all dependencies
% and utility functions are accessible.

    % Get the root directory of the project
    root_dir = fileparts(mfilename('fullpath'));

    % Define directories to add to path
    dirs_to_add = {
        fullfile(root_dir, 'src');
        fullfile(root_dir, 'src', 'utils');
        fullfile(root_dir, 'src', 'dynamics');
        fullfile(root_dir, 'lib');
        fullfile(root_dir, 'lib', 'ADBSat');
        fullfile(root_dir, 'lib', 'Gyroscope');
        fullfile(root_dir, 'lib', 'HPOP');
        fullfile(root_dir, 'lib', 'SGP4');
        fullfile(root_dir, 'lib', 'Kumar_Examples');
        fullfile(root_dir, 'tests');
        fullfile(root_dir, 'examples');
        fullfile(root_dir, 'data');
    };

    % Add directories
    fprintf('Setting up project paths...\n');
    for i = 1:length(dirs_to_add)
        if exist(dirs_to_add{i}, 'dir')
            addpath(dirs_to_add{i});
            fprintf('  Added: %s\n', dirs_to_add{i});
        else
            fprintf('  Warning: Directory not found: %s\n', dirs_to_add{i});
        end
    end

    % Optional: Add subdirectories of lib if they contain more code
    % addpath(genpath(fullfile(root_dir, 'lib'))); 
    % (Commented out to prevent adding .git folders or other clutter if using genpath on root)

    fprintf('Project setup complete.\n');
end
