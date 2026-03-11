function rootDir = project_root()
    persistent cachedRootDir

    if isempty(cachedRootDir) || ~isfolder(cachedRootDir)
        candidate = fileparts(mfilename('fullpath'));

        while ~isfile(fullfile(candidate, 'setup_project.m'))
            parentDir = fileparts(candidate);
            if strcmp(parentDir, candidate)
                error('VLEO_Slew_Maneuver:SetupNotFound', ...
                    'Could not locate setup_project.m from the project utilities package.');
            end
            candidate = parentDir;
        end

        cachedRootDir = candidate;
    end

    rootDir = cachedRootDir;
end
