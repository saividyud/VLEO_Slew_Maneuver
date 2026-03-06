function setup_project()
% SETUP_PROJECT Adds the project tree to the MATLAB path.
% It recursively adds folders from the repository root while skipping
% hidden/editor folders and personal workspaces.

    root_dir = fileparts(mfilename('fullpath'));
    excluded_dirs = {'.git', '.claude', '.vscode', 'workspaces', 'assets', 'docs'};

    path_entries = strsplit(genpath(root_dir), pathsep);

    fprintf('Setting up project paths...\n');
    for i = 1:numel(path_entries)
        candidate = path_entries{i};
        if isempty(candidate) || should_skip_path(candidate, root_dir, excluded_dirs)
            continue;
        end

        addpath(candidate);
        fprintf('  Added: %s\n', candidate);
    end

    fprintf('Project setup complete.\n');
end

function tf = should_skip_path(candidate, root_dir, excluded_dirs)
    if strcmp(candidate, root_dir)
        relative_parts = {};
    else
        prefix = [root_dir filesep];
        if startsWith(candidate, prefix)
            relative_path = candidate(numel(prefix) + 1:end);
        else
            relative_path = candidate;
        end
        relative_parts = strsplit(relative_path, filesep);
    end

    tf = any(ismember(relative_parts, excluded_dirs));
end
