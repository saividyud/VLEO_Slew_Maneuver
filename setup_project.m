function setup_project(verbose)
% SETUP_PROJECT Adds the project tree to the MATLAB path.
% It recursively adds folders from the repository root while skipping
% hidden/editor folders and personal workspaces.

    if nargin < 1 || isempty(verbose)
        verbose = false;
    elseif ~(isscalar(verbose) && (islogical(verbose) || isnumeric(verbose)))
        error('VLEO_Slew_Maneuver:InvalidSetupVerbosity', ...
            'setup_project verbosity flag must be a numeric or logical scalar.');
    else
        verbose = logical(verbose);
    end

    persistent cached_root_dir cached_path_entries

    root_dir = fileparts(mfilename('fullpath'));
    setup_marker_name = 'VLEO_Slew_Maneuver_SetupRoot';
    excluded_dirs = {'assets', 'data', 'docs', 'legacy', 'simulations', 'templates', 'tools', 'workspaces'};

    if isempty(cached_path_entries) || ~strcmp(cached_root_dir, root_dir)
        cached_path_entries = get_project_paths(root_dir, excluded_dirs);
        cached_root_dir = root_dir;
    end

    if is_setup_marker_valid(root_dir, setup_marker_name) && ...
            are_project_paths_available(cached_path_entries)
        return;
    end

    if are_project_paths_available(cached_path_entries)
        setappdata(0, setup_marker_name, root_dir);
        return;
    end

    if verbose
        fprintf('Setting up project paths...\n');
    end

    for i = 1:numel(cached_path_entries)
        candidate = cached_path_entries{i};
        if isempty(candidate) || is_path_entry_present(candidate)
            continue;
        end

        addpath(candidate);
        if verbose
            fprintf('  Added: %s\n', candidate);
        end
    end

    if verbose
        fprintf('Project setup complete.\n');
    end

    setappdata(0, setup_marker_name, root_dir);
end

function path_entries = get_project_paths(root_dir, excluded_dirs)
    raw_entries = strsplit(genpath(root_dir), pathsep);
    keep_mask = true(size(raw_entries));

    for i = 1:numel(raw_entries)
        candidate = raw_entries{i};
        keep_mask(i) = ~isempty(candidate) && ...
            ~should_skip_path(candidate, root_dir, excluded_dirs);
    end

    path_entries = raw_entries(keep_mask);
end

function tf = are_project_paths_available(path_entries)
    tf = true;

    for i = 1:numel(path_entries)
        if ~is_path_entry_present(path_entries{i})
            tf = false;
            return;
        end
    end
end

function tf = is_path_entry_present(candidate)
    current_path_entries = strsplit(path, pathsep);
    tf = any(strcmp(current_path_entries, candidate));
end

function tf = is_setup_marker_valid(root_dir, setup_marker_name)
    tf = isappdata(0, setup_marker_name) && ...
        strcmp(getappdata(0, setup_marker_name), root_dir);
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

    tf = any(startsWith(relative_parts, '.')) || ...
        any(startsWith(relative_parts, '+')) || ...
        any(startsWith(relative_parts, '@')) || ...
        any(ismember(relative_parts, excluded_dirs));
end
