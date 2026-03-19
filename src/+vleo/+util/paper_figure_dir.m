function dirPath = paper_figure_dir(varargin)
    dirPath = vleo.util.project_path( ...
        'docs', 'Preparation_of_Papers_for_AIAA_Technical_Conferences', 'figures', varargin{:});
    vleo.util.ensure_directory(dirPath);
end
