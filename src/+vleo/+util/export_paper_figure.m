function exportInfo = export_paper_figure(figHandle, fileStem, varargin)
    parser = inputParser();
    parser.addRequired('figHandle', @(x) ishghandle(x));
    parser.addRequired('fileStem', @(x) (ischar(x) || isstring(x)) && isscalar(string(x)));
    parser.addParameter('subdir', '', @(x) ischar(x) || isstring(x));
    parser.addParameter('resolution', 300, @(x) isscalar(x) && isnumeric(x) && isfinite(x) && x > 0);
    parser.parse(figHandle, fileStem, varargin{:});

    fileStem = char(string(parser.Results.fileStem));
    subdir = strtrim(char(string(parser.Results.subdir)));

    if isempty(subdir)
        figureDir = vleo.util.paper_figure_dir();
    else
        figureDir = vleo.util.paper_figure_dir(subdir);
    end

    pngPath = fullfile(figureDir, strcat(fileStem, '.png'));
    pdfPath = fullfile(figureDir, strcat(fileStem, '.pdf'));

    exportgraphics(figHandle, pngPath, 'Resolution', parser.Results.resolution);
    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector');

    exportInfo = struct();
    exportInfo.directory = figureDir;
    exportInfo.pngPath = pngPath;
    exportInfo.pdfPath = pdfPath;
end
