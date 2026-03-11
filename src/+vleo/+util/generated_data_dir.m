function dirPath = generated_data_dir()
    dirPath = vleo.util.project_path('data', 'generated');
    vleo.util.ensure_directory(dirPath);
end
