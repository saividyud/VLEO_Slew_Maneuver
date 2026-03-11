function dirPath = simulations_dir()
    dirPath = vleo.util.project_path('simulations');
    vleo.util.ensure_directory(dirPath);
end
