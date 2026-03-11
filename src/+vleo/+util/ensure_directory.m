function ensure_directory(dirPath)
    if ~isfolder(dirPath)
        mkdir(dirPath);
    end
end
