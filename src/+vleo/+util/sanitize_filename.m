function name = sanitize_filename(name)
    name = regexprep(char(name), '[^a-zA-Z0-9_\- ]', '_');
    name = strtrim(name);
    name = regexprep(name, '\s+', '_');
    if isempty(name)
        name = 'simulation';
    end
end
