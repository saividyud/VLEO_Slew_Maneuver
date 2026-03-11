function modelName = normalize_gsi_model(modelName)
    modelName = lower(strtrim(char(modelName)));
    if ~any(strcmp(modelName, {'cook', 'sentman'}))
        modelName = 'cook';
    end
end
