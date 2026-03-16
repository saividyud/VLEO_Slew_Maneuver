function hide_from_legend(graphicsHandle)
    if isprop(graphicsHandle, 'HandleVisibility')
        graphicsHandle.HandleVisibility = 'off';
    end

    try
        graphicsHandle.Annotation.LegendInformation.IconDisplayStyle = 'off';
    catch
    end
end
