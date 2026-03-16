function set_visible_domain(ax, tStartVisible, tEndVisible)
    if ~isnan(tStartVisible) && ~isnan(tEndVisible)
        xlim(ax, [tStartVisible, tEndVisible] / 60);
    end
end
