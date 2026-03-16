function visibilityInfo = analyze_visibility(tspan, visibilityHistory)
    visibilityChanges = diff([0; visibilityHistory(:)]);
    idxStartVisible = find(visibilityChanges == 1, 1, 'first');
    idxEndVisible = find(visibilityChanges == -1, 1, 'first');

    visibilityInfo = struct();
    visibilityInfo.visibility_changes = visibilityChanges;
    visibilityInfo.idx_start_visible = idxStartVisible;
    visibilityInfo.idx_end_visible = idxEndVisible;
    visibilityInfo.t_start_visible = NaN;
    visibilityInfo.t_end_visible = NaN;

    if ~isempty(idxStartVisible)
        visibilityInfo.t_start_visible = tspan(idxStartVisible);
    end
    if ~isempty(idxEndVisible)
        visibilityInfo.t_end_visible = tspan(idxEndVisible);
    end
end
