function qRows = align_quaternion_signs(qRows)
    if isempty(qRows)
        return;
    end

    for rowIdx = 2:size(qRows, 1)
        if dot(qRows(rowIdx - 1, :), qRows(rowIdx, :)) < 0
            qRows(rowIdx, :) = -qRows(rowIdx, :);
        end
    end
end
