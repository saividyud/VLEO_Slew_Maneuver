function qRows = normalize_quaternion_rows(qRows)
    for rowIdx = 1:size(qRows, 1)
        qNorm = norm(qRows(rowIdx, :));
        if qNorm > 0
            qRows(rowIdx, :) = qRows(rowIdx, :) / qNorm;
        else
            qRows(rowIdx, :) = [1, 0, 0, 0];
        end
    end
end
