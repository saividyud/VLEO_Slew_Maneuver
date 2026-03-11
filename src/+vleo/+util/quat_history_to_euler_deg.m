function eulerHistoryDeg = quat_history_to_euler_deg(qHistoryScalarFirst)
    nRows = size(qHistoryScalarFirst, 1);
    eulerHistoryDeg = zeros(nRows, 3);

    for rowIdx = 1:nRows
        eulerHistoryDeg(rowIdx, :) = vleo.util.quat_to_euler_deg(qHistoryScalarFirst(rowIdx, :));
    end
end
