function rHist = body_to_eci_history_from_quaternion_rows(qRows)
    nSteps = size(qRows, 1);
    rHist = zeros(3, 3, nSteps);

    for i = 1:nSteps
        rHist(:, :, i) = quat2dcm(qRows(i, :))';
    end
end
