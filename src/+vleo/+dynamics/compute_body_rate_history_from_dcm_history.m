function omegaBodyHistory = compute_body_rate_history_from_dcm_history(rHist, dt)
    nSteps = size(rHist, 3);
    omegaBodyHistory = zeros(nSteps, 3);
    if nSteps == 0
        return;
    end
    if nSteps == 1
        return;
    end

    rDot = zeros(3, 3, nSteps);
    for i = 2:nSteps-1
        rDot(:, :, i) = (rHist(:, :, i + 1) - rHist(:, :, i - 1)) / (2 * dt);
    end
    rDot(:, :, 1) = (rHist(:, :, 2) - rHist(:, :, 1)) / dt;
    rDot(:, :, nSteps) = (rHist(:, :, nSteps) - rHist(:, :, nSteps - 1)) / dt;

    for i = 1:nSteps
        skewW = rHist(:, :, i)' * rDot(:, :, i);
        omegaBodyHistory(i, 1) = skewW(3, 2);
        omegaBodyHistory(i, 2) = skewW(1, 3);
        omegaBodyHistory(i, 3) = skewW(2, 1);
    end
end
