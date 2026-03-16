function alphaBodyHistory = compute_angular_acceleration_history(omegaBodyHistory, dt)
    nSteps = size(omegaBodyHistory, 1);
    alphaBodyHistory = zeros(nSteps, 3);
    if nSteps == 0
        return;
    end
    if nSteps == 1
        return;
    end

    for i = 2:nSteps-1
        alphaBodyHistory(i, :) = (omegaBodyHistory(i + 1, :) - omegaBodyHistory(i - 1, :)) / (2 * dt);
    end
    alphaBodyHistory(1, :) = (omegaBodyHistory(2, :) - omegaBodyHistory(1, :)) / dt;
    alphaBodyHistory(end, :) = (omegaBodyHistory(end, :) - omegaBodyHistory(end - 1, :)) / dt;
end
