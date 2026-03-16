function omegaBodyHistory = eci_rates_to_body_history(omegaEciHistory, qRows)
    nSteps = size(omegaEciHistory, 1);
    omegaBodyHistory = zeros(nSteps, 3);

    for i = 1:nSteps
        rEciToBody = quat2dcm(qRows(i, :));
        omegaBodyHistory(i, :) = (rEciToBody * omegaEciHistory(i, :)')';
    end
end
