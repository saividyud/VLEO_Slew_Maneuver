function tauAeroHistory = estimate_aero_torque_history(tspan, XHistory, rBodyToEciHist, params)
    if numel(tspan) < 2
        tauAeroHistory = zeros(numel(tspan), 3);
        return;
    end

    dt = abs(tspan(2) - tspan(1));
    sampleStride = max(1, round(params.aeroSampleStep / dt));
    sampleIndices = unique([1:sampleStride:numel(tspan), numel(tspan)]);
    tauSamples = zeros(numel(sampleIndices), 3);

    for k = 1:numel(sampleIndices)
        idx = sampleIndices(k);
        rEciToBody = rBodyToEciHist(:, :, idx)';
        tauSamples(k, :) = vleo.dynamics.evaluate_aero_torque(tspan(idx), XHistory(idx, :)', rEciToBody, params)';
    end

    if numel(sampleIndices) == 1
        tauAeroHistory = repmat(tauSamples, numel(tspan), 1);
        return;
    end

    tauAeroHistory = zeros(numel(tspan), 3);
    for axisIdx = 1:3
        tauAeroHistory(:, axisIdx) = interp1(tspan(sampleIndices), tauSamples(:, axisIdx), tspan, 'pchip');
    end
end
