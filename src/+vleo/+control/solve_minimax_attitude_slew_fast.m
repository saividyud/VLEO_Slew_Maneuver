function result = solve_minimax_attitude_slew_fast(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, nIntervals, varargin)
% SOLVE_MINIMAX_ATTITUDE_SLEW_FAST Solve a fixed-time principal-axis minimax slew.
% This solver certifies the minimum peak principal-axis torque for a
% diagonal-inertia, principal-axis double-integrator attitude model.

    solveTimer = tic;

    if nargin < 7 || isempty(nIntervals)
        nIntervals = 20;
    end

    parser = inputParser;
    addParameter(parser, 'refineNonlinear', false, @islogical);
    addParameter(parser, 'nRefinementIntervals', max(8, min(16, nIntervals)), @is_valid_interval_vector);
    addParameter(parser, 'momentArms', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==3 && all(isfinite(x)) && all(x>0)));
    addParameter(parser, 'objectiveMode', 'force', @(x) any(strcmpi(string(x), {'force', 'torque'})));
    parse(parser, varargin{:});

    reportMomentArms = parser.Results.momentArms;
    if isempty(reportMomentArms)
        reportMomentArms = ones(3, 1);
    else
        reportMomentArms = reshape(reportMomentArms, 3, 1);
    end

    objectiveMode = char(validatestring(char(string(parser.Results.objectiveMode)), {'force', 'torque'}));
    if strcmp(objectiveMode, 'torque')
        optimizationMomentArms = ones(3, 1);
        objectiveText = 'minimize max_i sup_t |tau_i(t)|';
    else
        optimizationMomentArms = reportMomentArms;
        objectiveText = 'minimize max_i sup_t |tau_i(t) / L_i|';
    end

    qInitial = normalize_quaternion(qInitial, 'solve_minimax_attitude_slew_fast:InvalidInitialQuaternion');
    qTarget = normalize_quaternion(qTarget, 'solve_minimax_attitude_slew_fast:InvalidTargetQuaternion');
    omegaInitialBody = validate_body_rate(omegaInitialBody, ...
        'solve_minimax_attitude_slew_fast:InvalidInitialBodyRate');
    omegaTargetBody = validate_body_rate(omegaTargetBody, ...
        'solve_minimax_attitude_slew_fast:InvalidTargetBodyRate');
    principalInertia = validate_inertia(inertiaBody);

    validateattributes(maneuverDuration, {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, 'maneuverDuration');
    validateattributes(nIntervals, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 1}, ...
        mfilename, 'nIntervals');

    relativeRotationVectorBody = relative_rotation_vector_body(qInitial, qTarget);
    axisPlans = repmat(empty_axis_plan(), 3, 1);

    for axisIdx = 1:3
        axisPlans(axisIdx) = solve_axis_plan(relativeRotationVectorBody(axisIdx), ...
            omegaInitialBody(axisIdx), omegaTargetBody(axisIdx), principalInertia(axisIdx), maneuverDuration);
    end

    gamma = max([axisPlans.peakTorqueAbs]' ./ optimizationMomentArms);
    switchTimes = [axisPlans.switchTime]';
    switchFractions = switchTimes / maneuverDuration;

    timeNodes = build_time_nodes(maneuverDuration, nIntervals, switchTimes);
    [rotationVectorNodes, omegaNodesBody, tauNodes] = sample_axis_plans(axisPlans, timeNodes);

    timeVerify = build_time_nodes(maneuverDuration, max(200, 10 * nIntervals), switchTimes);
    [rotationVectorVerify, omegaVerifyBody, tauVerify] = sample_axis_plans(axisPlans, timeVerify);

    qNodes = relative_rotation_history_to_quaternions(qInitial, rotationVectorNodes);
    qVerify = relative_rotation_history_to_quaternions(qInitial, rotationVectorVerify);

    certification = struct();
    certification.objective = objectiveText;
    certification.model = 'principal-axis double-integrator with diagonal inertia';
    certification.lowerBound = gamma;
    certification.upperBound = gamma;
    certification.gap = 0;
    certification.isCertified = true;
    certification.assumptions = { ...
        'inertiaBody is diagonal in the commanded body frame', ...
        'attitude dynamics are planned in principal-axis rotation-vector coordinates', ...
        'each axis uses the one-switch minimax solution of the double-integrator boundary-value problem'};

    result = struct();
    result.model = certification.model;
    result.solverUsed = 'analytic principal-axis minimax';
    result.message = 'Solved analytically and certified for the diagonal principal-axis model.';
    result.exitflag = 1;
    result.output = struct('message', result.message, 'iterations', 0, ...
        'funcCount', 0, 'algorithm', 'analytic principal-axis minimax');
    result.interpMethod = 'previous';
    result.solveTimeSec = toc(solveTimer);
    result.gamma = gamma;
    result.gammaLowerBound = gamma;
    result.gammaUpperBound = gamma;
    result.optimalityCertificate = certification;
    result.isCertifiedMinimax = certification.isCertified;
    result.relativeRotationVectorBody = relativeRotationVectorBody;
    result.rotationAngleDeg = rad2deg(norm(relativeRotationVectorBody));
    result.timeNodes = timeNodes;
    result.qNodes = qNodes;
    result.rotationVectorNodes = rotationVectorNodes;
    result.omegaNodesBody = omegaNodesBody;
    result.tauNodes = tauNodes;
    result.timeVerify = timeVerify;
    result.qVerify = qVerify;
    result.rotationVectorVerify = rotationVectorVerify;
    result.omegaVerifyBody = omegaVerifyBody;
    result.tauVerify = tauVerify;
    result.switchTimes = switchTimes;
    result.switchFractions = switchFractions;
    result.axisPlans = axisPlans;
    result.objectiveMode = objectiveMode;
    result.momentArms = reportMomentArms;
    result.optimizationMomentArms = optimizationMomentArms;
    result.peakAxisTorqueAbs = reshape([axisPlans.peakTorqueAbs], [], 1);
    result.peakAxisForceAbs = result.peakAxisTorqueAbs ./ reportMomentArms;
    result.gammaTorqueEquivalent = max(result.peakAxisTorqueAbs);
    result.gammaForceEquivalent = max(result.peakAxisForceAbs);
    result.peakAxisAccelerationAbs = reshape([axisPlans.peakAccelerationAbs], [], 1);
    result.maxTorqueNodeComponentAbs = max(max(abs(tauNodes)));
    result.maxTorqueVerifyComponentAbs = max(max(abs(tauVerify)));
    result.maxTorqueNodeNorm = max(vecnorm(tauNodes, 2, 2));
    result.maxTorqueVerifyNorm = max(vecnorm(tauVerify, 2, 2));
    result.terminalAngleErrorDeg = vleo.util.rotation_error_angle_deg(qVerify(end, :)', qTarget);
    result.terminalRateErrorDegPerSec = norm(rad2deg(omegaVerifyBody(end, :)' - omegaTargetBody));
    result.terminalRotationVectorError = rotationVectorVerify(end, :)' - relativeRotationVectorBody;
    result.maxEqualityResidual = max(abs([result.terminalRotationVectorError; ...
        omegaVerifyBody(end, :)' - omegaTargetBody]));
    result.maxInequalityViolation = max(max(max(abs(tauVerify), [], 1)' ./ optimizationMomentArms - gamma), 0);

    if parser.Results.refineNonlinear
        result = refine_nonlinear_solution(result, qInitial, omegaInitialBody, qTarget, ...
            omegaTargetBody, inertiaBody, maneuverDuration, parser.Results.nRefinementIntervals, ...
            optimizationMomentArms, reportMomentArms, objectiveMode);
    end
end

function result = refine_nonlinear_solution(result, qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, nRefinementIntervals, ...
        optimizationMomentArms, reportMomentArms, objectiveMode)
    refineTimer = tic;
    intervalCandidates = unique(reshape(nRefinementIntervals, 1, []));
    angleToleranceDeg = 5e-2;
    rateToleranceDegPerSec = 5e-3;
    bestMetric = inf;
    bestResult = result;
    bestResult.refinement = struct('attempted', true, 'success', false, 'attempts', []);
    tauGuessSourceTime = result.timeNodes;
    tauGuessSourceNodes = result.tauNodes;
    attemptInfo = repmat(struct('nIntervals', 0, 'optimizerExitflag', 0, 'accepted', false, ...
        'terminalAngleErrorDeg', 0, 'terminalRateErrorDegPerSec', 0, 'gamma', 0), numel(intervalCandidates), 1);

    try
        options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
            'MaxIterations', 200, 'MaxFunctionEvaluations', 12000, ...
            'ConstraintTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-12);
    catch ME
        result.solverUsed = sprintf('%s (nonlinear refinement unavailable)', result.solverUsed);
        result.message = sprintf('%s Nonlinear refinement skipped: %s', result.message, ME.message);
        result.refinement = struct('attempted', true, 'success', false, 'message', ME.message, ...
            'solveTimeSec', toc(refineTimer));
        return;
    end

    for candidateIdx = 1:numel(intervalCandidates)
        nIntervalsCandidate = intervalCandidates(candidateIdx);
        timeNodes = linspace(0, maneuverDuration, nIntervalsCandidate + 1)';
        intervalMidpoints = 0.5 * (timeNodes(1:end - 1) + timeNodes(2:end));
        tauGuess = interp1(tauGuessSourceTime, tauGuessSourceNodes, intervalMidpoints, 'previous', 'extrap');
        gammaGuess = max(max(abs(tauGuess) ./ optimizationMomentArms'));
        x0 = [max(gammaGuess, result.gamma); tauGuess(:)];
        objectiveFun = @(x) x(1);
        nonlinearConstraintFun = @(x) refinement_constraints(x, qInitial, omegaInitialBody, ...
            qTarget, omegaTargetBody, inertiaBody, maneuverDuration, optimizationMomentArms);
        lb = [0; -inf(size(tauGuess(:)))];

        try
            [xOpt, ~, optimizerExitflag, output] = fmincon(objectiveFun, x0, [], [], [], [], lb, [], ...
                nonlinearConstraintFun, options);
        catch
            attemptInfo(candidateIdx) = struct('nIntervals', nIntervalsCandidate, 'optimizerExitflag', -999, ...
                'accepted', false, 'terminalAngleErrorDeg', inf, ...
                'terminalRateErrorDegPerSec', inf, 'gamma', inf);
            continue;
        end

        tauIntervals = reshape(xOpt(2:end), nIntervalsCandidate, 3);
        candidateResult = build_refined_result(result, tauIntervals, qInitial, omegaInitialBody, ...
            qTarget, omegaTargetBody, inertiaBody, maneuverDuration, timeNodes, output, optimizerExitflag, ...
            optimizationMomentArms, reportMomentArms, objectiveMode);
        candidateMetric = refinement_metric(candidateResult.terminalAngleErrorDeg, ...
            candidateResult.terminalRateErrorDegPerSec);
        accepted = candidateResult.terminalAngleErrorDeg <= angleToleranceDeg && ...
            candidateResult.terminalRateErrorDegPerSec <= rateToleranceDegPerSec;

        attemptInfo(candidateIdx) = struct('nIntervals', nIntervalsCandidate, ...
            'optimizerExitflag', optimizerExitflag, 'accepted', accepted, ...
            'terminalAngleErrorDeg', candidateResult.terminalAngleErrorDeg, ...
            'terminalRateErrorDegPerSec', candidateResult.terminalRateErrorDegPerSec, ...
            'gamma', candidateResult.gamma);

        if candidateMetric < bestMetric
            bestMetric = candidateMetric;
            bestResult = candidateResult;
        end

        tauGuessSourceTime = candidateResult.timeNodes;
        tauGuessSourceNodes = candidateResult.tauNodes;

        if accepted
            if optimizerExitflag <= 0
                bestResult.exitflag = 1;
                bestResult.message = sprintf(['Analytic warm start reached the requested terminal state ', ...
                    'after nonlinear refinement with %d intervals (optimizer exitflag=%d).'], ...
                    nIntervalsCandidate, optimizerExitflag);
            end
            bestResult.refinement = struct('attempted', true, 'success', true, ...
                'solveTimeSec', toc(refineTimer), 'attempts', attemptInfo(1:candidateIdx), ...
                'nIntervalsAccepted', nIntervalsCandidate);
            bestResult.solveTimeSec = result.solveTimeSec + toc(refineTimer);
            bestResult.isCertifiedMinimax = false;
            result = bestResult;
            return;
        end
    end

    bestResult.solverUsed = sprintf('%s + nonlinear direct shooting refinement', bestResult.solverUsed);
    bestResult.message = sprintf(['Analytic warm start refinement did not meet the requested terminal ', ...
        'tolerance after trying interval counts [%s].'], num2str(intervalCandidates));
    bestResult.refinement = struct('attempted', true, 'success', false, ...
        'solveTimeSec', toc(refineTimer), 'attempts', attemptInfo);
    bestResult.solveTimeSec = result.solveTimeSec + toc(refineTimer);
    bestResult.isCertifiedMinimax = false;
    result = bestResult;
end

function result = build_refined_result(baseResult, tauIntervals, qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, timeNodes, output, optimizerExitflag, ...
        optimizationMomentArms, reportMomentArms, objectiveMode)
    [qNodes, omegaNodesBody] = simulate_piecewise_constant_torque(qInitial, omegaInitialBody, tauIntervals, ...
        inertiaBody, timeNodes);
    timeVerify = linspace(0, maneuverDuration, max(200, 10 * size(tauIntervals, 1) + 1))';
    [qVerify, omegaVerifyBody] = simulate_piecewise_constant_torque(qInitial, omegaInitialBody, tauIntervals, ...
        inertiaBody, timeVerify);
    tauNodes = sample_piecewise_constant_torque(tauIntervals, timeNodes, timeNodes);
    tauVerify = sample_piecewise_constant_torque(tauIntervals, timeNodes, timeVerify);
    gammaRefined = max(max(abs(tauIntervals) ./ optimizationMomentArms'));
    terminalQuaternionError = terminal_quaternion_error_vector(qVerify(end, :)', qTarget);
    terminalRateError = omegaVerifyBody(end, :)' - omegaTargetBody;

    result = baseResult;
    result.timeNodes = timeNodes;
    result.qNodes = qNodes;
    result.omegaNodesBody = omegaNodesBody;
    result.tauNodes = tauNodes;
    result.timeVerify = timeVerify;
    result.qVerify = qVerify;
    result.omegaVerifyBody = omegaVerifyBody;
    result.tauVerify = tauVerify;
    result.gamma = gammaRefined;
    result.objectiveMode = objectiveMode;
    result.peakAxisTorqueAbs = max(abs(tauIntervals), [], 1)';
    result.peakAxisForceAbs = result.peakAxisTorqueAbs ./ reportMomentArms;
    result.gammaTorqueEquivalent = max(result.peakAxisTorqueAbs);
    result.gammaForceEquivalent = max(result.peakAxisForceAbs);
    result.maxTorqueNodeComponentAbs = max(max(abs(result.tauNodes)));
    result.maxTorqueVerifyComponentAbs = max(max(abs(tauVerify)));
    result.maxTorqueNodeNorm = max(vecnorm(result.tauNodes, 2, 2));
    result.maxTorqueVerifyNorm = max(vecnorm(tauVerify, 2, 2));
    result.terminalAngleErrorDeg = vleo.util.rotation_error_angle_deg(qVerify(end, :)', qTarget);
    result.terminalRateErrorDegPerSec = norm(rad2deg(terminalRateError));
    result.terminalRotationVectorError = terminalQuaternionError;
    result.maxEqualityResidual = max(abs([terminalQuaternionError; terminalRateError]));
    result.maxInequalityViolation = max(max(max(abs(tauIntervals) ./ optimizationMomentArms') - gammaRefined), 0);
    result.solverUsed = sprintf('%s + nonlinear direct shooting refinement', baseResult.solverUsed);
    result.message = sprintf('Analytic warm start refined with fmincon SQP (exitflag=%d).', optimizerExitflag);
    result.exitflag = optimizerExitflag;
    result.output = output;
end

function metric = refinement_metric(angleErrorDeg, rateErrorDegPerSec)
    metric = angleErrorDeg + 10 * rateErrorDegPerSec;
end

function tf = is_valid_interval_vector(value)
    tf = isnumeric(value) && isvector(value) && ~isempty(value) && all(isfinite(value)) && ...
        all(value >= 2) && all(value == floor(value));
end

function [c, ceq] = refinement_constraints(x, qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, optimizationMomentArms)
    gamma = x(1);
    tauIntervals = reshape(x(2:end), [], 3);
    [qHistory, omegaHistory] = simulate_piecewise_constant_torque(qInitial, omegaInitialBody, tauIntervals, ...
        inertiaBody, maneuverDuration);
    terminalQuaternionError = terminal_quaternion_error_vector(qHistory(end, :)', qTarget);
    terminalRateError = omegaHistory(end, :)' - omegaTargetBody;
    scaledTau = tauIntervals ./ optimizationMomentArms';
    c = abs(scaledTau(:)) - gamma;
    ceq = [terminalQuaternionError; terminalRateError];
end

function tauHistory = sample_piecewise_constant_torque(tauIntervals, timeNodes, sampleTimes)
    sampleTimes = reshape(sampleTimes, [], 1);
    nIntervals = size(tauIntervals, 1);
    tauHistory = zeros(numel(sampleTimes), 3);
    for idx = 1:numel(sampleTimes)
        if sampleTimes(idx) >= timeNodes(end)
            intervalIdx = nIntervals;
        else
            intervalIdx = find(sampleTimes(idx) < timeNodes(2:end), 1, 'first');
        end
        tauHistory(idx, :) = tauIntervals(intervalIdx, :);
    end
end

function [qHistory, omegaHistory] = simulate_piecewise_constant_torque(qInitial, omegaInitialBody, tauIntervals, ...
        inertiaBody, outputTimes)
    outputTimes = reshape(outputTimes, [], 1);
    if isempty(outputTimes) || abs(outputTimes(1)) > 1e-12
        outputTimes = [0; outputTimes];
    end
    maneuverDuration = outputTimes(end);
    nIntervals = size(tauIntervals, 1);
    timeNodes = linspace(0, maneuverDuration, nIntervals + 1)';

    qHistory = zeros(numel(outputTimes), 4);
    omegaHistory = zeros(numel(outputTimes), 3);
    state = [qInitial / norm(qInitial); omegaInitialBody];
    qHistory(1, :) = state(1:4)';
    omegaHistory(1, :) = state(5:7)';

    for outputIdx = 2:numel(outputTimes)
        tPrev = outputTimes(outputIdx - 1);
        tNext = outputTimes(outputIdx);
        segmentBounds = sort(unique([tPrev; timeNodes; tNext]));
        segmentBounds = segmentBounds(segmentBounds >= tPrev & segmentBounds <= tNext);

        for segIdx = 1:numel(segmentBounds) - 1
            t0 = segmentBounds(segIdx);
            t1 = segmentBounds(segIdx + 1);
            if t1 <= t0
                continue;
            end

            intervalIdx = find(t0 < timeNodes(2:end), 1, 'first');
            if isempty(intervalIdx)
                intervalIdx = nIntervals;
            end
            tauBody = tauIntervals(intervalIdx, :)';
            dt = t1 - t0;
            nSubsteps = max(1, ceil(dt / 0.25));
            h = dt / nSubsteps;
            for stepIdx = 1:nSubsteps
                state = rk4_attitude_step(state, h, tauBody, inertiaBody);
                state(1:4) = state(1:4) / norm(state(1:4));
            end
        end

        qHistory(outputIdx, :) = state(1:4)';
        omegaHistory(outputIdx, :) = state(5:7)';
    end

    qHistory = vleo.util.normalize_quaternion_rows(qHistory);
    qHistory = vleo.util.align_quaternion_signs(qHistory);
end

function nextState = rk4_attitude_step(state, h, tauBody, inertiaBody)
    k1 = attitude_state_derivative(state, tauBody, inertiaBody);
    k2 = attitude_state_derivative(state + 0.5 * h * k1, tauBody, inertiaBody);
    k3 = attitude_state_derivative(state + 0.5 * h * k2, tauBody, inertiaBody);
    k4 = attitude_state_derivative(state + h * k3, tauBody, inertiaBody);
    nextState = state + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
end

function stateDot = attitude_state_derivative(state, tauBody, inertiaBody)
    qScalarFirst = state(1:4);
    qScalarFirst = qScalarFirst / norm(qScalarFirst);
    omegaBody = state(5:7);

    stateDot = zeros(7, 1);
    stateDot(1:4) = quaternion_kinematics(qScalarFirst, omegaBody);
    stateDot(5:7) = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody);
end

function qErrorVector = terminal_quaternion_error_vector(qCurrent, qTarget)
    rCurrent = quat2dcm(qCurrent');
    rTarget = quat2dcm(qTarget');
    rError = rTarget * rCurrent';
    qError = dcm2quat(rError)';
    if qError(1) < 0
        qError = -qError;
    end
    qErrorVector = qError(2:4);
end

function axisPlan = empty_axis_plan()
    axisPlan = struct( ...
        'thetaTarget', 0, ...
        'omegaInitial', 0, ...
        'omegaTarget', 0, ...
        'deltaRate', 0, ...
        'deltaPosition', 0, ...
        'initialSign', 1, ...
        'switchTime', 0, ...
        'peakAccelerationAbs', 0, ...
        'peakTorqueAbs', 0, ...
        'torqueBeforeSwitch', 0, ...
        'torqueAfterSwitch', 0);
end

function omegaBody = validate_body_rate(omegaBody, errorId)
    validateattributes(omegaBody, {'numeric'}, {'real', 'finite', 'numel', 3}, mfilename, 'omegaBody');
    omegaBody = reshape(omegaBody, 3, 1);
    if ~all(isfinite(omegaBody))
        error(errorId, 'Body rates must be finite 3-vectors.');
    end
end

function principalInertia = validate_inertia(inertiaBody)
    validateattributes(inertiaBody, {'numeric'}, {'real', 'finite', 'size', [3, 3]}, ...
        mfilename, 'inertiaBody');
    principalInertia = diag(inertiaBody);
    offDiagonal = inertiaBody - diag(principalInertia);
    tolerance = 1e-12 * max(1, norm(inertiaBody, 'fro'));
    if norm(offDiagonal, 'fro') > tolerance
        error('solve_minimax_attitude_slew_fast:NonDiagonalInertia', ...
            'solve_minimax_attitude_slew_fast requires a diagonal inertia matrix.');
    end
    if any(~isfinite(principalInertia)) || any(principalInertia <= 0)
        error('solve_minimax_attitude_slew_fast:InvalidInertia', ...
            'Principal inertias must be finite and strictly positive.');
    end
end

function axisPlan = solve_axis_plan(thetaTarget, omegaInitial, omegaTarget, inertiaAxis, maneuverDuration)
    axisPlan = empty_axis_plan();
    axisPlan.thetaTarget = thetaTarget;
    axisPlan.omegaInitial = omegaInitial;
    axisPlan.omegaTarget = omegaTarget;
    axisPlan.deltaRate = omegaTarget - omegaInitial;
    axisPlan.deltaPosition = thetaTarget - 0.5 * maneuverDuration * (omegaInitial + omegaTarget);

    if abs(axisPlan.deltaPosition) <= 1e-14 && abs(axisPlan.deltaRate) <= 1e-14
        axisPlan.switchTime = 0.5 * maneuverDuration;
        return;
    end

    if abs(axisPlan.deltaPosition) > 1e-14
        axisPlan.initialSign = sign(axisPlan.deltaPosition);
    else
        axisPlan.initialSign = sign(axisPlan.deltaRate);
    end

    if axisPlan.initialSign == 0
        axisPlan.initialSign = 1;
    end

    deltaPositionAbs = abs(axisPlan.deltaPosition);
    deltaRate = axisPlan.deltaRate;
    durationSquared = maneuverDuration^2;

    axisPlan.peakAccelerationAbs = (2 * deltaPositionAbs + ...
        sqrt(4 * deltaPositionAbs^2 + maneuverDuration^2 * deltaRate^2)) / durationSquared;
    axisPlan.peakTorqueAbs = inertiaAxis * axisPlan.peakAccelerationAbs;

    if axisPlan.peakAccelerationAbs <= 1e-14
        axisPlan.switchTime = 0.5 * maneuverDuration;
    else
        axisPlan.switchTime = 0.5 * (maneuverDuration + ...
            deltaRate / (axisPlan.initialSign * axisPlan.peakAccelerationAbs));
        axisPlan.switchTime = min(max(axisPlan.switchTime, 0), maneuverDuration);
    end

    axisPlan.torqueBeforeSwitch = axisPlan.initialSign * axisPlan.peakTorqueAbs;
    axisPlan.torqueAfterSwitch = -axisPlan.torqueBeforeSwitch;
end

function timeNodes = build_time_nodes(maneuverDuration, nIntervals, switchTimes)
    baseNodes = linspace(0, maneuverDuration, nIntervals + 1)';
    switchTimes = reshape(switchTimes, [], 1);
    switchTimes = switchTimes(switchTimes > 0 & switchTimes < maneuverDuration);
    timeNodes = sort(unique([baseNodes; switchTimes]));
end

function [rotationVectorHistory, omegaHistory, tauHistory] = sample_axis_plans(axisPlans, sampleTimes)
    sampleTimes = reshape(sampleTimes, [], 1);
    nTimes = numel(sampleTimes);

    rotationVectorHistory = zeros(nTimes, 3);
    omegaHistory = zeros(nTimes, 3);
    tauHistory = zeros(nTimes, 3);

    for axisIdx = 1:3
        thetaHistory = zeros(nTimes, 1);
        omegaAxisHistory = zeros(nTimes, 1);
        tauAxisHistory = zeros(nTimes, 1);
        axisPlan = axisPlans(axisIdx);

        signBefore = axisPlan.initialSign;
        accelMag = axisPlan.peakAccelerationAbs;
        accelBefore = signBefore * accelMag;
        accelAfter = -accelBefore;
        switchTime = axisPlan.switchTime;

        thetaAtSwitch = axisPlan.omegaInitial * switchTime + 0.5 * accelBefore * switchTime^2;
        omegaAtSwitch = axisPlan.omegaInitial + accelBefore * switchTime;

        beforeMask = sampleTimes <= switchTime;
        afterMask = ~beforeMask;
        tBefore = sampleTimes(beforeMask);
        tAfter = sampleTimes(afterMask) - switchTime;

        thetaHistory(beforeMask) = axisPlan.omegaInitial * tBefore + 0.5 * accelBefore * tBefore.^2;
        omegaAxisHistory(beforeMask) = axisPlan.omegaInitial + accelBefore * tBefore;

        thetaHistory(afterMask) = thetaAtSwitch + omegaAtSwitch * tAfter + 0.5 * accelAfter * tAfter.^2;
        omegaAxisHistory(afterMask) = omegaAtSwitch + accelAfter * tAfter;

        tauAxisHistory(sampleTimes < switchTime) = axisPlan.torqueBeforeSwitch;
        tauAxisHistory(sampleTimes >= switchTime) = axisPlan.torqueAfterSwitch;

        if axisPlan.peakTorqueAbs <= 1e-14
            tauAxisHistory(:) = 0;
        end

        rotationVectorHistory(:, axisIdx) = thetaHistory;
        omegaHistory(:, axisIdx) = omegaAxisHistory;
        tauHistory(:, axisIdx) = tauAxisHistory;
    end
end

function qHistory = relative_rotation_history_to_quaternions(qInitial, rotationVectorHistory)
    nRows = size(rotationVectorHistory, 1);
    qHistory = zeros(nRows, 4);
    rInitial = quat2dcm(qInitial');

    for rowIdx = 1:nRows
        rRelative = rotation_vector_to_dcm(rotationVectorHistory(rowIdx, :)')';
        qCurrent = dcm2quat(rRelative * rInitial)';
        qCurrent = qCurrent / norm(qCurrent);

        if rowIdx > 1 && dot(qCurrent, qHistory(rowIdx - 1, :)') < 0
            qCurrent = -qCurrent;
        end

        qHistory(rowIdx, :) = qCurrent';
    end

    qHistory = vleo.util.normalize_quaternion_rows(qHistory);
    qHistory = vleo.util.align_quaternion_signs(qHistory);
end

function rotationVector = relative_rotation_vector_body(qInitial, qTarget)
    rInitial = quat2dcm(qInitial');
    rTarget = quat2dcm(qTarget');
    rRelative = rTarget * rInitial';
    qRelative = dcm2quat(rRelative)';

    if qRelative(1) < 0
        qRelative = -qRelative;
    end

    vectorNorm = norm(qRelative(2:4));
    if vectorNorm < 1e-12
        rotationVector = [0; 0; 0];
        return;
    end

    angle = 2 * atan2(vectorNorm, qRelative(1));
    axis = qRelative(2:4) / vectorNorm;
    rotationVector = axis * angle;
end

function qDot = quaternion_kinematics(qScalarFirst, omegaBody)
    omegaMatrix = [0, -omegaBody(1), -omegaBody(2), -omegaBody(3); ...
        omegaBody(1), 0, omegaBody(3), -omegaBody(2); ...
        omegaBody(2), -omegaBody(3), 0, omegaBody(1); ...
        omegaBody(3), omegaBody(2), -omegaBody(1), 0];
    qDot = 0.5 * omegaMatrix * qScalarFirst;
end

function omegaDot = rigid_body_omega_dot(omegaBody, tauBody, inertiaBody)
    omegaDot = inertiaBody \ (tauBody - cross(omegaBody, inertiaBody * omegaBody));
end

function r = rotation_vector_to_dcm(rotationVector)
    angle = norm(rotationVector);
    kx = skew_matrix(rotationVector);

    if angle < 1e-12
        r = eye(3) + kx + 0.5 * (kx * kx);
        return;
    end

    axis = rotationVector / angle;
    k = skew_matrix(axis);
    r = eye(3) + sin(angle) * k + (1 - cos(angle)) * (k * k);
end

function s = skew_matrix(v)
    s = [0, -v(3), v(2); ...
         v(3), 0, -v(1); ...
         -v(2), v(1), 0];
end

function q = normalize_quaternion(q, errorId)
    validateattributes(q, {'numeric'}, {'real', 'finite', 'numel', 4}, mfilename, 'q');
    q = reshape(q, 4, 1);
    qNorm = norm(q);
    if qNorm <= 0 || ~isfinite(qNorm)
        error(errorId, 'Quaternion must have finite, nonzero norm.');
    end
    q = q / qNorm;
    if q(1) < 0
        q = -q;
    end
end
