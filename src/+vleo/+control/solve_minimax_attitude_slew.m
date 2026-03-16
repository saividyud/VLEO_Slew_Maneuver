function result = solve_minimax_attitude_slew(qInitial, omegaInitialBody, qTarget, ...
        omegaTargetBody, inertiaBody, maneuverDuration, nIntervals)
    nNodes = nIntervals + 1;
    timeNodes = linspace(0, maneuverDuration, nNodes)';

    [tauPositiveGuess, ~, switchFractionGuess] = bang_bang_initial_guess( ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody, maneuverDuration);

    x0 = [tauPositiveGuess; switchFractionGuess];
    residualFun = @(x) fast_bang_bang_terminal_residual(x, qInitial, omegaInitialBody, ...
        qTarget, omegaTargetBody, inertiaBody, maneuverDuration);

    solverUsed = 'analytic bang-bang approximation';
    exitflag = 0;
    output = struct('message', 'Used analytic bang-bang approximation for fast slew generation.');
    xOpt = x0;

    try
        fsOpts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 400, ...
            'MaxFunctionEvaluations', 4000, 'FunctionTolerance', 1e-10, ...
            'StepTolerance', 1e-14, 'OptimalityTolerance', 1e-10);
        [xSolved, ~, fsFlag, fsOutput] = fsolve(residualFun, x0, fsOpts);
        xSolved(4:6) = min(max(xSolved(4:6), 0.01), 0.99);

        if fsFlag > 0
            xOpt = xSolved;
            solverUsed = sprintf('bang-bang with fsolve refinement (exitflag=%d)', fsFlag);
            exitflag = fsFlag;
            output = fsOutput;
        else
            residualOpt = norm(residualFun(xSolved));
            residualGuess = norm(residualFun(x0));
            if residualOpt < residualGuess
                xOpt = xSolved;
                solverUsed = sprintf('bang-bang with partial fsolve refinement (exitflag=%d)', fsFlag);
                exitflag = fsFlag;
                output = fsOutput;
            end
        end
    catch ME
        warning('control_test4:fsolveUnavailable', ...
            'fsolve unavailable (%s); using analytic guess only.', ME.message);
    end

    tauPositive = xOpt(1:3);
    switchFractions = xOpt(4:6);
    tauNodes = build_bang_bang_torque_nodes(timeNodes, tauPositive, switchFractions);
    [qNodes, omegaNodes] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, timeNodes);
    ceqOpt = [terminal_quaternion_error_vector(qNodes(end, :)', qTarget); ...
        omegaNodes(end, :)' - omegaTargetBody];
    cOpt = vecnorm(tauNodes, 2, 2) - norm(tauPositive);
    gamma = norm(tauPositive);

    timeVerify = linspace(0, maneuverDuration, max(200, 8 * nIntervals + 1))';
    [qVerify, omegaVerifyBody] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, timeVerify);

    result = struct();
    result.timeNodes = timeNodes;
    result.qNodes = qNodes;
    result.omegaNodesBody = omegaNodes;
    result.tauNodes = tauNodes;
    result.interpMethod = 'previous';
    result.solverUsed = solverUsed;
    result.gamma = gamma;
    result.exitflag = exitflag;
    result.output = output;
    result.message = output.message;
    result.maxTorqueNodeNorm = max(vecnorm(tauNodes, 2, 2));
    result.maxEqualityResidual = max(abs(ceqOpt));
    result.maxInequalityViolation = max(max(cOpt), 0);
    result.timeVerify = timeVerify;
    result.qVerify = qVerify;
    result.omegaVerifyBody = omegaVerifyBody;
    result.terminalAngleErrorDeg = vleo.util.rotation_error_angle_deg(qVerify(end, :)', qTarget);
    result.terminalRateErrorDegPerSec = norm(rad2deg(omegaVerifyBody(end, :)' - omegaTargetBody));
end

function [tauPositiveGuess, tauNegativeGuess, switchFractionGuess] = bang_bang_initial_guess( ...
        qInitial, omegaInitialBody, qTarget, omegaTargetBody, inertiaBody, maneuverDuration)
    relativeRotationVector = relative_rotation_vector_body(qInitial, qTarget);
    principalInertia = diag(inertiaBody);
    tauPositiveGuess = zeros(3, 1);
    tauNegativeGuess = zeros(3, 1);
    switchFractionGuess = 0.5 * ones(3, 1);

    for axisIdx = 1:3
        [alphaAmplitude, switchTime] = solve_bang_bang_1d(0, omegaInitialBody(axisIdx), ...
            relativeRotationVector(axisIdx), omegaTargetBody(axisIdx), maneuverDuration);
        tauPositiveGuess(axisIdx) = principalInertia(axisIdx) * alphaAmplitude;
        tauNegativeGuess(axisIdx) = -principalInertia(axisIdx) * alphaAmplitude;
        switchFractionGuess(axisIdx) = vleo.util.clamp_scalar(switchTime / maneuverDuration, 0.05, 0.95);
    end
end

function residual = fast_bang_bang_terminal_residual(x, qInitial, omegaInitialBody, ...
        qTarget, omegaTargetBody, inertiaBody, maneuverDuration)
    tauPositive = x(1:3);
    switchFractions = x(4:6);
    [qFinal, omegaFinal] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, tauPositive, ...
        switchFractions, inertiaBody, maneuverDuration, maneuverDuration);
    residual = [terminal_quaternion_error_vector(qFinal(end, :)', qTarget); ...
        omegaFinal(end, :)' - omegaTargetBody];
end

function [qHistory, omegaHistory] = simulate_piecewise_bang_bang(qInitial, omegaInitialBody, ...
        tauPositive, switchFractions, inertiaBody, maneuverDuration, outputTimes)
    outputTimes = reshape(outputTimes, [], 1);
    outputTimes = sort(outputTimes);
    if isempty(outputTimes) || abs(outputTimes(1)) > 1e-12
        outputTimes = [0; outputTimes];
    end
    if abs(outputTimes(end) - maneuverDuration) > 1e-12
        outputTimes = [outputTimes; maneuverDuration];
    end

    switchTimes = sort(unique(maneuverDuration * switchFractions(:)));
    segmentBoundaries = sort(unique([0; switchTimes; outputTimes; maneuverDuration]));

    qHistory = zeros(numel(outputTimes), 4);
    omegaHistory = zeros(numel(outputTimes), 3);
    state = [qInitial / norm(qInitial); omegaInitialBody];
    outputIdx = 1;
    qHistory(outputIdx, :) = state(1:4)';
    omegaHistory(outputIdx, :) = state(5:7)';
    outputIdx = outputIdx + 1;

    for segIdx = 1:numel(segmentBoundaries) - 1
        t0 = segmentBoundaries(segIdx);
        t1 = segmentBoundaries(segIdx + 1);
        if t1 <= t0
            continue;
        end

        tauSegment = symmetric_bang_bang_torque(0.5 * (t0 + t1), tauPositive, switchFractions, maneuverDuration);
        dtSegment = t1 - t0;
        nSubsteps = max(1, ceil(dtSegment / 0.25));
        h = dtSegment / nSubsteps;

        for stepIdx = 1:nSubsteps
            state = rk4_attitude_step(state, h, tauSegment, inertiaBody);
            state(1:4) = state(1:4) / norm(state(1:4));
        end

        while outputIdx <= numel(outputTimes) && abs(outputTimes(outputIdx) - t1) < 1e-10
            qHistory(outputIdx, :) = state(1:4)';
            omegaHistory(outputIdx, :) = state(5:7)';
            outputIdx = outputIdx + 1;
        end
    end

    qHistory = vleo.util.normalize_quaternion_rows(qHistory);
    qHistory = vleo.util.align_quaternion_signs(qHistory);
end

function tauBody = symmetric_bang_bang_torque(t, tauPositive, switchFractions, maneuverDuration)
    switchTimes = maneuverDuration * switchFractions(:);
    tauBody = tauPositive(:);
    for axisIdx = 1:3
        if t > switchTimes(axisIdx)
            tauBody(axisIdx) = -tauBody(axisIdx);
        end
    end
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

function [amplitude, switchTime] = solve_bang_bang_1d(x0, v0, xf, vf, tau)
    if abs(vf - v0) < 1e-12
        switchTime = tau / 2;
        amplitude = 4 * (xf - x0 - v0 * tau) / tau^2;
        return;
    end

    candidateGuesses = [0.3, 0.5, 0.7, 0.25, 0.75] * tau;
    amplitude = 0;
    switchTime = tau / 2;

    for guess = candidateGuesses
        try
            root = fzero(@(ts) bang_bang_switch_residual(ts, x0, v0, xf, vf, tau), guess);
            if root > 0 && root < tau
                denominator = 2 * root - tau;
                if abs(denominator) > 1e-12
                    amplitude = (vf - v0) / denominator;
                    switchTime = root;
                    return;
                end
            end
        catch
        end
    end

    denominator = 2 * switchTime - tau;
    if abs(denominator) > 1e-12
        amplitude = (vf - v0) / denominator;
    end
end

function residual = bang_bang_switch_residual(switchTime, x0, v0, xf, vf, tau)
    denominator = 2 * switchTime - tau;
    if abs(denominator) < 1e-12
        residual = 1e10;
        return;
    end

    amplitude = (vf - v0) / denominator;
    shapeTerm = 2 * switchTime * tau - switchTime^2 - 0.5 * tau^2;
    xTau = x0 + v0 * tau + amplitude * shapeTerm;
    residual = xTau - xf;
end

function tauNodes = build_bang_bang_torque_nodes(timeNodes, tauPositive, switchFractions)
    nNodes = numel(timeNodes);
    maneuverDuration = timeNodes(end);
    switchTimes = maneuverDuration * switchFractions(:)';
    tauNodes = zeros(nNodes, 3);

    for nodeIdx = 1:nNodes
        t = timeNodes(nodeIdx);
        for axisIdx = 1:3
            if t <= switchTimes(axisIdx)
                tauNodes(nodeIdx, axisIdx) = tauPositive(axisIdx);
            else
                tauNodes(nodeIdx, axisIdx) = -tauPositive(axisIdx);
            end
        end
    end
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
