function layout = corner_thruster_layout(params, layoutMode, varargin)
    if nargin < 2 || isempty(layoutMode)
        layoutMode = "standard";
    end

    parser = inputParser();
    parser.addParameter('minTranslationForcePerAxis', 1.5, @(x) isscalar(x) && isfinite(x) && x > 0);
    parser.addParameter('gridResolution', 300, @(x) isscalar(x) && isfinite(x) && x >= 10);
    parser.parse(varargin{:});

    layoutMode = normalize_layout_mode(layoutMode);

    dims_m = reshape(params.bodyDims, 1, []);
    if numel(dims_m) ~= 3
        error('VLEO_Slew_Maneuver:InvalidBodyDims', ...
            'params.bodyDims must contain exactly 3 elements.')
    end

    inertiaDiag = diag(params.I_CB).';
    if numel(inertiaDiag) ~= 3 || any(~isfinite(inertiaDiag)) || any(inertiaDiag <= 0)
        error('VLEO_Slew_Maneuver:InvalidInertia', ...
            'params.I_CB must be a positive 3x3 inertia matrix.')
    end

    a = dims_m(1) / 2;
    b = dims_m(2) / 2;
    c = dims_m(3) / 2;

    signs = [+1 +1 +1;  +1 +1 -1;  +1 -1 +1;  +1 -1 -1; ...
             -1 +1 +1;  -1 +1 -1;  -1 -1 +1;  -1 -1 -1];

    if layoutMode == "standard"
        exhaustDirectionUnit = [1, 1, 1] / sqrt(3);
        score = min(4 * [abs(c * exhaustDirectionUnit(2) - b * exhaustDirectionUnit(3)) / inertiaDiag(1), ...
                         abs(a * exhaustDirectionUnit(3) - c * exhaustDirectionUnit(1)) / inertiaDiag(2), ...
                         abs(b * exhaustDirectionUnit(1) - a * exhaustDirectionUnit(2)) / inertiaDiag(3)]);
    else
        [exhaustDirectionUnit, score] = optimize_direction(a, b, c, inertiaDiag, ...
            parser.Results.minTranslationForcePerAxis, parser.Results.gridResolution);
    end

    pos = signs .* [a, b, c];
    exhaustDirections = signs .* exhaustDirectionUnit;
    forceDirections = -exhaustDirections;

    nThrusters = size(signs, 1);
    wrenchMatrix = zeros(6, nThrusters);
    for idx = 1:nThrusters
        wrenchMatrix(:, idx) = [forceDirections(idx, :).'; cross(pos(idx, :), forceDirections(idx, :)).'];
    end

    torqueArmsPerThruster = [abs(c * exhaustDirectionUnit(2) - b * exhaustDirectionUnit(3)), ...
                             abs(a * exhaustDirectionUnit(3) - c * exhaustDirectionUnit(1)), ...
                             abs(b * exhaustDirectionUnit(1) - a * exhaustDirectionUnit(2))];

    layout = struct();
    layout.mode = char(layoutMode);
    layout.signs = signs;
    layout.positions_m = pos;
    layout.exhaust_direction_unit = exhaustDirectionUnit(:);
    layout.exhaust_directions = exhaustDirections;
    layout.force_directions = forceDirections;
    layout.wrench_matrix = wrenchMatrix;
    layout.translation_net_per_unit_thrust = 4 * exhaustDirectionUnit(:);
    layout.torque_arms_per_thruster_m = torqueArmsPerThruster(:);
    layout.torque_net_per_unit_thrust_Nm = 4 * torqueArmsPerThruster(:);
    layout.angular_accel_per_unit_thrust = layout.torque_net_per_unit_thrust_Nm ./ inertiaDiag(:);
    layout.min_translation_force_per_axis = parser.Results.minTranslationForcePerAxis;
    layout.optimization_grid_resolution = parser.Results.gridResolution;
    layout.minimax_score = score;
end

function layoutMode = normalize_layout_mode(layoutMode)
    if islogical(layoutMode)
        if layoutMode
            layoutMode = "optimized";
        else
            layoutMode = "standard";
        end
        return;
    end

    layoutMode = lower(string(layoutMode));
    if ~ismember(layoutMode, ["standard", "optimized"])
        error('VLEO_Slew_Maneuver:InvalidThrusterLayoutMode', ...
            'layoutMode must be ''standard'' or ''optimized''.')
    end
end

function [exhaustDirectionUnit, score] = optimize_direction(a, b, c, inertiaDiag, minTranslationForcePerAxis, gridResolution)
    azimuth = linspace(0, pi / 2, gridResolution);
    elevation = linspace(0, pi / 2, gridResolution);
    [azimuthGrid, elevationGrid] = meshgrid(azimuth, elevation);

    dx = cos(elevationGrid) .* cos(azimuthGrid);
    dy = cos(elevationGrid) .* sin(azimuthGrid);
    dz = sin(elevationGrid);

    alphaX = 4 * abs(c .* dy - b .* dz) / inertiaDiag(1);
    alphaY = 4 * abs(a .* dz - c .* dx) / inertiaDiag(2);
    alphaZ = 4 * abs(b .* dx - a .* dy) / inertiaDiag(3);

    scoreGrid = min(min(alphaX, alphaY), alphaZ);
    validMask = (4 * dx >= minTranslationForcePerAxis) & ...
        (4 * dy >= minTranslationForcePerAxis) & ...
        (4 * dz >= minTranslationForcePerAxis);
    scoreGrid(~validMask) = -inf;

    [score, bestIndex] = max(scoreGrid(:));
    if ~isfinite(score)
        error('VLEO_Slew_Maneuver:ThrusterOrientationOptimizationFailed', ...
            'No feasible thruster orientation satisfies the translation constraint.')
    end

    exhaustDirectionUnit = [dx(bestIndex), dy(bestIndex), dz(bestIndex)];
end
