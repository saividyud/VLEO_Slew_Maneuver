% thruster_placement_analysis.m
% Minimum thruster count and optimal placement for full 6-DOF control
% of the 6U CubeSat (10x20x30 cm, 12 kg).
%
% Key result:
%   Theoretical minimum: 7 unidirectional thrusters (Steinitz theorem)
%   Practical optimum:   8 thrusters at box corners, firing along
%                        inward body diagonal.  Each DOF mode fires
%                        exactly 4 thrusters at equal thrust.
%
% Run: matlab -batch "run('scripts/thruster_placement_analysis.m')"

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();
clc; close all;

%% ---- Vehicle -----------------------------------------------------------
params = vleo.dynamics.default_control_test_params(true);
dims_m  = reshape(params.bodyDims, 1, []);
mass_kg = params.mass;
I_body  = params.I_CB;

a = dims_m(1) / 2;   % 0.05 m
b = dims_m(2) / 2;   % 0.10 m
c = dims_m(3) / 2;   % 0.15 m

fprintf('================================================================\n');
fprintf('  THRUSTER PLACEMENT ANALYSIS — 6U CubeSat (%.0fx%.0fx%.0f cm)\n', dims_m*100);
fprintf('================================================================\n');
fprintf('Mass: %.0f kg   I = diag([%.2f, %.2f, %.2f]) kg*m^2\n\n', ...
    mass_kg, I_body(1,1), I_body(2,2), I_body(3,3));

%% ---- Theory ------------------------------------------------------------
fprintf('MINIMUM THRUSTER COUNT\n');
fprintf('  Physical thrusters are unilateral: F >= 0 (push only).\n');
fprintf('  Each thruster produces a 6-D wrench w = [force; torque].\n');
fprintf('  Full 6-DOF requires the wrench set to positively span R^6.\n');
fprintf('  Steinitz theorem: min vectors to positively span R^n = n+1.\n');
fprintf('  => Theoretical minimum = 7 thrusters.\n');
fprintf('  => Practical optimum   = 8 (symmetric, decoupled modes).\n\n');

%% ---- 8-thruster corner placement ---------------------------------------
%  One thruster at each of the 8 box corners.
%  Reaction force on satellite: inward along computed optimal diagonal.

layout = vleo.control.corner_thruster_layout(params, 'optimized');
signs = layout.signs;
nT = size(signs, 1);
pos = layout.positions_m;
directionUnit = layout.exhaust_direction_unit.';
dx = directionUnit(1);
dy = directionUnit(2);
dz = directionUnit(3);
exhaust_dir = layout.exhaust_directions;
force_dir = layout.force_directions;

fprintf('8-THRUSTER CORNER LAYOUT\n');
fprintf('  ID   Position [cm]             Exhaust direction (Plume)\n');
for k = 1:nT
    fprintf('  T%d   (%+5.1f, %+5.1f, %+5.1f)     (%+.3f, %+.3f, %+.3f)\n', ...
        k, pos(k,:)*100, exhaust_dir(k,:));
end
fprintf('\n');

%% ---- Wrench matrix W (6 x 8) ------------------------------------------
%  w_k = [ d_k ; r_k x d_k ]  (using force_dir)

W = layout.wrench_matrix;

fprintf('WRENCH MATRIX  (rows 1-3 = force, rows 4-6 = torque)\n');
fprintf('         T1       T2       T3       T4       T5       T6       T7       T8\n');
labels = {'Fx','Fy','Fz','tx','ty','tz'};
for r = 1:6
    fprintf('  %s ', labels{r});
    fprintf(' %+7.4f', W(r,:));
    fprintf('\n');
end
fprintf('  rank(W) = %d  (need 6 for full 6-DOF)\n\n', rank(W));

%% ---- Firing patterns for every mode ------------------------------------
% Translation: fire the 4 thrusters on the OPPOSITE half of the box.
%   +Fx -> thrusters with sx<0
%
% Rotation: fire the 4 thrusters where the relevant sign-product is correct.
%   +Rx -> sy*sz > 0   (torque_x = (c*dy - b*dz) > 0)
%   +Ry -> sx*sz < 0   (torque_y = (a*dz - c*dx) > 0)
%   +Rz -> sx*sy > 0   (torque_z = (b*dx - a*dy) > 0)

modeLabel = {'+Tx','-Tx', '+Ty','-Ty', '+Tz','-Tz', ...
             '+Rx','-Rx', '+Ry','-Ry', '+Rz','-Rz'};
modeActive = {
    [5 6 7 8]      % +Tx
    [1 2 3 4]      % -Tx
    [3 4 7 8]      % +Ty
    [1 2 5 6]      % -Ty
    [2 4 6 8]      % +Tz
    [1 3 5 7]      % -Tz
    [1 4 5 8]      % +Rx
    [2 3 6 7]      % -Rx
    [2 4 5 7]      % +Ry
    [1 3 6 8]      % -Ry
    [1 2 7 8]      % +Rz
    [3 4 5 6]      % -Rz
};

fprintf('FIRING PATTERNS  (all 4 thrusters at equal thrust F0)\n');
fprintf('  Mode   Active        Net Force (x F0)         Net Torque (x F0) [m]\n');
fprintf('  ----   ----------    ---------------------     -------------------------\n');
for m = 1:12
    idx = modeActive{m};
    F = zeros(nT, 1); F(idx) = 1;
    wr = W * F;
    fS = wr(1:3);
    tS = wr(4:6);
    idStr = sprintf('T%d,', idx); idStr(end) = [];
    fprintf('  %s   %-12s  (%+6.3f, %+6.3f, %+6.3f)     (%+6.3f, %+6.3f, %+6.3f)\n', ...
        modeLabel{m}, idStr, fS, tS);
end
fprintf('\n');

%% ---- Demonstrate the physics of each mode ------------------------------
fprintf('HOW EACH MODE WORKS\n\n');

fprintf('  TRANSLATIONS (example: +Fx = fire T5,T6,T7,T8)\n');
fprintf('  Each active thruster contributes dx*F0 along +x.\n');
fprintf('  The y and z force components cancel: 2 thrusters push +y, 2 push -y.\n');
fprintf('  All torques cancel by the same pairing symmetry.\n');
fprintf('  The opposite 4 thrusters give the reverse direction.\n\n');

fprintf('  ROTATIONS (example: +Rx = fire T1,T4,T5,T8)\n');
fprintf('  These 4 corners share sy*sz > 0, so all have positive x-torque.\n');
fprintf('  Forces cancel: T1+T4 push -x, T5+T8 push +x (net zero).\n');
fprintf('  Similarly: y-forces and z-forces each cancel in pairs.\n');
fprintf('  Result: pure torque about +x, no net force.\n\n');

%% ---- Control authority per unit thrust ---------------------------------
fprintf('CONTROL AUTHORITY (per individual thruster thrust F0)\n\n');

fprintf('  Translation: Fx = %.3f * F0, Fy = %.3f * F0, Fz = %.3f * F0\n\n', 4*dx, 4*dy, 4*dz);

torqueArm = layout.torque_arms_per_thruster_m.';
torqueNet = layout.torque_net_per_unit_thrust_Nm.';
alphaNet  = layout.angular_accel_per_unit_thrust.';

fprintf('  Axis   Torque/F0 [Nm]   Ang. accel/F0 [rad/s^2]\n');
fprintf('  ----   --------------   -----------------------\n');
for ax = 1:3
    fprintf('   %s       %.5f             %.4f\n', char('x'+ax-1), torqueNet(ax), alphaNet(ax));
end
fprintf('\n');
fprintf('  Effective torque arms = (%.2f, %.2f, %.2f) cm\n', torqueArm*100);
fprintf('  Angular acceleration is optimized to minimize the bottleneck axis.\n\n');

%% ---- Seven-thruster feasibility (null-space method) --------------------
fprintf('SEVEN-THRUSTER FEASIBILITY (remove one at a time)\n\n');

testDirs = [eye(6); -eye(6)];
nDirs    = size(testDirs, 1);
sevenOK  = false(nT, 1);

for removeK = 1:nT
    keep  = setdiff(1:nT, removeK);
    W7    = W(:, keep);
    N7    = null(W7);                    % 7x1 null vector (rank 6, 7 cols)
    allOK = true;

    for d = 1:nDirs
        target = testDirs(d, :).';
        Fp     = W7 \ target;            % particular solution (7x1)

        % Find t range such that Fp + t*N7 >= 0
        tLo = -inf;  tHi = inf;
        feasible = true;
        for j = 1:7
            if N7(j) > 1e-12
                tLo = max(tLo, -Fp(j) / N7(j));
            elseif N7(j) < -1e-12
                tHi = min(tHi, -Fp(j) / N7(j));
            else
                if Fp(j) < -1e-10
                    feasible = false;
                    break;
                end
            end
        end
        if ~feasible || tLo > tHi + 1e-10
            allOK = false;
            break;
        end
    end

    sevenOK(removeK) = allOK;
    if allOK
        fprintf('  Remove T%d: FEASIBLE  (7 thrusters still span all 6-DOF)\n', removeK);
    else
        fprintf('  Remove T%d: NOT feasible\n', removeK);
    end
end

fprintf('\n  Result: %d / %d subsets achieve 6-DOF with 7 thrusters.\n', sum(sevenOK), nT);
if any(sevenOK)
    fprintf('  The theoretical minimum (7) IS achievable on this geometry,\n');
    fprintf('  but firing levels become unequal and authority is asymmetric.\n');
end
fprintf('\n');

%% ---- Demonstrate 7- vs 8-thruster firing for one mode ------------------
if any(sevenOK)
    remIdx   = find(sevenOK, 1, 'first');
    keep7    = setdiff(1:nT, remIdx);
    W7       = W(:, keep7);
    N7       = null(W7);

    fprintf('COMPARISON: 8 vs 7 thrusters for +Rx mode (remove T%d)\n\n', remIdx);
    fprintf('  8-thruster solution: T1,T4,T5,T8 each at F0 (equal thrust)\n');

    % Solve for +Rx with 7 thrusters
    target = [0; 0; 0; 1; 0; 0];
    Fp = W7 \ target;
    % Find t that maximizes the minimum component
    tLo = -inf; tHi = inf;
    for j = 1:7
        if N7(j) > 1e-12
            tLo = max(tLo, -Fp(j) / N7(j));
        elseif N7(j) < -1e-12
            tHi = min(tHi, -Fp(j) / N7(j));
        end
    end
    tOpt = (tLo + tHi) / 2;
    F7 = Fp + tOpt * N7;

    fprintf('  7-thruster solution (normalized to same net torque):\n');
    for j = 1:7
        fprintf('    T%d: F = %+.4f\n', keep7(j), F7(j));
    end
    fprintf('  Verification: W7*F = [');
    fprintf('%.4f ', W7 * F7);
    fprintf(']\n');
    fprintf('  Note: unequal thrust levels -> more complex allocation.\n\n');
end

%% ==== VISUALIZATION ====================================================

%% ---- Figure 1: thruster layout -----------------------------------------
figure('Color', 'w', 'Name', 'Thruster Layout', 'Position', [50 50 750 650]);

draw_cubesat(a, b, c);
hold on;

arrowLen = 0.07;
for k = 1:nT
    p = pos(k, :);
    d = exhaust_dir(k, :);
    plot3(p(1), p(2), p(3), 'ko', 'MarkerSize', 9, 'MarkerFaceColor', [0.2 0.55 1]);
    quiver3(p(1), p(2), p(3), d(1)*arrowLen, d(2)*arrowLen, d(3)*arrowLen, 0, ...
        'Color', [0.85 0.15 0.15], 'LineWidth', 2, 'MaxHeadSize', 0.7);
    off = -0.012 * signs(k, :);
    text(p(1)+off(1), p(2)+off(2), p(3)+off(3), sprintf('T%d', k), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

% Body axes
axLen = 0.22;
quiver3(0,0,0, axLen,0,0, 0, 'r', 'LineWidth', 1.5);
quiver3(0,0,0, 0,axLen,0, 0, 'Color', [0 0.6 0], 'LineWidth', 1.5);
quiver3(0,0,0, 0,0,axLen, 0, 'b', 'LineWidth', 1.5);
text(axLen*1.1, 0, 0, 'X', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
text(0, axLen*1.1, 0, 'Y', 'Color', [0 0.6 0], 'FontSize', 12, 'FontWeight', 'bold');
text(0, 0, axLen*1.1, 'Z', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');

hold off;
axis equal; grid on;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('6U CubeSat — 8 Corner Thrusters (diagonal canting)', 'FontSize', 14);
view([1, 1, 1]);

%% ---- Figure 2: 6-DOF firing modes (positive direction only) -----------
posModeIdx   = [1, 3, 5, 7, 9, 11];      % +Tx, +Ty, +Tz, +Rx, +Ry, +Rz
posModeLabel = {'+T_x', '+T_y', '+T_z', '+R_x', '+R_y', '+R_z'};
activeColor  = [0.1 0.75 0.2];

% Color per axis: X=red, Y=green, Z=blue
axClr = [0.85 0.1 0.1; 0 0.6 0; 0.1 0.3 0.9];

figure('Color', 'w', 'Name', '6-DOF Firing Modes', ...
    'Position', [80 80 1500 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for panel = 1:6
    nexttile;
    m = posModeIdx(panel);
    activeIdx = modeActive{m};

    draw_cubesat(a, b, c);
    hold on;

    % Body axes in every panel for reference
    refLen = 0.18;
    quiver3(0,0,0, refLen,0,0, 0, 'Color', axClr(1,:), 'LineWidth', 1.2);
    quiver3(0,0,0, 0,refLen,0, 0, 'Color', axClr(2,:), 'LineWidth', 1.2);
    quiver3(0,0,0, 0,0,refLen, 0, 'Color', axClr(3,:), 'LineWidth', 1.2);
    text(refLen*1.08, 0, 0, 'X', 'Color', axClr(1,:), 'FontSize', 9, 'FontWeight', 'bold');
    text(0, refLen*1.08, 0, 'Y', 'Color', axClr(2,:), 'FontSize', 9, 'FontWeight', 'bold');
    text(0, 0, refLen*1.08, 'Z', 'Color', axClr(3,:), 'FontSize', 9, 'FontWeight', 'bold');

    % Draw all thrusters
    for k = 1:nT
        p = pos(k,:);
        d = exhaust_dir(k,:);
        isOn = ismember(k, activeIdx);

        if isOn
            plot3(p(1), p(2), p(3), 'ko', 'MarkerSize', 9, 'MarkerFaceColor', activeColor);
            quiver3(p(1), p(2), p(3), d(1)*arrowLen, d(2)*arrowLen, d(3)*arrowLen, 0, ...
                'Color', activeColor, 'LineWidth', 2.5, 'MaxHeadSize', 0.7);
            off = -0.010 * signs(k,:);
            text(p(1)+off(1), p(2)+off(2), p(3)+off(3), sprintf('T%d', k), ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', [0 0.5 0], ...
                'HorizontalAlignment', 'center');
        else
            plot3(p(1), p(2), p(3), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7]);
        end
    end

    % Determine which axis this mode acts on (1=x, 2=y, 3=z)
    axIdx = mod(panel - 1, 3) + 1;
    modeClr = axClr(axIdx, :);

    % Show net result
    Fvec = zeros(nT, 1); Fvec(activeIdx) = 1;
    wr = W * Fvec;
    netF = wr(1:3);
    netT = wr(4:6);

    if panel <= 3
        % Translation: thick force arrow from origin
        uF = netF / norm(netF);
        quiver3(0, 0, 0, uF(1)*0.18, uF(2)*0.18, uF(3)*0.18, 0, ...
            'Color', modeClr, 'LineWidth', 3.5, 'MaxHeadSize', 0.35);
    else
        % Rotation: arc + thick axis line, color-matched to axis
        uT = netT / norm(netT);
        draw_torque_arc(uT, 0.10, modeClr);
        axSpan = 0.22;
        plot3([-1 1]*uT(1)*axSpan, [-1 1]*uT(2)*axSpan, [-1 1]*uT(3)*axSpan, ...
            '-', 'Color', modeClr, 'LineWidth', 3);
    end

    hold off;
    axis equal; grid on;
    view([1, 1, 1]);
    title(posModeLabel{panel}, 'FontSize', 16, 'FontWeight', 'bold');
end

title(tl, 'Green = active thrusters | Colored arrows = body axes (X red, Y green, Z blue)', ...
    'FontSize', 13);

fprintf('Done.  Two figures generated.\n');

%% ==== LOCAL HELPERS =====================================================

function draw_cubesat(a, b, c)
    V = [-a -b -c; +a -b -c; +a +b -c; -a +b -c;
         -a -b +c; +a -b +c; +a +b +c; -a +b +c];
    E = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
    for e = 1:size(E,1)
        plot3(V(E(e,:),1), V(E(e,:),2), V(E(e,:),3), 'k-', 'LineWidth', 1.2);
        if e == 1, hold on; end
    end
    faces = {[1 2 3 4],[5 6 7 8],[1 2 6 5],[3 4 8 7],[1 4 8 5],[2 3 7 6]};
    for f = 1:6
        patch('Vertices', V, 'Faces', faces{f}, ...
            'FaceColor', [0.85 0.90 0.95], 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end
end

function draw_torque_arc(axis, radius, clr)
    % Draw a partial-circle arc around the given axis to indicate torque
    theta = linspace(0, 1.4 * pi, 50);

    % Build two perpendicular vectors to the axis
    if abs(axis(3)) > 0.9
        ref = [1; 0; 0];
    else
        ref = [0; 0; 1];
    end
    v1 = cross(axis, ref); v1 = v1 / norm(v1);
    v2 = cross(axis, v1);

    pts = radius * (cos(theta) .* v1 + sin(theta) .* v2);
    plot3(pts(1,:), pts(2,:), pts(3,:), '-', 'Color', clr, 'LineWidth', 3);

    % Arrowhead at the end
    endPt = pts(:, end);
    tang  = pts(:, end) - pts(:, end-3);
    tang  = tang / norm(tang);
    quiver3(endPt(1), endPt(2), endPt(3), tang(1)*0.015, tang(2)*0.015, tang(3)*0.015, 0, ...
        'Color', clr, 'LineWidth', 3, 'MaxHeadSize', 4);
end
