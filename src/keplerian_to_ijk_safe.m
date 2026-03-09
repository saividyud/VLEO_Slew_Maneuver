function [rIJK, vIJK] = keplerian_to_ijk_safe(a, ecc, inclDeg, raanDeg, argpDeg, nuDeg, varargin)
% KEPLERIAN_TO_IJK_SAFE Converts Keplerian elements for all orbit classes.
% It maps circular and equatorial special cases onto the alternate angular
% inputs expected by MATLAB's keplerian2ijk.

    persistent projectSetupChecked
    if isempty(projectSetupChecked)
        run(fullfile(fileparts(mfilename('fullpath')), 'ensure_project_setup.m'));
        projectSetupChecked = true;
    end

    eccTolerance = 1e-10;
    inclTolerance = 1e-12;

    isCircular = abs(ecc) < eccTolerance;
    isEquatorial = abs(sind(inclDeg)) < inclTolerance;

    if isCircular && isEquatorial
        trueLonDeg = normalize_angle_deg(raanDeg + argpDeg + nuDeg);
        [rIJK, vIJK] = keplerian2ijk(a, ecc, inclDeg, 0, 0, 0, ...
            'truelon', trueLonDeg, varargin{:});
        return;
    end

    if isCircular
        argLatDeg = normalize_angle_deg(argpDeg + nuDeg);
        [rIJK, vIJK] = keplerian2ijk(a, ecc, inclDeg, raanDeg, 0, 0, ...
            'arglat', argLatDeg, varargin{:});
        return;
    end

    if isEquatorial
        lonPerDeg = normalize_angle_deg(raanDeg + argpDeg);
        [rIJK, vIJK] = keplerian2ijk(a, ecc, inclDeg, 0, 0, nuDeg, ...
            'lonper', lonPerDeg, varargin{:});
        return;
    end

    [rIJK, vIJK] = keplerian2ijk(a, ecc, inclDeg, raanDeg, argpDeg, nuDeg, varargin{:});
end

function angleDeg = normalize_angle_deg(angleDeg)
    angleDeg = mod(angleDeg, 360);
    angleDeg(angleDeg < 0) = angleDeg(angleDeg < 0) + 360;
end
