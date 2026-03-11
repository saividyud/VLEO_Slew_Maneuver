function [rEci, vEci] = keplerian_to_eci_safe(a, ecc, inclDeg, raanDeg, argpDeg, nuDeg, varargin)
    eccTolerance = 1e-10;
    inclTolerance = 1e-12;

    isCircular = abs(ecc) < eccTolerance;
    isEquatorial = abs(sind(inclDeg)) < inclTolerance;

    if isCircular && isEquatorial
        trueLonDeg = normalize_angle_deg(raanDeg + argpDeg + nuDeg);
        [rEci, vEci] = keplerian2ijk(a, ecc, inclDeg, 0, 0, 0, ...
            'truelon', trueLonDeg, varargin{:});
        return;
    end

    if isCircular
        argLatDeg = normalize_angle_deg(argpDeg + nuDeg);
        [rEci, vEci] = keplerian2ijk(a, ecc, inclDeg, raanDeg, 0, 0, ...
            'arglat', argLatDeg, varargin{:});
        return;
    end

    if isEquatorial
        lonPerDeg = normalize_angle_deg(raanDeg + argpDeg);
        [rEci, vEci] = keplerian2ijk(a, ecc, inclDeg, 0, 0, nuDeg, ...
            'lonper', lonPerDeg, varargin{:});
        return;
    end

    [rEci, vEci] = keplerian2ijk(a, ecc, inclDeg, raanDeg, argpDeg, nuDeg, varargin{:});
end

function angleDeg = normalize_angle_deg(angleDeg)
    angleDeg = mod(angleDeg, 360);
    angleDeg(angleDeg < 0) = angleDeg(angleDeg < 0) + 360;
end
