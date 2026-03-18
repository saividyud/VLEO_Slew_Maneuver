function q_ref = Ref(t,X)
    longitudeDeg = 15;
    latitudeDeg = 0;
    refPoint0 = [6378e3 * cosd(longitudeDeg) * cosd(latitudeDeg), ...
        6378e3 * sind(longitudeDeg) * cosd(latitudeDeg), ...
        6378e3 * sind(latitudeDeg)];

    earthRate = 7.29e-5;
    refPoint = refPoint0;
    rEcefToEci = [cos(earthRate * t), -sin(earthRate * t), 0; ...
        sin(earthRate * t), cos(earthRate * t), 0; ...
        0, 0, 1];
    refPoint = (rEcefToEci * refPoint0')';
    vPoint = refPoint - X(1:3)';
    v1 = vPoint(1:3);
    v1 = v1/norm(v1);
    up = [0, 0, 1];
    v2 = cross(up, v1);
    v2 = v2 / norm(v2);
    v3 = cross(v1, v2);
    v3 = v3 / norm(v3);
    dcm = [v1', v2', v3'];
    q_ref = dcm2quat(dcm)';
end