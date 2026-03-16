function rBodyToEci = build_tracking_frame(rSat, vSat, rVolcano, vVolcano)
    rRel = rVolcano - rSat;
    vRel = vVolcano - vSat;

    zBody = rRel / norm(rRel);
    dist = norm(rRel);
    zBodyDot = (vRel / dist) - zBody * (dot(rRel, vRel) / dist^2);

    yBody = cross(zBody, zBodyDot);
    if norm(yBody) < 1e-10
        yBody = cross(zBody, vSat);
        if norm(yBody) < 1e-10
            yBody = cross(zBody, rSat);
        end
    end
    yBody = yBody / norm(yBody);

    xBody = cross(yBody, zBody);
    xBody = xBody / norm(xBody);
    yBody = cross(zBody, xBody);
    yBody = yBody / norm(yBody);

    rBodyToEci = [xBody, yBody, zBody];
end
