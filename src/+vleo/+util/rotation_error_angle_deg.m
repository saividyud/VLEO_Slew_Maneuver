function angleErrorDeg = rotation_error_angle_deg(qCurrent, qTarget)
    rCurrent = quat2dcm(qCurrent');
    rTarget = quat2dcm(qTarget');
    rError = rCurrent * rTarget';
    traceValue = vleo.util.clamp_scalar((trace(rError) - 1) / 2, -1, 1);
    angleErrorDeg = rad2deg(acos(traceValue));
end
