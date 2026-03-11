function eulerDeg = quat_to_euler_deg(qScalarFirst)
    qScalarFirst = reshape(qScalarFirst, 1, []);
    qScalarFirst = qScalarFirst / norm(qScalarFirst);
    [yaw, pitch, roll] = quat2angle(qScalarFirst, 'ZYX');
    eulerDeg = rad2deg([roll, pitch, yaw]);
end
