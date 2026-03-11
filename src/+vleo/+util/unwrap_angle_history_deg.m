function angleHistoryDeg = unwrap_angle_history_deg(angleHistoryDeg)
    angleHistoryDeg = rad2deg(unwrap(deg2rad(angleHistoryDeg), [], 1));
end
