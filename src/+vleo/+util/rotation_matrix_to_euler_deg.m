function eulerDeg = rotation_matrix_to_euler_deg(rBodyToEci)
    eulerDeg = vleo.util.quat_to_euler_deg(dcm2quat(rBodyToEci'));
end
