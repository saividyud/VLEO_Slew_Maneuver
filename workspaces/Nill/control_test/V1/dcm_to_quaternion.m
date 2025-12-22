%% Helper Function: DCM to Quaternion
function q = dcm_to_quaternion(R)
    % (Function code is the same as original)
    trace_R = trace(R);
    if trace_R > 0
        s = 0.5 / sqrt(trace_R + 1.0);
        qw = 0.25 / s;
        qx = (R(3,2) - R(2,3)) * s;
        qy = (R(1,3) - R(3,1)) * s;
        qz = (R(2,1) - R(1,2)) * s;
    elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        s = 2.0 * sqrt(1.0 + R(1,1) - R(2,2) - R(3,3));
        qw = (R(3,2) - R(2,3)) / s;
        qx = 0.25 * s;
        qy = (R(1,2) + R(2,1)) / s;
        qz = (R(1,3) + R(3,1)) / s;
    elseif R(2,2) > R(3,3)
        s = 2.0 * sqrt(1.0 + R(2,2) - R(1,1) - R(3,3));
        qw = (R(1,3) - R(3,1)) / s;
        qx = (R(1,2) + R(2,1)) / s;
        qy = 0.25 * s;
        qz = (R(2,3) + R(3,2)) / s;
    else
        s = 2.0 * sqrt(1.0 + R(3,3) - R(1,1) - R(2,2));
        qw = (R(2,1) - R(1,2)) / s;
        qx = (R(1,3) + R(3,1)) / s;
        qy = (R(2,3) + R(3,2)) / s;
        qz = 0.25 * s;
    end
    q = [qx; qy; qz; qw];
    q = q / norm(q);
end