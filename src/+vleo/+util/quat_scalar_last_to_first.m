function qScalarFirst = quat_scalar_last_to_first(qScalarLast)
    qScalarLast = reshape(qScalarLast, [], 1);
    if numel(qScalarLast) ~= 4
        error('VLEO_Slew_Maneuver:InvalidQuaternion', ...
            'Expected a 4-element quaternion.');
    end

    qScalarFirst = [qScalarLast(4); qScalarLast(1:3)];
end
