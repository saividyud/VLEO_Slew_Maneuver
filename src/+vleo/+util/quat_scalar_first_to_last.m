function qScalarLast = quat_scalar_first_to_last(qScalarFirst)
    qScalarFirst = reshape(qScalarFirst, [], 1);
    if numel(qScalarFirst) ~= 4
        error('VLEO_Slew_Maneuver:InvalidQuaternion', ...
            'Expected a 4-element quaternion.');
    end

    qScalarLast = [qScalarFirst(2:4); qScalarFirst(1)];
end
