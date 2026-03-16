function clampedValue = clamp_scalar(value, lowerBound, upperBound)
    clampedValue = min(max(value, lowerBound), upperBound);
end
