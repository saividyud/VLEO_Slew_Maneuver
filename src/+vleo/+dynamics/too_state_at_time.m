function [rTOO, vTOO] = too_state_at_time(t, rTOO0, omegaEarth)
    cosWt = cos(omegaEarth * t);
    sinWt = sin(omegaEarth * t);
    rTOO = [cosWt * rTOO0(1) - sinWt * rTOO0(2); ...
            sinWt * rTOO0(1) + cosWt * rTOO0(2); ...
            rTOO0(3)];
    vTOO = cross([0; 0; omegaEarth], rTOO);
end
