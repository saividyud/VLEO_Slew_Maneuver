function [rVolcano, vVolcano] = volcano_state_at_time(t, rVolcano0, omegaEarth)
    cosWt = cos(omegaEarth * t);
    sinWt = sin(omegaEarth * t);
    rVolcano = [cosWt * rVolcano0(1) - sinWt * rVolcano0(2); ...
                sinWt * rVolcano0(1) + cosWt * rVolcano0(2); ...
                rVolcano0(3)];
    vVolcano = cross([0; 0; omegaEarth], rVolcano);
end
