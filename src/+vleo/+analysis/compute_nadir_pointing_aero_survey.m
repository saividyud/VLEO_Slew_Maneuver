function survey = compute_nadir_pointing_aero_survey(altitude_km, varargin)
    parser = inputParser();
    parser.addRequired('altitude_km', @(x) isscalar(x) && isfinite(x) && x > 0);
    parser.addParameter('nOrbits', 3, @(x) isscalar(x) && isfinite(x) && x >= 1);
    parser.addParameter('samplesPerOrbit', 120, @(x) isscalar(x) && isfinite(x) && x >= 10);
    parser.addParameter('gsi_model', 'cook', @(x) ischar(x) || isstring(x));
    parser.addParameter('thrusterLayoutMode', 'standard', @(x) islogical(x) || ischar(x) || isstring(x));
    parser.addParameter('Isp_s', 760, @(x) isscalar(x) && isfinite(x) && x > 0);
    parser.addParameter('usableFuelMass_kg', 1.5, @(x) isscalar(x) && isfinite(x) && x > 0);
    parser.addParameter('showProgress', false, @(x) islogical(x) && isscalar(x));
    parser.parse(altitude_km, varargin{:});

    c = vleo.util.constants();

    altitude_m = parser.Results.altitude_km * 1e3;
    orbitRadius_m = c.R_earth + altitude_m;
    meanMotion = sqrt(c.mu_earth / orbitRadius_m^3);
    orbitalPeriod_s = 2 * pi / meanMotion;
    totalTime_s = parser.Results.nOrbits * orbitalPeriod_s;

    nSamples = parser.Results.nOrbits * parser.Results.samplesPerOrbit + 1;
    tSec = linspace(0, totalTime_s, nSamples).';

    time0 = struct('year', 2002, 'dayOfYear', 106, 'UTseconds', 43200);
    objFile = vleo.util.geometry_asset_path('6U CubeSat.obj');
    aeroOpts = struct( ...
        'f107Average', 150, ...
        'f107Daily', 150, ...
        'magneticIndex', ones(1, 7) * 3, ...
        'anomalousOxygen', 0, ...
        'gsi_model', char(parser.Results.gsi_model), ...
        'alpha', 1, ...
        'Tw', 300, ...
        'enableShadow', 0, ...
        'enableSolar', 0 ...
    );

    attitude = struct('aoa', 180, 'aos', 0);

    forceBody_N = zeros(nSamples, 3);
    torqueBody_Nm = zeros(nSamples, 3);
    density_kg_m3 = zeros(nSamples, 1);

    for idx = 1:nSamples
        theta = meanMotion * tSec(idx);
        posEci = orbitRadius_m * [cos(theta); sin(theta); 0];

        timeNow = time0;
        timeNow.UTseconds = time0.UTseconds + tSec(idx);

        aeroResult = vleo.aero.compute_aero_forces(objFile, ...
            struct('positionECI', posEci), attitude, timeNow, aeroOpts);

        dynamicPressure = 0.5 * aeroResult.rho * aeroResult.vinf^2;
        forceBody_N(idx, :) = (dynamicPressure * aeroResult.AreaRef * aeroResult.Cf_b).';
        torqueBody_Nm(idx, :) = aeroResult.M_aero.';
        density_kg_m3(idx) = aeroResult.rho;

        if parser.Results.showProgress && mod(idx, 60) == 0
            fprintf('  %d / %d  (%.0f%%)\n', idx, nSamples, 100 * idx / nSamples);
        end
    end

    params = vleo.dynamics.default_control_test_params(true);
    layout = vleo.control.corner_thruster_layout(params, parser.Results.thrusterLayoutMode);
    requiredWrenchHistory = -[forceBody_N, torqueBody_Nm];
    thrusterForces_N = vleo.control.allocate_thruster_forces(requiredWrenchHistory, params, parser.Results.thrusterLayoutMode);
    totalThrusterForce_N = sum(thrusterForces_N, 2);
    averageTotalThrusterForce_N = mean(totalThrusterForce_N);

    g0 = 9.80665;
    massFlowRate_kg_s = averageTotalThrusterForce_N / (parser.Results.Isp_s * g0);
    if massFlowRate_kg_s > 0
        missionLife_s = parser.Results.usableFuelMass_kg / massFlowRate_kg_s;
    else
        missionLife_s = inf;
    end

    survey = struct();
    survey.altitude_km = parser.Results.altitude_km;
    survey.nOrbits = parser.Results.nOrbits;
    survey.samplesPerOrbit = parser.Results.samplesPerOrbit;
    survey.gsi_model = char(parser.Results.gsi_model);
    survey.thruster_layout_mode = layout.mode;
    survey.time0 = time0;
    survey.tSec = tSec;
    survey.timeMin = tSec / 60;
    survey.orbitalPeriod_s = orbitalPeriod_s;
    survey.totalTime_s = totalTime_s;
    survey.forceBody_N = forceBody_N;
    survey.torqueBody_Nm = torqueBody_Nm;
    survey.density_kg_m3 = density_kg_m3;
    survey.requiredWrenchHistory = requiredWrenchHistory;
    survey.thrusterForces_N = thrusterForces_N;
    survey.totalThrusterForce_N = totalThrusterForce_N;
    survey.averageTotalThrusterForce_N = averageTotalThrusterForce_N;
    survey.Isp_s = parser.Results.Isp_s;
    survey.usableFuelMass_kg = parser.Results.usableFuelMass_kg;
    survey.massFlowRate_kg_s = massFlowRate_kg_s;
    survey.massFlowRate_g_day = massFlowRate_kg_s * 86400 * 1e3;
    survey.missionLife_s = missionLife_s;
    survey.missionLife_days = missionLife_s / 86400;
    survey.missionLife_years = missionLife_s / (86400 * 365.25);
    survey.layout = layout;
end
