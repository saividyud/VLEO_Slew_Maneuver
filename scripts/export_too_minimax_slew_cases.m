projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup_project();

if ~exist('caseSeeds', 'var')
    caseSeeds = 1:12;
end

clearvars -except caseSeeds
clc;

validateattributes(caseSeeds, {'numeric'}, {'vector', 'real', 'finite', 'nonnegative'}, mfilename, 'caseSeeds');
caseSeeds = unique(double(caseSeeds(:)'));
if any(caseSeeds ~= floor(caseSeeds))
    error('export_control_test4_cases:InvalidCaseSeeds', ...
        'caseSeeds must contain only nonnegative integers.');
end

nCases = numel(caseSeeds);
cases = repmat(vleo.analysis.build_too_minimax_slew_case(caseSeeds(1)), nCases, 1);

fprintf('Exporting %d TOO minimax slew solver cases...\n', nCases);
for caseIdx = 1:nCases
    cases(caseIdx) = vleo.analysis.build_too_minimax_slew_case(caseSeeds(caseIdx));
    fprintf('  [%02d/%02d] %s  duration=%.2f s  angle=%.2f deg\n', ...
        caseIdx, nCases, cases(caseIdx).caseId, cases(caseIdx).maneuverDurationSec, ...
        cases(caseIdx).relativeAngleDeg);
end

generatedDir = vleo.util.generated_data_dir();
matPath = fullfile(generatedDir, 'too_minimax_slew_solver_cases.mat');
csvPath = fullfile(generatedDir, 'too_minimax_slew_solver_cases.csv');

metadata = struct();
metadata.sourceScript = 'scripts/export_too_minimax_slew_cases.m';
metadata.createdAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
metadata.caseSeeds = caseSeeds;
metadata.caseCount = nCases;

save(matPath, 'cases', 'metadata');
writetable(build_case_table(cases), csvPath);

fprintf('\nSaved MAT case bundle: %s\n', matPath);
fprintf('Saved CSV summary: %s\n', csvPath);

function caseTable = build_case_table(cases)
    nCases = numel(cases);
    caseId = strings(nCases, 1);
    scenarioSeed = zeros(nCases, 1);
    maneuverDurationSec = zeros(nCases, 1);
    relativeAngleDeg = zeros(nCases, 1);
    rollInitialDeg = zeros(nCases, 1);
    pitchInitialDeg = zeros(nCases, 1);
    yawInitialDeg = zeros(nCases, 1);
    rollTargetDeg = zeros(nCases, 1);
    pitchTargetDeg = zeros(nCases, 1);
    yawTargetDeg = zeros(nCases, 1);
    wxInitialDegPerSec = zeros(nCases, 1);
    wyInitialDegPerSec = zeros(nCases, 1);
    wzInitialDegPerSec = zeros(nCases, 1);
    wxTargetDegPerSec = zeros(nCases, 1);
    wyTargetDegPerSec = zeros(nCases, 1);
    wzTargetDegPerSec = zeros(nCases, 1);

    for caseIdx = 1:nCases
        caseId(caseIdx) = string(cases(caseIdx).caseId);
        scenarioSeed(caseIdx) = cases(caseIdx).scenarioSeed;
        maneuverDurationSec(caseIdx) = cases(caseIdx).maneuverDurationSec;
        relativeAngleDeg(caseIdx) = cases(caseIdx).relativeAngleDeg;
        rollInitialDeg(caseIdx) = cases(caseIdx).eulerInitialDeg(1);
        pitchInitialDeg(caseIdx) = cases(caseIdx).eulerInitialDeg(2);
        yawInitialDeg(caseIdx) = cases(caseIdx).eulerInitialDeg(3);
        rollTargetDeg(caseIdx) = cases(caseIdx).eulerTargetDeg(1);
        pitchTargetDeg(caseIdx) = cases(caseIdx).eulerTargetDeg(2);
        yawTargetDeg(caseIdx) = cases(caseIdx).eulerTargetDeg(3);
        wxInitialDegPerSec(caseIdx) = cases(caseIdx).omegaInitialDegPerSec(1);
        wyInitialDegPerSec(caseIdx) = cases(caseIdx).omegaInitialDegPerSec(2);
        wzInitialDegPerSec(caseIdx) = cases(caseIdx).omegaInitialDegPerSec(3);
        wxTargetDegPerSec(caseIdx) = cases(caseIdx).omegaTargetDegPerSec(1);
        wyTargetDegPerSec(caseIdx) = cases(caseIdx).omegaTargetDegPerSec(2);
        wzTargetDegPerSec(caseIdx) = cases(caseIdx).omegaTargetDegPerSec(3);
    end

    caseTable = table(caseId, scenarioSeed, maneuverDurationSec, relativeAngleDeg, ...
        rollInitialDeg, pitchInitialDeg, yawInitialDeg, ...
        rollTargetDeg, pitchTargetDeg, yawTargetDeg, ...
        wxInitialDegPerSec, wyInitialDegPerSec, wzInitialDegPerSec, ...
        wxTargetDegPerSec, wyTargetDegPerSec, wzTargetDegPerSec);
end
