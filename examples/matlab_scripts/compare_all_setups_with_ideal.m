clear; close all; clc;

%% ============================================================
%  USER SETTINGS
% ============================================================

% Ideal/reference reconstruction file.
% This file is used as the ground truth for all comparisons.
idealFile = 'C:\SCIOSPEC_Soft\ideal_reconstruction_arrays.mat';

% Folder containing the reconstructed data files from real measurements.
realFolder = 'C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\reconstruction_results_34';

% Folder where the comparison table and plots will be saved.
outputFolder = 'C:\Users\Acer\Documents\GitHub\sciopy\examples\Data\Comparison_results';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Current injection patterns to compare.
patterns = {'square', 'triangle'};

% Skip values used during the measurements.
skipValues = 0:16;

% If skip15 is missing or unreliable, uncomment this line:
% skipValues(skipValues == 15) = [];

%% ============================================================
%  LOAD IDEAL RECONSTRUCTION
% ============================================================

if ~isfile(idealFile)
    error('Ideal file not found: %s', idealFile);
end

idealStruct = load(idealFile);
idealData = extract_reconstruction_vector(idealStruct);

%% ============================================================
%  COMPARE ALL SETUPS
% ============================================================

results = table();

for p = 1:numel(patterns)

    patternName = patterns{p};

    for s = 1:numel(skipValues)

        skip = skipValues(s);

        fprintf('\nProcessing %s skip%02d...\n', patternName, skip);

        % Expected file name, for example:
        % square_reconstruction_skip00.mat
        realFile = fullfile(realFolder, ...
            sprintf('%s_reconstruction_skip%02d.mat', patternName, skip));

        if ~isfile(realFile)
            warning('Real file not found: %s', realFile);
            continue;
        end

        % Load reconstruction from the current real measurement file.
        realStruct = load(realFile);
        realData = extract_reconstruction_vector(realStruct);

        % Both reconstructions must have the same number of elements.
        if numel(idealData) ~= numel(realData)
            warning('Size mismatch for %s skip%02d. Ideal: %d, Real: %d', ...
                patternName, skip, numel(idealData), numel(realData));
            continue;
        end

        % Remove invalid values before calculating metrics.
        validMask = isfinite(idealData) & isfinite(realData);

        idealVec = idealData(validMask);
        realVec  = realData(validMask);

        if isempty(idealVec)
            warning('No valid data for %s skip%02d.', patternName, skip);
            continue;
        end

        % Calculate comparison metrics.
        metrics = calculate_metrics(realVec, idealVec);

        % Add one complete row to the results table.
        newRow = table( ...
            string(patternName), ...
            skip, ...
            string(sprintf('%s_skip%02d', patternName, skip)), ...
            metrics.RMSE, ...
            metrics.NRMSE, ...
            metrics.MAE, ...
            metrics.Pearson, ...
            metrics.CosineSimilarity, ...
            metrics.AmplitudeError, ...
            'VariableNames', { ...
                'Pattern', ...
                'Skip', ...
                'Setup', ...
                'RMSE', ...
                'NRMSE', ...
                'MAE', ...
                'Pearson', ...
                'CosineSimilarity', ...
                'AmplitudeError' ...
            } ...
        );

        results = [results; newRow];

    end
end

%% ============================================================
%  CHECK RESULTS
% ============================================================

if isempty(results)
    error('No valid results were produced. Check file names and folders.');
end

%% ============================================================
%  SORT RESULTS
% ============================================================

% Sort by NRMSE as the main error metric.
% Lower NRMSE means better agreement with the ideal reconstruction.
results = sortrows(results, 'NRMSE', 'ascend');

%% ============================================================
%  SAVE RESULTS
% ============================================================

save(fullfile(outputFolder, 'comparison_results.mat'), 'results');
writetable(results, fullfile(outputFolder, 'comparison_results.xlsx'));

disp(' ');
disp('============================================================');
disp(' BEST SETUPS');
disp('============================================================');

disp(results(:, {'Setup', 'Pattern', 'Skip', ...
    'RMSE', 'NRMSE', 'MAE', 'AmplitudeError', ...
    'Pearson', 'CosineSimilarity'}));

bestSetup = results(1,:);

fprintf('\nBest setup: %s\n', bestSetup.Setup);
fprintf('Pattern: %s\n', bestSetup.Pattern);
fprintf('Skip: %d\n', bestSetup.Skip);
fprintf('NRMSE: %.6f\n', bestSetup.NRMSE);
fprintf('Pearson: %.6f\n', bestSetup.Pearson);
fprintf('Cosine similarity: %.6f\n', bestSetup.CosineSimilarity);
fprintf('Amplitude error: %.6f\n', bestSetup.AmplitudeError);
fprintf('\nBest setup according to NRMSE: %s\n', bestSetup.Setup);
fprintf('Pattern: %s\n', bestSetup.Pattern);
fprintf('Skip: %d\n', bestSetup.Skip);
fprintf('RMSE: %.6f\n', bestSetup.RMSE);
fprintf('NRMSE: %.6f\n', bestSetup.NRMSE);
fprintf('MAE: %.6f\n', bestSetup.MAE);
fprintf('Amplitude error: %.6f\n', bestSetup.AmplitudeError);
fprintf('Pearson correlation: %.6f\n', bestSetup.Pearson);
fprintf('Cosine similarity: %.6f\n', bestSetup.CosineSimilarity);

%% ============================================================
%  PLOTS
% ============================================================

% Error metrics: lower values are better.
plot_metric(results.Setup, results.RMSE, ...
    'RMSE', ...
    'RMSE comparison', ...
    fullfile(outputFolder, 'comparison_RMSE.png'));

plot_metric(results.Setup, results.NRMSE, ...
    'NRMSE', ...
    'NRMSE comparison', ...
    fullfile(outputFolder, 'comparison_NRMSE.png'));

plot_metric(results.Setup, results.MAE, ...
    'MAE', ...
    'MAE comparison', ...
    fullfile(outputFolder, 'comparison_MAE.png'));

plot_metric(results.Setup, results.AmplitudeError, ...
    'Amplitude error', ...
    'Amplitude error comparison', ...
    fullfile(outputFolder, 'comparison_amplitude_error.png'));

% Similarity metrics: higher values are better.
plot_metric(results.Setup, results.Pearson, ...
    'Pearson correlation', ...
    'Pearson correlation with ideal reconstruction', ...
    fullfile(outputFolder, 'comparison_Pearson.png'));

plot_metric(results.Setup, results.CosineSimilarity, ...
    'Cosine similarity', ...
    'Cosine similarity with ideal reconstruction', ...
    fullfile(outputFolder, 'comparison_CosineSimilarity.png'));
%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function vec = extract_reconstruction_vector(S)
% Extract reconstruction values from a loaded .mat file.
%
% The function intentionally avoids rec_elem_centers, because this variable
% contains element coordinates, not reconstructed conductivity values.

    fieldNames = fieldnames(S);

    % First, look for known numeric reconstruction variables.
    preferredNames = {
        'rec_elem_data_ideal'
        'rec_elem_data_real'
        'rec_elem_data'
        'elem_data'
        'reconstruction'
        'reconstruction_data'
    };

    for i = 1:numel(preferredNames)

        name = preferredNames{i};

        if isfield(S, name) && isnumeric(S.(name))
            fprintf('Using numeric variable: %s\n', name);
            vec = S.(name)(:);
            return;
        end

    end

    % If the file contains an EIDORS image structure, use elem_data or node_data.
    for i = 1:numel(fieldNames)

        name = fieldNames{i};
        value = S.(name);

        if isstruct(value)

            if isfield(value, 'elem_data')
                fprintf('Using structure variable: %s.elem_data\n', name);
                vec = value.elem_data(:);
                return;
            end

            if isfield(value, 'node_data')
                fprintf('Using structure variable: %s.node_data\n', name);
                vec = value.node_data(:);
                return;
            end

        end

    end

    % rec_elem_centers is geometry information and must not be used here.
    if isfield(S, 'rec_elem_centers')
        error(['The file contains rec_elem_centers, but this is geometry, ', ...
               'not reconstruction data. Check the reconstruction variable name.']);
    end

    error('No valid reconstruction data variable found.');

end

function metrics = calculate_metrics(realVec, idealVec)
% Calculate numerical similarity metrics between real and ideal data.

    errorVec = realVec - idealVec;

    metrics.RMSE = sqrt(mean(errorVec.^2));
    metrics.NRMSE = norm(errorVec) / norm(idealVec);
    metrics.MAE = mean(abs(errorVec));

    metrics.Pearson = manual_pearson(realVec, idealVec);

    metrics.CosineSimilarity = dot(realVec, idealVec) / ...
        (norm(realVec) * norm(idealVec));

    metrics.AmplitudeError = abs(norm(realVec) - norm(idealVec)) / ...
        norm(idealVec);

end

function r = manual_pearson(x, y)
% Calculate Pearson correlation without using the Statistics Toolbox.

    x = x(:);
    y = y(:);

    x = x - mean(x);
    y = y - mean(y);

    denominator = sqrt(sum(x.^2) * sum(y.^2));

    if denominator == 0
        r = NaN;
    else
        r = sum(x .* y) / denominator;
    end

end

function y = normalize_minmax(x)
% Normalize a vector to the range 0...1.

    x = x(:);

    xmin = min(x);
    xmax = max(x);

    if xmax == xmin
        y = zeros(size(x));
    else
        y = (x - xmin) ./ (xmax - xmin);
    end

end

function plot_metric(setupNames, values, yLabelText, plotTitle, outputFile)
% Create and save a bar plot for one comparison metric.

    figure;
    bar(values);
    xticks(1:numel(values));
    xticklabels(setupNames);
    xtickangle(45);
    ylabel(yLabelText);
    title(plotTitle);
    grid on;

    saveas(gcf, outputFile);

end